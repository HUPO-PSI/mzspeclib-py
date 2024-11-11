#!/usr/bin/env python3
import os
import pathlib

from typing import Any, Iterator, Optional, Set, Type, List, Union
from mzspeclib.attributes import Attribute, AttributeManagedProperty, AttributeManager
from mzspeclib.backends.base import LIBRARY_DESCRIPTION, LIBRARY_NAME, LIBRARY_URI, LIBRARY_VERSION
from mzspeclib.cluster import SpectrumCluster

from mzspeclib.spectrum_library_index import SpectrumLibraryIndex
from mzspeclib.spectrum import Spectrum
from mzspeclib.index import IndexBase
from mzspeclib.backends import guess_implementation, SpectralLibraryBackendBase, SpectralLibraryWriterBase, FormatInferenceFailure


class SpectrumLibrary:
    """
    Read, write, and search through a spectrum library.

    This type will attempt to infer the correct spectrum library reader from the input
    file, but may need to be explicitly prompted if there are ambiguities.

    Attributes
    ----------
    identifier: str
        A unique identifier string assigned to this library by a spectral library host or
        provider.
    filename: str or Path
        A location on the local file system where the spectral library is stored
    format : str
        The name of the format for the current encoding of the library.
    backend: :class:`~.SpectralLibraryBackendBase`
        The implementation used to parse the file

    Parameters
    ----------
    identifier : str, optional
        A universal identifier for a hosted spectral library to fetch.
    filename : str, os.PathLike, or io.IOBase, optional
        A path-like or file-like object that holds a spectral library to read.
    format : string
        Name of the format for the current encoding of the library.
    index_type : Type[:class:`~.mzlib.index.base.IndexBase`]
        The type of index to preferentially construct.
    create_index : bool
        Whether to construct an index over the library if one does not exist
        already. This limits the library to sequential iteration but does not
        incur an expensive end-to-end parsing step upon opening.
    """

    backend: SpectralLibraryBackendBase
    filename: str
    identifier: str
    format: str
    index_type: Type[IndexBase]

    name = AttributeManagedProperty[str](LIBRARY_NAME)
    description = AttributeManagedProperty[str](LIBRARY_DESCRIPTION)
    uri = AttributeManagedProperty[str](LIBRARY_URI)
    library_version = AttributeManagedProperty[str](LIBRARY_VERSION)

    _identifier = None
    _filename = None
    _create_index: bool

    def __repr__(self):
        components = []

        if self.identifier:
            components.append(f"identifier={self.identifier!r}")
        if self.filename:
            components.append(f"filename={self.filename!r}")
        if self.backend is not None:
            components.append(f"backend={self.backend.__class__.__name__}")

        return f"{self.__class__.__name__}({', '.join(components)})"

    def __init__(
        self,
        identifier: Optional[str] = None,
        filename: Optional[str] = None,
        format: Optional[Type[SpectralLibraryBackendBase]]=None,
        index_type: Optional[Type[IndexBase]]=None,
        create_index: bool=True,
    ):
        self.backend = None
        self._create_index = create_index
        self._format = format
        self.identifier = identifier
        self.index_type = index_type
        self.filename = filename

    def _init_from_identifier(self):
        if os.path.exists(self.identifier) and self._filename is None and not self._backend_initialized():
            # This in turn initializes the backend
            self.filename = self.identifier

    def _init_from_filename(self, index_type: Optional[Type[IndexBase]]=None):
        if index_type is None:
            index_type = self.index_type
        if self.format is None:
            self.backend = guess_implementation(self.filename, index_type)
            self._format = self.backend.format_name
        else:
            if callable(self.format):
                backend_type = self.format
            else:
                backend_type = SpectralLibraryBackendBase.type_for_format(self.format)
            if backend_type is None:
                raise ValueError(
                    f"Could not find an implementation for {self.format}")
            self.backend = backend_type(
                self.filename, index_type=index_type, create_index=self._create_index)
            self._format = self.backend.format_name
        self._identifier = self.backend.identifier

    def _backend_initialized(self):
        return self.backend is not None

    def _requires_backend(self):
        if not self._backend_initialized():
            raise ValueError(
                "Cannot read library data, library parser not yet initialized")

    @classmethod
    def from_backend(cls, backend: SpectralLibraryBackendBase) -> "SpectrumLibrary":
        """
        Wrap a format-specific backend in a :class:`SpectrumLibrary`.

        This is useful because the respective backends do not support all
        operations.
        """
        self = cls()
        self.backend = backend
        self._format = backend.__class__
        self.filename = self.backend.filename
        self.index_type = self.backend.index.__class__
        return self

    @property
    def spectrum_attribute_sets(self):
        """The spectrum attribute sets of the spectral library"""
        self._requires_backend()
        return self.backend.spectrum_attribute_sets

    @property
    def analyte_attribute_sets(self):
        """The analyte attribute sets of the spectral library"""
        self._requires_backend()
        return self.backend.analyte_attribute_sets

    @property
    def interpretation_attribute_sets(self):
        """The interpretation attribute sets of the spectral library"""
        self._requires_backend()
        return self.backend.interpretation_attribute_sets

    @property
    def cluster_attribute_sets(self):
        """The spectrum cluster attribute sets of the spectral library"""
        self._requires_backend()
        return self.backend.cluster_attribute_sets

    #### Define getter/setter for attribute identifier
    @property
    def identifier(self) -> Optional[str]:
        """
        The spectrum library's identifier, such an accession number provided by a host or repository.

        This attribute may not exist for libraries loaded from files on disk. It relies on the presence
        of the ``MS:1003187|library identifier`` attribute in the library header.
        """
        if self._identifier is None:
            if self._backend_initialized():
                return self.backend.identifier
        return self._identifier

    @identifier.setter
    def identifier(self, identifier: Optional[str]):
        self._identifier = identifier
        if self.backend is not None:
            self.backend.identifier = identifier

    #### Define getter/setter for attribute filename
    @property
    def filename(self) -> Optional[str]:
        """
        The path to the library file on disk.

        This may also be a file-like object.
        """
        return self._filename

    @filename.setter
    def filename(self, filename):
        self._filename = filename
        if filename is not None:
            self._init_from_filename()

    #### Define getter/setter for attribute format
    @property
    def format(self):
        return self._format

    @property
    def index(self) -> IndexBase:
        if self._backend_initialized():
            return self.backend.index
        return None

    @property
    def attributes(self) -> AttributeManager:
        """The library level attributes"""
        if self._backend_initialized():
            return self.backend.attributes
        return None

    def read_header(self) -> bool:
        """
        Read just the header of the whole library

        Returns
        -------
        bool:
            Whether the operation was successful
        """
        self._requires_backend()
        return self.backend.read_header()

    def read(self):
        """
        Create a sequential iterator over the spectrum library entries.

        Yields
        ------
        Spectrum or SpectrumCluster
        """
        self._requires_backend()
        return self.backend.read()

    def write(self, destination, format: str=None, **kwargs):
        """
        Write the library to disk.

        Parameters
        ----------
        destination : str, os.PathLike, or io.IOBase
            The path or stream to write the library to.
        format : str, Type, or Callable
            The name of the format or a callable object that returns
            a :class:`~.SpectrumLibraryWriterBase`.
        **kwargs
            Passed to implementation.
        """
        filename = destination
        if not isinstance(filename, (str, pathlib.Path, os.PathLike)):
            filename = getattr(destination, "name", None)

        if format is None and filename is not None:
            basename = os.path.basename(filename)
            tokens = basename.rsplit(".", 2)
            if len(tokens) == 3:
                writer_type = SpectralLibraryWriterBase.type_for_format('.'.join(tokens[1:]))
            else:
                raise ValueError(f"Could not guess format from file name {filename}")
        else:
            writer_type = SpectralLibraryWriterBase.type_for_format(format)
        if writer_type is None:
            raise ValueError(
                f"Could not find a format writer from file name {filename} or format {format}")
        writer = writer_type(destination, **kwargs)
        if self._backend_initialized():
            with writer:
                writer.write_library(self.backend)
        else:
            print("Library not initialized")
            writer.close()

    def get_spectrum(self, spectrum_number: int=None, spectrum_name: str=None) -> Spectrum:
        """
        Retrieve a single spectrum from the library.

        Parameters
        ----------
        spectrum_number : int, optional
            The index of the specturm in the library
        spectrum_name : str, optional
            The name of the spectrum in the library

        Returns
        -------
        :class:`~.Spectrum`
        """
        self._requires_backend()
        return self.backend.get_spectrum(spectrum_number, spectrum_name)

    def get_cluster(self, cluster_number: int) -> SpectrumCluster:
        """
        Retrieve a single spectrum cluster from the library.

        Parameters
        ----------
        cluster_number : int, optional
            The index of the cluster in the library

        Returns
        -------
        :class:`~.SpectrumCluster`
        """
        self._requires_backend()
        return self.backend.get_cluster(cluster_number)

    def find_spectra(self, specification, **query_keys) -> List[Spectrum]:
        """
        Return a list of spectra given query constraints

        Returns
        -------
        List[Spectrum]
        """
        self._requires_backend()
        return self.backend.find_spectra(specification, **query_keys)

    def __getitem__(self, i: int) -> Spectrum:
        self._requires_backend()
        return self.backend[i]

    def __len__(self):
        if self._backend_initialized():
            return len(self.backend)
        return 0

    def __iter__(self) -> Iterator[Spectrum]:
        if self._backend_initialized():
            return iter(self.backend)
        return iter([])

    def add_attribute(self, key, value, group_identifier=None):
        """
        Add an attribute to the library level attributes store.

        Parameters
        ----------
        key : str
            The name of the attribute to add
        value : object
            The value of the attribute to add
        group_identifier : str, optional
            The attribute group identifier to use, if any. If not provided,
            no group is assumed.
        """
        self._requires_backend()
        return self.backend.add_attribute(key, value, group_identifier=group_identifier)

    def get_attribute(
        self, key: str, group_identifier: Optional[str] = None, raw: bool = False
    ) -> Union[Any, List[Any], Attribute, List[Attribute]]:
        """
        Get the value or values associated with a given
        attribute key.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.
        raw : bool
            Whether to return the :class:`Attribute` object or unwrap the value

        Returns
        -------
        attribute_value: object or list[object]
            Returns single or multiple values for the requested attribute.
        """
        self._requires_backend()
        return self.backend.get_attribute(key, group_identifier=group_identifier, raw=raw)

    def remove_attribute(self, key, group_identifier=None):
        """
        Remove the value or values associated with a given
        attribute key from the library level attribute store.

        This rebuilds the entire store, which may be expensive.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        """
        self._requires_backend()
        return self.backend.remove_attribute(key, group_identifier=group_identifier)

    def has_attribute(self, key):
        """
        Test for the presence of a given attribute in the library
        level store.

        Parameters
        ----------
        key : str
            The attribute to test for

        Returns
        -------
        bool
        """
        self._requires_backend()
        return self.backend.has_attribute(key)

    def summarize_parsing_errors(self):
        """Retrieve a free-form description of parsing errors"""
        return self.backend.summarize_parsing_errors()

    @staticmethod
    def supported_file_extensions() -> Set[str]:
        """
        Get the set of file extensions that are currently recognized as spectral
        library formats
        """
        return set(SpectralLibraryBackendBase._file_extension_to_implementation)


#### Example using this class
def example():

    #### Create a new RTXFeedback object
    spectrum_library = SpectrumLibrary()
    spectrum_library.filename = "../refData/sigmaups1_consensus_final_true_lib.msp"
    spectrum_library.read_header()

    spectrum_buffer = spectrum_library.get_spectrum(spectrum_index_number=2000)
    spectrum = Spectrum()
    spectrum.parse(spectrum_buffer)
    buffer = spectrum.write(format="text")
    print(buffer)
    print()

    return()


#### If this class is run from the command line, perform a short little test to see if it is working correctly
def main():

    #### Run an example
    example()
    return()


if __name__ == "__main__": main()

