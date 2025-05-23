"""A set of base types, functions and constants for implementing spectral library formats."""
import io
import csv
import enum
import logging
import os
import warnings

from typing import Any, Callable, Dict, Iterable, Optional, Union, List, Type, Iterator, TYPE_CHECKING
from pathlib import Path


from psims.controlled_vocabulary import Entity
from psims.controlled_vocabulary.controlled_vocabulary import (
    load_uo, load_unimod, load_psims)
from mzspeclib.cluster import SpectrumCluster

from mzspeclib.index import MemoryIndex, SQLIndex, IndexBase
from mzspeclib.spectrum import Spectrum
from mzspeclib.analyte import Analyte, Interpretation, InterpretationMember
from mzspeclib.attributes import Attributed, AttributedEntity, AttributeSet, AttributeManagedProperty
from mzspeclib.ontology import _VocabularyResolverMixin
from mzspeclib.const import (
    FORMAT_VERSION,
    LIBRARY_NAME,
    LIBRARY_IDENTIFIER,
    LIBRARY_VERSION,
    LIBRARY_URI,
    LIBRARY_DESCRIPTION,
    ANALYTE_MIXTURE,
    LIBRARY_SPECTRUM_INDEX,
    LIBRARY_SPECTRUM_KEY
)

from .utils import open_stream, _LineBuffer

if TYPE_CHECKING:
    from mzspeclib import SpectrumLibrary

logger = logging.getLogger(__name__.rsplit(".", 1)[0])
logger.addHandler(logging.NullHandler())

ANALYTE_MIXTURE_CURIE = ANALYTE_MIXTURE.split("|")[0]

DEFAULT_VERSION = '1.0'


class AttributeSetTypes(enum.Enum):
    """Attribute set type tags used as keys and section constants"""

    spectrum = enum.auto()
    analyte = enum.auto()
    interpretation = enum.auto()
    cluster = enum.auto()


class FormatInferenceFailure(ValueError):
    """Indicates that we failed to infer the format type for a spectral library."""


class SubclassRegisteringMetaclass(type):
    """
    A metaclass that registers instances instance types on the following attributes:
        - :attr:`file_format`
        - :attr:`format_name`

    It will also create :attr:`_file_extension_to_implementation` on the instance type
    if it does not already exist in the type's inheritance hierarchy.

    This is meant to be used to automate registration of new backends so as soon as the
    type is imported, it is registered.
    """

    def __new__(mcs, name, parents, attrs):
        new_type = type.__new__(mcs, name, parents, attrs)
        if not hasattr(new_type, "_file_extension_to_implementation"):
            new_type._file_extension_to_implementation = dict()
        if not hasattr(new_type, "_file_format_to_implementation"):
            new_type._file_format_to_implementation = dict()

        file_extension = attrs.get("file_format")
        if file_extension is not None:
            if isinstance(file_extension, list):
                for ext in file_extension:
                    new_type._file_extension_to_implementation[ext] = new_type
            else:
                new_type._file_extension_to_implementation[file_extension] = new_type

        format_name = attrs.get("format_name")
        if format_name is not None:
            new_type._file_extension_to_implementation[format_name] = new_type
            new_type._file_format_to_implementation[format_name] = new_type
        else:
            attrs['format_name'] = file_extension
            new_type._file_format_to_implementation[file_extension] = new_type
        return new_type

    def type_for_format(cls, format_or_extension: str) -> Type:
        """Get the implementing type for a specific format or file extension, if it exists"""
        return cls._file_extension_to_implementation.get(format_or_extension)


class _LibraryViewMixin:

    name = AttributeManagedProperty[str](LIBRARY_NAME)
    identifier = AttributeManagedProperty[str](LIBRARY_IDENTIFIER)
    description = AttributeManagedProperty[str](LIBRARY_DESCRIPTION)
    uri = AttributeManagedProperty[str](LIBRARY_URI)
    library_version = AttributeManagedProperty[str](LIBRARY_VERSION)

    @property
    def format_version(self):
        try:
            value = self.get_attribute(FORMAT_VERSION)
            return value
        except KeyError:
            value = DEFAULT_VERSION
            self.add_attribute(FORMAT_VERSION, value)
            return value


class SpectralLibraryBackendBase(AttributedEntity, _VocabularyResolverMixin, _LibraryViewMixin, metaclass=SubclassRegisteringMetaclass):
    """A base class for all spectral library formats."""

    file_format = None

    _file_extension_to_implementation: Dict[str,
                                            Type['SpectralLibraryBackendBase']] = {}
    _format_name_to_implementation: Dict[str,
                                         Type['SpectralLibraryBackendBase']] = {}

    index: IndexBase

    spectrum_attribute_sets: Dict[str, AttributeSet]
    analyte_attribute_sets: Dict[str, AttributeSet]
    interpretation_attribute_sets: Dict[str, AttributeSet]
    cluster_attribute_sets: Dict[str, AttributeSet]

    @classmethod
    def guess_from_filename(cls, filename: Union[str, Path, io.FileIO]) -> bool:
        """
        Guess if the file is of this type by inspecting the file's name and extension.

        Parameters
        ----------
        filename : str
            The path to the file to inspect.

        Returns
        -------
        bool:
            Whether this is an appropriate backend for that file.
        """
        if hasattr(filename, "name"):
            filename = filename.name
        if not isinstance(filename, (str, Path)):
            return False
        if filename.endswith(".gz"):
            filename = filename[:-3]
        if isinstance(cls.file_format, list):
            return any(filename.endswith(ext) for ext in cls.file_format)
        return filename.endswith(cls.file_format)

    @classmethod
    def guess_from_header(cls, filename: Union[str, Path, io.FileIO]) -> bool:
        """
        Guess if the file is of this type by inspecting the file's header section

        Parameters
        ----------
        filename : str
            The path to the file to open.

        Returns
        -------
        bool:
            Whether this is an appropriate backend for that file.
        """
        return False

    @classmethod
    def guess_implementation(cls, filename: Union[str, Path, io.FileIO], index_type=None,
                             **kwargs) -> 'SpectralLibraryBackendBase':
        """
        Guess the backend implementation to use with this file format.

        Parameters
        ----------
        filename : str
            The path to the spectral library file to open.
        index_type : type, optional
            The :class:`~.IndexBase` derived type to use for this file. If
            :const:`None` is provided, the instance will decide based upon
            :meth:`has_index_preference`.
        **kwargs
            Passed to implementation

        Returns
        -------
        SpectralLibraryBackendBase
        """
        for key, impl in cls._file_extension_to_implementation.items():
            try:
                if impl.guess_from_filename(filename):
                    return impl(filename, index_type=index_type, **kwargs)
            except TypeError:
                pass
            try:
                if impl.guess_from_header(filename):
                    return impl(filename, index_type=index_type, **kwargs)
            except (TypeError, UnicodeDecodeError):
                pass
        raise FormatInferenceFailure(f"Could not guess backend implementation for {filename}")

    def __init__(self, filename: Union[str, Path, io.FileIO]):
        self.filename = filename
        self.index = MemoryIndex()

        self.spectrum_attribute_sets = {
            "all": AttributeSet("all", [])
        }
        self.analyte_attribute_sets = {
            "all": AttributeSet("all", [])
        }
        self.interpretation_attribute_sets = {
            "all": AttributeSet("all", [])
        }
        self.cluster_attribute_sets = {
            "all": AttributeSet("all", [])
        }

        super().__init__(None)

    def _infer_lib_name(self) -> Optional[str]:
        if hasattr(self.filename, "name"):
            name = self.filename.name.replace(".gz", "").rsplit(".", 1)[0].split(os.sep)[-1]
        elif isinstance(self.filename, str):
            name = self.filename.replace(".gz", "").rsplit(".", 1)[0].split(os.sep)[-1]
        elif isinstance(self.filename, bytes):
            name = self.filename.decode("utf8").replace(".gz", "").rsplit(".", 1)[0].split(os.sep)[-1]
        else:
            name = None
        return name

    def read_header(self) -> bool:
        """
        Read just the header of the whole library

        Returns
        -------
        bool
        """
        raise NotImplementedError()

    def _new_spectrum(self) -> Spectrum:
        spec = Spectrum()
        attr_set = self.spectrum_attribute_sets.get("all")
        if attr_set:
            attr_set.apply(spec)
        return spec

    def _new_interpretation(self, id=None) -> Interpretation:
        interp = Interpretation(id)
        attr_set = self.interpretation_attribute_sets.get('all')
        if attr_set:
            attr_set.apply(interp)
        return interp

    def _new_interpretation_member(self, id=None) -> InterpretationMember:
        return InterpretationMember(id)

    def _new_analyte(self, id=None) -> Analyte:
        analyte = Analyte(id)
        attr_set = self.analyte_attribute_sets.get('all')
        if attr_set:
            attr_set.apply(analyte)
        return analyte

    def _new_cluster(self) -> SpectrumCluster:
        cluster = SpectrumCluster()
        attr_set = self.cluster_attribute_sets.get('all')
        if attr_set:
            attr_set.apply(cluster)
        return cluster

    def _analyte_interpretation_link(self, spectrum: Spectrum,
                                     interpretation: Interpretation):
        if (interpretation.has_attribute(ANALYTE_MIXTURE) and
            not interpretation.analytes):
            analyte_ids = interpretation.get_attribute(ANALYTE_MIXTURE)
            if isinstance(analyte_ids, str):
                term = self.find_term_for(ANALYTE_MIXTURE_CURIE)
                analyte_ids = term.value_type(analyte_ids)

            # TODO: Enforce this attribute is a string at the CV level
            # if isinstance(analyte_ids_term, int):
            #     analyte_ids = [analyte_ids_term]
            #     interpretation.replace_attribute(ANALYTE_MIXTURE_TERM, str(analyte_ids_term))
            # else:
            #     analyte_ids = analyte_ids_term.split(',')
            for analyte_id in analyte_ids:
                interpretation.add_analyte(spectrum.get_analyte(analyte_id))
        return interpretation

    def _default_interpretation_to_analytes(self, spectrum: Spectrum):
        for interpretation in spectrum.interpretations.values():
            if not interpretation.analytes:
                for analyte in spectrum.analytes.values():
                    interpretation.add_analyte(analyte)

    def _is_analyte_defined(self, analyte: Analyte) -> bool:
        """
        Check if the :class:`~.Analyte` has been sufficiently described
        by the attributes we've seen so far, either for a peptide or some
        other kind of molecule.

        If not, we'll do other work in :meth:`_hoist_analyte_attributes_on_rejection`
        to move information out of the :class:`~.Analyte` up to the :class:`~.Spectrum`
        and omit the analyte entirely.
        """
        peptide_attr_root = self.find_term_for("MS:1003050")
        molecular_attr_root = self.find_term_for("MS:1003033")
        for attrib in analyte:
            acc, _name = attrib.key.split("|")
            term = self.find_term_for(acc)
            if term.is_of_type(peptide_attr_root):
                return True
            if term.is_of_type(molecular_attr_root):
                return True
        return False

    def _hoist_analyte_attributes_on_rejection(self, analyte: Analyte, spectrum: Spectrum):
        """
        When an :class:`~.Analyte` is not defined well enough, move any terms
        that are ion-related up to the :class:`~.Spectrum`.
        """
        ion_selection_root = self.find_term_for("MS:1000455")
        attribute_groups_to_hoist = {}
        for attrib in analyte:
            acc, _name = attrib.key.split("|")
            term = self.find_term_for(acc)
            if attrib.group_id in attribute_groups_to_hoist:
                attribute_groups_to_hoist[attrib.group_id].append(attrib)
            if term.is_of_type(ion_selection_root):
                if attrib.group_id is None:
                    spectrum.add_attribute(attrib.key, attrib.value)
                else:
                    attribute_groups_to_hoist[attrib.group_id] = [attrib]

        for group in attribute_groups_to_hoist.values():
            spectrum.add_attribute_group(group)

    def get_spectrum(self, spectrum_number: int=None,
                     spectrum_name: str=None) -> Spectrum:
        """
        Retrieve a single spectrum from the library.

        Parameters
        ----------
        spectrum_number : int, optional
            The index of the spectrum in the library
        spectrum_name : str, optional
            The name of the spectrum in the library

        Returns
        -------
        :class:`~.Spectrum`
        """
        raise NotImplementedError()

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
        raise NotImplementedError()

    def find_spectra(self, specification, **query_keys):
        raise NotImplementedError()

    def create_index(self) -> int:
        """
        Populate the spectrum index.

        This method may produce a large amount of file I/O.

        Returns
        -------
        n_spectra: int
            The number of entries read
        """
        raise NotImplementedError()

    def __iter__(self) -> Iterator[Spectrum]:
        if self.index:
            for record in self.index:
                yield self.get_spectrum(record.number)
        else:
            return self.read()

    def __len__(self):
        return len(self.index)

    def __getitem__(self, i) -> Union[Spectrum, List[Spectrum]]:
        record = self.index[i]
        if isinstance(record, list):
            result = [self.get_spectrum(rec.number) for rec in record]
        else:
            result = self.get_spectrum(record.number)
        return result

    @classmethod
    def has_index_preference(cls, filename: Union[str, Path, io.FileIO]) -> Type[IndexBase]:
        """
        Check if this backend prefers a particular index for this file.

        The base implementation checks to see if there is a SQL index
        for the filename provided, and if so, prefers :class:`~.SQLIndex`.
        Otherwise, prefers :class:`~.MemoryIndex`.

        Parameters
        ----------
        filename: str
            The name of the file to open.

        Returns
        -------
        index_type: type
            Returns a :class:`~.IndexBase` derived type which this backend
            would prefer to use.
        """
        try:
            if SQLIndex.exists(filename):
                return SQLIndex
            return MemoryIndex
        except Exception:
            return MemoryIndex

    def read(self) -> Iterator[Union[Spectrum, SpectrumCluster]]:
        """
        Create an sequential iterator over the spectrum library.

        Yields
        ------
        entry :  Union[:class:`~.Spectrum`, :class:`~.SpectrumCluster`]
        """
        raise NotImplementedError()

    def _add_attribute_set(self, attribute_set: AttributeSet,
                           attribute_set_type: AttributeSetTypes):
        if attribute_set_type == AttributeSetTypes.spectrum:
            self.spectrum_attribute_sets[attribute_set.name] = attribute_set
        elif attribute_set_type == AttributeSetTypes.analyte:
            self.analyte_attribute_sets[attribute_set.name] = attribute_set
        elif attribute_set_type == AttributeSetTypes.interpretation:
            self.interpretation_attribute_sets[attribute_set.name] = attribute_set
        elif attribute_set_type == AttributeSetTypes.cluster:
            self.cluster_attribute_sets[attribute_set.name] = attribute_set
        else:
            raise ValueError(f"Could not map {attribute_set_type}")

    def summarize_parsing_errors(self) -> Dict:
        """Retrieve a free-form description of parsing errors"""
        return {}


guess_implementation = SpectralLibraryBackendBase.guess_implementation


class _PlainTextSpectralLibraryBackendBase(SpectralLibraryBackendBase):

    def __init__(self, filename: Union[str, Path, io.FileIO], index_type=None,
                 read_metadata: bool=True, create_index: bool=True):
        if index_type is None and create_index:
            index_type = self.has_index_preference(filename)

        super(_PlainTextSpectralLibraryBackendBase, self).__init__(filename)

        if index_type is not None:
            self.index, was_initialized = index_type.from_filename(filename)
            if not was_initialized and create_index:
                self.create_index()
        if read_metadata:
            self.read_header()

    def _coerce_handle(self, filename_or_stream):
        if hasattr(filename_or_stream, 'read'):
            self.handle = filename_or_stream
        else:
            self.handle = open_stream(filename_or_stream, 'rt')

    def _buffer_from_stream(self, stream: io.IOBase) -> List:
        """
        Collect data from the readable stream until
        a complete spectrum entry has been observed.

        Parameters
        ----------
        stream: file-like
            Theinput file stream to read from.

        Returns
        -------
        line_buffer: List[str]
            A list of lines read from the input stream.
        """
        raise NotImplementedError()

    def read(self) -> Iterator[Spectrum]:
        with open_stream(self.filename, 'rb') as stream:
            i = 0
            match, offset = self._parse_header_from_stream(stream)
            if not match:
                raise ValueError("Could not locate valid header")
            else:
                stream.seek(offset)
            buffering_stream = _LineBuffer(stream, encoding="utf8")
            while True:
                # Will clip the first line of the next spectrum. Needs work
                buffer = self._buffer_from_stream(buffering_stream)

                # If the buffer is only a single line, then we must have reached
                # the end, so we're done. We're done because the buffering stream
                # will contain exactly one line (the buffered line)
                if len(buffer) <= 1:
                    break
                buffering_stream.push_line()
                yield self._parse(buffer, i)
                i += 1

    def _get_lines_for(self, offset: int) -> List[str]:
        filename = self.filename
        is_file_like_object = isinstance(filename, io.IOBase)

        infile = open_stream(filename, 'r')

        infile.seek(offset)
        spectrum_buffer = self._buffer_from_stream(infile)
        #### We will end up here if this is the last spectrum in the file
        if not is_file_like_object:
            infile.close()
        else:
            infile.detach()
        return spectrum_buffer

    def _parse(self, buffer: Iterable, spectrum_index: int=None):
        raise NotImplementedError()

    def search(self, specification, **query_keys) -> List[Spectrum]:
        records = self.index.search(specification, **query_keys)
        if not isinstance(records, list):
            records = [records]
        spectra = []
        for record in records:
            buffer = self._get_lines_for(record.offset)
            spectrum = self._parse(buffer, record.number)
            spectra.append(spectrum)
        return spectra


class _CSVSpectralLibraryBackendBase(SpectralLibraryBackendBase):
    _delimiter: str
    _header: List[str]

    _required_columns: List[str] = None

    @classmethod
    def guess_from_header(cls, filename) -> bool:
        with open_stream(filename, 'rt') as fh:
            line = fh.readline()
            tokens = line.split('\t')
            if len(tokens) > 3:
                warnings.warn(
                    "This file looks like a TSV file, but it's not possible to say definitively what"
                    " type from this alone. Please explicitly specify the format.")
                return False
        return False

    def __init__(self, filename: Union[str, Path, io.FileIO], index_type=None, delimiter='\t',
                 read_metadata: bool=True, create_index: bool = True, ** kwargs):
        if index_type is None:
            index_type = self.has_index_preference(filename)
        self._delimiter = delimiter
        self._headers = None
        super().__init__(filename)
        self.filename = filename

        self.read_header()
        self.index, was_initialized = index_type.from_filename(filename)
        if not was_initialized and create_index:
            self.create_index()

    def read_header(self) -> bool:
        self._read_header_line()
        return True

    def _read_header_line(self):
        headers = None
        with open_stream(self.filename) as stream:
            reader = csv.reader(stream, delimiter=self._delimiter)
            headers = next(reader)
            stream.seek(0)
        self._headers = headers
        if headers and self._required_columns:
            missing_required = set(self._required_columns) - set(headers)
            if missing_required:
                raise TypeError(
                    f"{self.format_name} requires column{'s' if len(missing_required) > 1 else ''} "
                    f"{', '.join(missing_required)}, but were not found."
                )

    def _open_reader(self, stream: io.TextIOBase) -> Union[Iterator[Dict[str, Any]], csv.DictReader]:
        return csv.DictReader(stream, fieldnames=self._headers, delimiter='\t')

    def get_spectrum(self, spectrum_number: int = None, spectrum_name: str = None) -> Spectrum:
        # keep the two branches separate for the possibility that this is not possible with all
        # index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError("Provide only one of spectrum_number or spectrum_name")
            offset = self.index.offset_for(spectrum_number)
        elif spectrum_name is not None:
            offset = self.index.offset_for(spectrum_name)
        else:
            raise ValueError("Must provide either spectrum_number or spectrum_name argument")
        buffer = self._get_lines_for(offset)
        spectrum = self._parse_from_buffer(buffer, spectrum_number)
        return spectrum

    def _batch_rows(self, iterator: Iterator[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
        """
        Gather successive rows by some shared key into a batch to be parsed into a single
        :class:`~.Spectrum` instance.

        This assumes there are no interleaved rows across spectra.

        Parameters
        ----------
        iterator : Iterator[Dict[str, Any]]
            An iterator over rows of the underlying CSV as dictionaries

        Yields
        ------
        batch : List[Dict[str, Any]]
            The next batch of rows corresponding to a single spectrum
        """
        raise NotImplementedError()

    def _get_lines_for(self, offset: int) -> List[Dict[str, Any]]:
        with open_stream(self.filename, 'r') as infile:
            infile.seek(offset)
            reader = self._open_reader(infile)
            spectrum_buffer = next(self._batch_rows(reader))
            #### We will end up here if this is the last spectrum in the file
        return spectrum_buffer

    def read(self) -> Iterator[Spectrum]:
        with open_stream(self.filename, 'rt') as stream:
            i = 0
            reader = self._open_reader(stream)
            if self._headers:
                # Skip the header line if we've already parsed them
                _ = next(reader)
            buffering_reader = self._batch_rows(reader)
            for i, buffer in enumerate(buffering_reader):
                yield self._parse_from_buffer(buffer, i)


class SpectralLibraryWriterBase(_VocabularyResolverMixin, metaclass=SubclassRegisteringMetaclass):
    """
    A base type for spectral library writers.

    This type implements the context manager protocol, controlling the closing of the
    enclosed IO stream.

    Attributes
    ----------
    filename : str, :class:`pathlib.Path`, or :class:`io.IOBase`
    """

    _already_started_writing: bool = False

    def __init__(self, filename, **kwargs):
        self.filename = filename
        super().__init__(**kwargs)

    def _filter_attributes(self, attributes: Attributed,
                           filter_fn: Callable) -> Iterable:
        if isinstance(attributes, AttributedEntity):
            attributes = attributes.attributes
        for attrib in attributes:
            if filter_fn(attrib):
                yield attrib

    def _not_analyte_mixture_term(self, attrib):
        if attrib:
            key = attrib[0]
            if key == ANALYTE_MIXTURE:
                return False
        return True

    def _not_entry_index(self, attrib):
        if attrib:
            key = attrib[0]
            if key == LIBRARY_SPECTRUM_INDEX:
                return False
        return True

    def _not_entry_key_or_index(self, attrib):
        if attrib:
            key = attrib[0]
            if key in (LIBRARY_SPECTRUM_INDEX, LIBRARY_SPECTRUM_KEY):
                return False
        return True

    def _coerce_handle(self, filename_or_stream):
        if hasattr(filename_or_stream, 'write'):
            self.handle = filename_or_stream
        else:
            self.handle = open(filename_or_stream, 'wt')

    def write_library(self, library: Union[SpectralLibraryBackendBase, "SpectrumLibrary"]):
        """
        Write out the entire library.

        Parameters
        ----------
        library : :class:`SpectralLibraryBackendBase` or :class:`SpectrumLibrary`
            The library to write out.

        Raises
        ------
        :class:`ValueError` :
            If the writer has already started writing one library, an error will be
            raised.
        """
        if self._already_started_writing:
            raise ValueError("Cannot write a new library, already started writing")
        self._already_started_writing = True

        self.write_header(library)
        n = len(library)
        step = max(min(n // 100, 5000), 1)
        ident = ''
        i = 0
        for i, entry in enumerate(library.read()):
            if i % step == 0 and i:
                if isinstance(entry, SpectrumCluster):
                    tag = "cluster "
                else:
                    tag = ""
                try:
                    ident = f"{tag}{entry.key}:{entry.name}"
                except Exception:
                    ident = f"{tag}{entry.key}"
                logger.info(f"Wrote {ident} {i}/{n} ({i / n * 100.0:0.2f}%)")
            if isinstance(entry, Spectrum):
                self.write_spectrum(entry)
            elif isinstance(entry, SpectrumCluster):
                self.write_cluster(entry)
            else:
                raise TypeError(f"Don't know how to save {entry.__class__}")

        i = n
        logger.info(f"Wrote {n} spectra")

    def write_spectrum(self, spectrum: Spectrum):
        """
        Write out a :class:`~.Spectrum` and all of its
        components.

        Parameters
        ----------
        spectrum : :class:`~.Spectrum`
            The spectrum to write.
        """
        raise NotImplementedError()

    def write_cluster(self, cluster: SpectrumCluster):
        """
        Write out a :class:`~.SpectrumCluster` and all of its
        components.

        Parameters
        ----------
        cluster : :class:`~.SpectrumCluster`
            The spectrum cluster to write.
        """
        raise NotImplementedError()

    def __enter__(self) -> 'SpectralLibraryWriterBase':
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        """
        Close the library writer, performing any necessary finalization.

        This is called automatically when :meth:`__exit__` is called.
        """
        pass


class LibraryIterator(AttributedEntity, _LibraryViewMixin, Iterator[Spectrum]):
    """An iterator wrapper for a library source that doesn't permit random access"""

    backend: SpectralLibraryBackendBase
    attributes: Attributed
    iter: Iterator[Spectrum]
    _buffer: Optional[Spectrum]

    def __init__(self, backend: SpectralLibraryBackendBase) -> None:
        self.backend = backend
        self.attributes = backend
        self.iter = backend.read()
        try:
            self._buffer = next(self.iter)
        except StopIteration:
            self._buffer = None

    def __iter__(self):
        return self

    def __next__(self) -> Spectrum:
        if self._buffer is not None:
            result = self._buffer
            self._buffer = None
            return result
        return next(self.iter)
