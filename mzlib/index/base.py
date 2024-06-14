"""Abstract base types for spectral library indices"""
import io
import pathlib
import warnings

from typing import Collection, Iterator, Optional, Tuple, Union, Any, List, NamedTuple


class IndexInitializedRecord(NamedTuple):
    """
    An instantiated index object, and whether it was initialized already or the caller
    needs to populate it.

    Attributes
    ----------
    index : :class:`~.IndexBase`
        A newly index object
    initialized : :class:`bool`
        Whether the index was populated or not
    """

    index: 'IndexBase'
    initialized: bool


class IndexRecordBase:
    """
    A base type describing index entries

    Attributes
    ----------
    number : int
        The "key" of a library entry
    index : int
        The sequential index for the entry
    offset : int
        The offset into the file for the entry. Usually in bytes but
        for some file types which use a container format that abstracts
        the IO, this might just be the index again.
    """

    __slots__ = ()

    number: int
    index: int
    offset: int
    name: str


class IndexBase(Collection):
    """
    A base type for spectral indices.

    Retrieve information about entries' identifiers and any associated
    metadata.
    """

    @classmethod
    def from_filename(cls, filename: Union[str, pathlib.Path, io.FileIO], library=None) -> IndexInitializedRecord:
        """
        Get a file path for an index file, given the library filename.

        Returns
        -------
        str or None
        """
        raise NotImplementedError()

    @classmethod
    def exists(cls, filename: Union[str, pathlib.Path, io.FileIO]) -> bool:
        """Check if an index file exists"""
        return False

    def offset_for(self, record_label) -> int:
        """Retrieve the byte offset of a spectrum identifier"""
        record = self.record_for(record_label)
        return record.offset

    def offset_for_cluster(self, record_label) -> int:
        """Retrieve the byte offset of a cluster identifier"""
        record = self.record_for_cluster(record_label)
        return record.offset

    def record_for(self, record_label: Union[int, str]) -> IndexRecordBase:
        """Retrieve a an index record for a spectrum identifier"""
        record = self.search(record_label)
        if isinstance(record, list):
            warnings.warn(
                f"Multiple records found for {record_label}, using the first")
            record = record[0]
        return record

    def record_for_cluster(self, record_label: int) -> IndexRecordBase:
        """Retrieve a an index record for a cluster identifier"""
        record = self.search_clusters(record_label)
        if isinstance(record, list):
            warnings.warn(
                f"Multiple records found for {record_label}, using the first")
            record = record[0]
        return record

    def search(self, i: Union[str, int, slice], **kwargs) -> Union[IndexRecordBase, List[IndexRecordBase]]:
        """Search for one or more spectrum records by index, slice or identifier"""
        raise NotImplementedError()

    def search_clusters(self, i: Optional[Union[int, slice]]=None, **kwargs) -> Union[IndexRecordBase, List[IndexRecordBase]]:
        """Search for one or more cluster records by index, slice or identifier"""
        raise NotImplementedError()

    def add(self, number: int, offset: int, name: str, analyte: Any, attributes=None):
        """
        Add a new entry to the spectrum index.

        Parameters
        ----------
        number : int
            A numerical identifier for this spectrum.
        offset : int
            The offset in the file to reach the spectrum (in bytes if appropriate)
        name : str,
            A text identifier for this spectrum.
        analyte : str, optional
            A text representation of the analyte for that record
        attributes : Dict[str, Any], optional
            A key-value pair collection of this record, currently not supported.
        """
        raise NotImplementedError()

    def add_cluster(self, number: int, offset: int, attributes=None):
        """
        Add a new entry to the cluster index.

        Parameters
        ----------
        number : int
            A numerical identifier for this spectrum.
        offset : int
            The offset in the file to reach the spectrum (in bytes if appropriate)
        attributes : Dict[str, Any], optional
            A key-value pair collection of this record, currently not supported.
        """
        raise NotImplementedError()

    def commit(self):
        """
        Commit any index state to disk, if this index supports persistence.

        Has no effect on index types that do not have a persistence functionality.
        """
        raise NotImplementedError()

    def iter_clusters(self) -> Iterator[IndexRecordBase]:
        """Iterate over cluster records"""
        raise NotImplementedError()

    def iter_spectra(self) -> Iterator[IndexRecordBase]:
        """Iterate over peptide records"""
        for i in range(len(self)):
            yield self[i]

    def _get_by_index(self, i: Union[int, slice]) -> Union[IndexRecordBase, List[IndexRecordBase]]:
        raise NotImplementedError()

    def __iter__(self):
        return self.iter_spectra()

    def __getitem__(self, i: Union[int, str, slice]):
        return self.search(i)

    def __len__(self):
        raise NotImplementedError()

    def __contains__(self, key) -> bool:
        try:
            hit = self.search(key)
            return True
        except (KeyError, IndexError, ValueError):
            return False

    def check_names_unique(self) -> bool:
        """
        Check that all indexed spectra have unique ``spectrum name`` parameters.

        Returns
        -------
        bool:
            Whether the spectrum names in the index are unique.
        """
        seen = set()
        for record in self:
            if record.name in seen:
                return False
            seen.add(record.name)
        return True
