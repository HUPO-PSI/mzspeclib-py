from typing import Iterator, List, Optional, Union

from mzspeclib.cluster import SpectrumCluster
from mzspeclib.index.memory import MemoryIndex, IndexRecord, ClusterIndexRecord
from mzspeclib.spectrum import Spectrum

from .base import SpectralLibraryBackendBase


class InMemorySpectrumLibrary(SpectralLibraryBackendBase):
    """
    An in-memory spectrum library backend that holds all its entries
    in memory.

    This can be used when generating a library in-silico.
    """

    spectra: List[Spectrum]
    clusters: List[SpectrumCluster]
    index: MemoryIndex

    @classmethod
    def guess_from_filename(cls, filename) -> bool:
        return filename == ":memory:"

    def __init__(self, spectra: Optional[List[Spectrum]]=None, clusters: Optional[List[SpectrumCluster]]=None, create_index: bool = True, **kwargs):
        super().__init__(":memory:")
        self.spectra = spectra or []
        self.clusters = clusters or []

        if create_index:
            self.create_index()

    def read_header(self) -> bool:
        return False

    def read(self) -> Iterator[Union[Spectrum, SpectrumCluster]]:
        yield from self.clusters
        yield from self.spectra

    def get_spectrum(self, spectrum_number: Optional[int] = None, spectrum_name: Optional[str] = None) -> Spectrum:
        # keep the two branches separate for the possibility that this is not
        # possible with all index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError("Provide only one of spectrum_number or spectrum_name")
            index_record = self.index.record_for(spectrum_number)
            offset = index_record.offset
        elif spectrum_name is not None:
            index_record = self.index.record_for(spectrum_name)
            offset = index_record.offset
            spectrum_number = index_record.number
        else:
            raise ValueError("Must provide either spectrum_number or spectrum_name argument")
        spectrum = self.spectra[offset]
        return spectrum

    def get_cluster(self, cluster_number: int) -> SpectrumCluster:
        offset = self.index.offset_for_cluster(cluster_number)
        cluster = self.clusters[offset]
        return cluster

    def create_index(self) -> int:
        index = MemoryIndex()
        for i, spec in enumerate(self.spectra):
            index.add(IndexRecord(spec.key, i, spec.name, None, i, None))
        for i, clus in enumerate(self.clusters):
            index.add_cluster(
                ClusterIndexRecord(
                    clus.key, i,None
                )
            )
        index.commit()
        self.index = index
        return len(index)

    def append(self, spectrum: Spectrum):
        """Append a :class:`~.Spectrum` to the in-memory buffer"""
        i = len(self.spectra)
        spectrum.index = i
        self.spectra.append(spectrum)
        if self.index:
            self.index.add(IndexRecord(spectrum.key, i, spectrum.name, None, i, None))
            self.index.commit()

    def append_cluster(self, cluster: SpectrumCluster):
        """Append a :class:`~.SpectrumCluster` to the in-memory buffer"""
        i = len(self.clusters)
        self.clusters.append(cluster)
        if self.index:
            self.index.add_cluster(ClusterIndexRecord(cluster.key, i, None))