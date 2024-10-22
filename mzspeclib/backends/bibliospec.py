"""
Read BiblioSpec 2 SQLite3 spectral libraries.

For more information, see `BiblioSpec's website <https://skyline.ms/project/home/software/BiblioSpec/begin.view>`_.
"""
import sqlite3
import zlib
from dataclasses import dataclass

from typing import Iterator, List, Mapping, Tuple, Iterable, Type

import numpy as np

from pyteomics import proforma

from mzspeclib import annotation
from mzspeclib.analyte import FIRST_ANALYTE_KEY, FIRST_INTERPRETATION_KEY, Analyte
from mzspeclib.spectrum import Spectrum, SPECTRUM_NAME, CHARGE_STATE
from mzspeclib.attributes import AttributeManager, Attributed
from mzspeclib.const import (
    CHARGE_STATE,
    LIBRARY_IDENTIFIER,
    PROFORMA_ION,
    LIBRARY_NAME,
    SCAN_NUMBER,
    SOURCE_FILE,
    STRIPPED_PEPTIDE_SEQ,
    RETENTION_TIME,
    SW_VERSION,
    THEORETICAL_MASS,
    NUMBER_OF_REPLICATE_SPECTRA_AVAILABLE,
    NUMBER_OF_REPLICATE_SPECTRA_USED,
    LIBRARY_CREATION_SW,
)

from mzspeclib.backends.base import SpectralLibraryBackendBase, FORMAT_VERSION, DEFAULT_VERSION

from mzspeclib.index.base import IndexBase


class BibliospecBase:
    """Minimal base class for interacting with BiblioSpec 2 libraries"""

    connection: sqlite3.Connection

    def _correct_modifications_in_sequence(self, row: Mapping) -> proforma.ProForma:
        """
        Correct the modifications in BiblioSpec's modified peptide sequence.

        Bibliospec only stores modifications as delta masses.
        """
        mods = self.connection.execute("SELECT * FROM Modifications WHERE RefSpectraID = ?", (row['id'], )).fetchall()
        peptide = proforma.ProForma.parse(row["peptideModSeq"])
        for mod in mods:
            position = mod['position'] - 1
            mass = mod['mass']
            peptide[position][1][0] = proforma.MassModification(mass)
        return peptide




@dataclass
class BlibIndexRecord:
    """An index record read from a BiblioSpec 2 spectral library"""

    number: int
    precursor_mz: float
    precursor_charge: int
    peptide: str


class BlibIndex(BibliospecBase, IndexBase):
    """An index implementation backed by the BiblioSpec 2 SQLite3 database"""

    connection: sqlite3.Connection

    def __init__(self, connection):
        self.connection = connection

    def __getitem__(self, i):
        if isinstance(i, int):
            return self.search(i + 1)
        elif isinstance(i, slice):
            return [self.search(j + 1) for j in range(i.start or 0, i.stop or len(self), i.step or 1)]
        else:
            raise TypeError(f"Cannot index {self.__class__.__name__} with {i}")

    def _record_from(self, row: Mapping) -> BlibIndexRecord:
        peptide_sequence = str(self._correct_modifications_in_sequence(row))
        return BlibIndexRecord(row['id'], row['precursorMZ'], row['precursorCharge'], peptide_sequence)

    def search(self, i):
        if isinstance(i, int):
            info = self.connection.execute("SELECT * FROM RefSpectra WHERE id = ?", (i, )).fetchone()
            return self._record_from(info)
        elif isinstance(i, str):
            raise NotImplementedError()

    def __iter__(self):
        return map(self._record_from, self.connection.execute("SELECT * FROM RefSpectra ORDER BY id").fetchall())

    def __len__(self):
        return self.connection.execute("SELECT count(*) FROM RefSpectra;").fetchone()[0]


class BibliospecSpectralLibrary(BibliospecBase, SpectralLibraryBackendBase):
    """Read BiblioSpec 2 SQLite3 spectral library files."""

    connection: sqlite3.Connection

    file_format = "blib"
    format_name = "bibliospec"

    @classmethod
    def has_index_preference(cls, filename) -> Type[IndexBase]:
        return BlibIndex

    def __init__(self, filename, **kwargs):
        super().__init__(filename)
        self.connection = sqlite3.connect(filename)
        self.connection.row_factory = sqlite3.Row
        self.index = BlibIndex(self.connection)
        self.read_header()

    def read_header(self) -> bool:
        attribs = AttributeManager()
        attribs.add_attribute(FORMAT_VERSION, DEFAULT_VERSION)
        attribs.add_attribute(LIBRARY_CREATION_SW, "Bibliospec")

        info = self.connection.execute("SELECT * FROM LibInfo;").fetchone()
        library_id = info['libLSID']
        _, pfx_name = library_id.split("bibliospec:")
        _, name = pfx_name.split(":", 1)
        attribs.add_attribute(LIBRARY_NAME, name)
        attribs.add_attribute(LIBRARY_IDENTIFIER, library_id)
        attribs.add_attribute(SW_VERSION, f"{info['majorVersion']}.{info['minorVersion']}")
        self.attributes = attribs
        return True

    def _populate_analyte(self, analyte: Analyte, row: Mapping):
        """
        Fill an analyte with details describing a peptide sequence and inferring
        from context its traits based upon the assumptions Bibliospec makes.

        Bibliospec only stores modifications as delta masses.
        """
        peptide = self._correct_modifications_in_sequence(row)
        peptide.charge_state = row['precursorCharge']
        analyte.add_attribute(PROFORMA_ION, str(peptide))
        analyte.add_attribute(THEORETICAL_MASS, peptide.mass)
        analyte.add_attribute(STRIPPED_PEPTIDE_SEQ, row['peptideSeq'])

    def get_spectrum(self, spectrum_number: int = None, spectrum_name: str = None):
        """
        Read a spectrum from the spectrum library.

        Bibliospec does not support alternative labeling of spectra with a
        plain text name so looking up by `spectrum_name` is not supported.
        """
        if spectrum_number is None:
            raise ValueError("Only spectrum number queries are supported. spectrum_number must have an integer value")
        try:
            spectrum_number = int(spectrum_number)
        except (ValueError, TypeError):
            raise ValueError(f"spectrum_number must have an integer value, received {spectrum_number!r}") from None
        info = self.connection.execute("SELECT * FROM RefSpectra WHERE id = ?", (spectrum_number, )).fetchone()
        spectrum = self._new_spectrum()
        spectrum.key = info['id']
        spectrum.index = info['id'] - 1
        spectrum.precursor_mz = info['precursorMZ']
        spectrum.add_attribute(CHARGE_STATE, info['precursorCharge'])
        try:
            spectrum.add_attribute(RETENTION_TIME, info['retentionTime'])
        except KeyError:
            pass

        try:
            spectrum.add_attribute(NUMBER_OF_REPLICATE_SPECTRA_AVAILABLE, info['copies'])
            spectrum.add_attribute(NUMBER_OF_REPLICATE_SPECTRA_USED, 1)
        except KeyError:
            pass

        n_peaks = info['numPeaks']
        spectrum.add_attribute("MS:1003059|number of peaks", n_peaks)

        try:
            spectrum.add_attribute(
                SOURCE_FILE,
                self.connection.execute("SELECT fileName FROM SpectrumSourceFiles WHERE id = ?",
                                        (info['fileID'], )).fetchone()['fileName']
            )
        except KeyError:
            pass
        spectrum.add_attribute(
            SCAN_NUMBER,
            info["SpecIDinFile"]
        )

        analyte = self._new_analyte(1)
        self._populate_analyte(analyte, info)

        spectrum.add_analyte(analyte)

        interp = self._new_interpretation(1)
        interp.add_analyte(analyte)
        spectrum.add_interpretation(interp)

        peak_data = self.connection.execute("SELECT * FROM RefSpectraPeaks WHERE RefSpectraID = ?", (spectrum_number, )).fetchone()
        try:
            mz_array = np.frombuffer(zlib.decompress(peak_data['peakMZ']), dtype=np.float64)
        except zlib.error:
            mz_array = np.frombuffer(peak_data['peakMZ'], dtype=np.float64)
        if mz_array.size != n_peaks:
            raise ValueError(f"m/z array does not match expected peak count ({mz_array.size} != {n_peaks})")

        try:
            intensity_array = np.frombuffer(zlib.decompress(peak_data['peakIntensity']), dtype=np.float32)
        except zlib.error:
            intensity_array = np.frombuffer(peak_data['peakIntensity'], dtype=np.float32)
        if intensity_array.size != n_peaks:
            raise ValueError(f"m/z array does not match expected peak count ({intensity_array.size} != {n_peaks})")

        peak_list = []
        for i, mz in enumerate(mz_array):
            row = (mz, intensity_array[i], [], '')
            peak_list.append(row)
        spectrum.peak_list = peak_list
        return spectrum

    def read(self) -> Iterator[Spectrum]:
        for rec in self.index:
            yield self.get_spectrum(rec.number)

