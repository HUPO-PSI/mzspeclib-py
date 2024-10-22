"""Reader for Biognosys Spectronaut TSV spectral libraries"""

import json
import os

from typing import List, Tuple, Dict, Iterator, Any, Deque, Union

from pyteomics import proforma

from mzspeclib import annotation
from mzspeclib.analyte import Analyte
from mzspeclib.backends.base import LIBRARY_NAME, _CSVSpectralLibraryBackendBase, FORMAT_VERSION, DEFAULT_VERSION
from mzspeclib.backends.utils import open_stream, urlify, try_cast
from mzspeclib.spectrum import Spectrum, SPECTRUM_NAME
from mzspeclib import const as c
from mzspeclib.const import (
    CHARGE_STATE,
    SOURCE_FILE,
    STRIPPED_PEPTIDE_SEQ as STRIPPED_PEPTIDE_TERM,
    PROFORMA_SEQ as PROFORMA_PEPTIDE_TERM,
    CUSTOM_ATTRIBUTE_NAME,
    CUSTOM_ATTRIBUTE_VALUE,
)


SELECTED_ION_MZ = "MS:1003208|experimental precursor monoisotopic m/z"

ID_SEP = "/"
_ID_SEP = ID_SEP.encode("ascii")

NO_LOSS = 'noloss'


def _rewrite_modified_peptide_as_proforma(sequence: str) -> str:
    if sequence.startswith("_"):
        sequence = sequence.strip("_")
    buffer: Deque[str] = Deque()
    last_paren = None
    for i, c in enumerate(sequence):
        if c == ']':
            # Erase any text in parentheses as these indicate the modification
            # rule and not the modificatin name. We could look at the modification
            # rule to infer N-term and C-term rules, but we don't have enough examples
            if last_paren is not None:
                k = i - last_paren
                for _ in range(k + 1):
                    buffer.pop()
                last_paren = None
                buffer.append(c)
        elif c == '(':
            last_paren = i
            buffer.append(c)
        else:
            buffer.append(c)
    pf_seq = ''.join(buffer)
    # A peptide with an N-terminal modification will start with a square brace
    # but needs to have a "-" added to be well-formed ProForma
    if pf_seq.startswith("["):
        i = pf_seq.find(']') + 1
        if i == 0:
            raise ValueError(f"Malformed peptide sequence {sequence}")
        pf_seq = f"{pf_seq[:i]}-{pf_seq[i:]}"
    return pf_seq


def _parse_value(value: str) -> Union[float, int, str, bool]:
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        lower = value.lower()
        if lower == "true":
            return True
        elif lower == "false":
            return False
        return value


class SpectronautTSVSpectralLibrary(_CSVSpectralLibraryBackendBase):
    """Read Spectronaut TSV spectral libraries."""

    format_name = "spectronaut.tsv"

    _custom_spectrum_keys = [
        "ExcludeFromAssay",
        "BGSInferenceId",
        "AllowForNormalization",
        "Workflow",
    ]

    _custom_analyte_keys = [
        "IsProteotypic",
        "FASTAName",
        "Database",
        "ProteinGroups",
    ]

    _key_columns = ['ModifiedPeptide', 'PrecursorCharge']

    _required_columns = ['PrecursorMz', 'PrecursorCharge',
                         'ModifiedPeptide', 'StrippedPeptide',
                         'FragmentMz', 'RelativeIntensity',
                         'FragmentType', 'FragmentNumber',
                         'FragmentCharge', 'FragmentLossType']

    def __init__(self, filename: str, index_type=None, **kwargs):
        super().__init__(filename, index_type=index_type, delimiter='\t', **kwargs)

    def _spectrum_origin_type(self):
        key = c.SPECTRUM_ORIGIN_TYPE
        value = c.SELECTED_FRAGMENT_THEORETICAL_OBSERVED_INTENSITY_SPECTRUM
        return key, value

    def _spectrum_aggregation_type(self):
        key = c.SPECTRUM_AGGREGATION_TYPE
        value = c.CONSENSUS_SPECTRUM
        return key, value

    def read_header(self) -> bool:
        result = super().read_header()
        self.add_attribute(FORMAT_VERSION, DEFAULT_VERSION)
        if hasattr(self.filename, 'name'):
            name = self.filename.name.replace(".gz", '').rsplit('.', 1)[0].split(os.sep)[-1]
        else:
            name = self.filename.replace(".gz", '').rsplit(".", 1)[0].split(os.sep)[-1]
        self.add_attribute(LIBRARY_NAME, name)
        self.add_attribute(c.LIBRARY_CREATION_SW, "MS:1001327|Spectronaut")
        return result

    def _batch_rows(self, iterator: Iterator[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
        group_key = None
        group = []
        i = 0
        for row in iterator:
            key = [row[k_i] for k_i in self._key_columns]
            if group_key is None:
                group_key = key
                group.append(row)
            elif group_key == key:
                group.append(row)
            else:
                i += 1
                yield group
                group = [row]
                group_key = key
        if group:
            yield group

    def create_index(self):
        key = None
        delimiter = self._delimiter.encode("ascii")
        with open_stream(self.filename, 'rb') as stream:
            header = stream.readline()
            header_cols = header.split(delimiter)
            column_keys = (
                header_cols.index(b'ModifiedPeptide'),
                header_cols.index(b'PrecursorCharge'),
            )
            offset = stream.tell()

            line = stream.readline()
            tokens = line.split(delimiter)
            key = (tokens[column_keys[0]].strip(b"_"), tokens[column_keys[1]])
            n = 0
            self.index.add(
                number=n,
                offset=offset,
                name=_ID_SEP.join(key).decode("utf8"),
                analyte=None
            )
            n += 1

            # To hold previous values
            last_offset = None
            while line:
                tokens = line.split(delimiter)
                next_key = (tokens[column_keys[0]].strip(b"_"), tokens[column_keys[1]])
                if next_key != key:
                    key = next_key
                    self.index.add(
                        number=n,
                        offset=offset,
                        name=_ID_SEP.join(key).decode("utf8"),
                        analyte=None
                    )
                    n += 1
                    last_offset = stream.tell()
                else:
                    offset = stream.tell()
                line = stream.readline()
        n += 1
        self.index.commit()
        return n

    def _generate_peaks(self, batch: List[Dict[str, Any]]) -> List[Tuple[float, float, List[annotation.IonAnnotationBase], List]]:
        peaks = []
        for row in batch:
            mz = float(row['FragmentMz'])
            intensity = float(row['RelativeIntensity'])

            series = row['FragmentType']
            ordinal = int(row['FragmentNumber'])
            charge = int(row['FragmentCharge'])

            loss_type = row['FragmentLossType']
            if loss_type != NO_LOSS:
                loss_type = annotation.NeutralName.parse('-' + loss_type)
            else:
                loss_type = None

            annot = annotation.PeptideFragmentIonAnnotation(
                series, ordinal, neutral_losses=loss_type, charge=charge,
                mass_error=annotation.MassError(0, 'Da')
            )

            peak = [
                mz, intensity, [annot], []
            ]
            peaks.append(peak)
        return peaks

    def _build_analyte(self, description: Dict[str, Any], analyte: Analyte) -> Analyte:
        pf_seq = _rewrite_modified_peptide_as_proforma(description['ModifiedPeptide'])
        peptide = proforma.ProForma.parse(pf_seq)
        peptide.charge_state = int(description['PrecursorCharge'])
        analyte.add_attribute(STRIPPED_PEPTIDE_TERM, description['StrippedPeptide'])
        analyte.add_attribute(c.PROFORMA_ION, str(peptide))
        analyte.add_attribute(c.THEORETICAL_MASS, peptide.mass)


        protein_group_id = analyte.get_next_group_identifier()


        if 'UniProtIds' in description:
            analyte.add_attribute(
                c.PROTEIN_ACCESSION,
                description['UniProtIds'],
                group_identifier=protein_group_id
            )
        if 'Protein Name' in description:
            analyte.add_attribute(
                c.PROTEIN_NAME,
                description["Protein Name"],
                group_identifier=protein_group_id
            )
        if 'ProteinDescription' in description:
            analyte.add_attribute(
                c.PROTEIN_DESCRIPTION,
                description['ProteinDescription'],
                group_identifier=protein_group_id
            )

        if "OrganismId" in description and description["OrganismId"] is not None:
            analyte.add_attribute_group([
                [c.TAXONOMY_NCBI_TAX_ID, try_cast(description['OrganismId'])],
                [c.TAXONOMY_SCIENTIFIC_NAME, description['Organisms']],
            ])

        for key in self._custom_analyte_keys:
            if key in description:
                analyte.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(description[key])]
                ])

        return analyte

    def _parse_from_buffer(self, buffer: List[Dict[str, Any]], spectrum_index: int = None) -> Spectrum:
        spec = self._new_spectrum()
        descr = buffer[0]

        key = (descr['ModifiedPeptide'].strip("_"), descr['PrecursorCharge'])

        spec.add_attribute(SPECTRUM_NAME, ID_SEP.join(key))
        spec.add_attribute(SELECTED_ION_MZ, float(descr['PrecursorMz']))
        spec.add_attribute(CHARGE_STATE, int(descr['PrecursorCharge']))
        spec.add_attribute(SOURCE_FILE, urlify(descr['ReferenceRun']))
        spec.add_attribute(*self._spectrum_origin_type())
        spec.add_attribute(*self._spectrum_aggregation_type())

        spec.add_attribute_group([
            [CUSTOM_ATTRIBUTE_NAME, "LabeledPeptide"],
            [CUSTOM_ATTRIBUTE_VALUE, descr['LabeledPeptide']]
        ])

        if 'IonMobility' in descr:
            spec.add_attribute("MS:1002476|ion mobility drift time", float(descr['IonMobility']))
        if 'CV' in descr:
            spec.add_attribute("MS:1001581|FAIMS compensation voltage", float(descr['CV']))
        if 'iRT' in descr:
            group_id = spec.get_next_group_identifier()
            spec.add_attribute("MS:1000896|normalized retention time", float(descr['iRT']), group_identifier=group_id)
            spec.add_attribute(
                c.UNIT,
                c.MINUTE,
                group_identifier=group_id,
            )

        analyte = self._new_analyte('1')
        self._build_analyte(descr, analyte)
        spec.add_analyte(analyte)

        for key in self._custom_spectrum_keys:
            if key in descr:
                spec.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(descr[key])]
                ])

        spec.peak_list = self._generate_peaks(buffer)
        spec.add_attribute(c.NUM_PEAKS, len(spec.peak_list))

        if spectrum_index is not None:
            spec.index = spectrum_index
        else:
            spec.index = -1
        spec.key = spec.index + 1
        return spec
