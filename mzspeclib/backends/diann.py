"""
Read DIA-NN TSV spectral libraries.

For more information, see `DIA-NN's website <https://github.com/vdemichev/DiaNN`_
"""
import json
import os

from typing import List, Optional, Tuple, Dict, Iterator, Any, Union

from pyteomics import proforma

from mzspeclib import annotation
from mzspeclib.backends.base import DEFAULT_VERSION, FORMAT_VERSION, LIBRARY_NAME, _CSVSpectralLibraryBackendBase
from mzspeclib.backends.utils import open_stream, urlify
from mzspeclib.spectrum import Spectrum, SPECTRUM_NAME
from mzspeclib.const import (
    PROFORMA_SEQ as PROFORMA_PEPTIDE_TERM,
    PROFORMA_ION,
    STRIPPED_PEPTIDE_SEQ as STRIPPED_PEPTIDE_TERM,
    SELECTED_ION_MZ as SPECTRUM_SELECTED_ION_MZ,
    CHARGE_STATE,
    SOURCE_FILE,
    CUSTOM_ATTRIBUTE_NAME,
    CUSTOM_ATTRIBUTE_VALUE
)


def _rewrite_unimod_peptide_as_proforma(sequence: str) -> str:
    return sequence.replace("(", '[').replace(')', ']').replace("UniMod", "UNIMOD")


THEO_SELECTED_ION_MZ = "MS:1003053|theoretical monoisotopic m/z"

NO_LOSS = 'noloss'


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


class DIANNTSVSpectralLibrary(_CSVSpectralLibraryBackendBase):
    """Reader for DIA-NN TSV spectral libraries."""

    format_name = "dia-nn.tsv"

    _custom_spectrum_keys = [
        "ExcludeFromAssay",
        "AllowForNormalization",
    ]

    _custom_analyte_keys = [
        "Proteotypic",
        "ProteinGroup",
    ]

    _required_columns = ['transition_group_id', 'PrecursorMz', 'PrecursorCharge',
                         'FullUniModPeptideName', 'ProductMz', 'LibraryIntensity',
                         'FragmentType', 'FragmentSeriesNumber',
                         'FragmentCharge', 'FragmentLossType']

    def __init__(self, filename: str, index_type=None, **kwargs):
        super().__init__(filename, index_type=index_type, delimiter='\t', **kwargs)

    def _spectrum_origin_type(self):
        key = "MS:1003072|spectrum origin type"
        value = "MS:1003074|predicted spectrum"
        return key, value

    def _spectrum_aggregation_type(self):
        key = "MS:1003065|spectrum aggregation type"
        value = "MS:1003074|predicted spectrum"
        return key, value

    def read_header(self) -> bool:
        result = super().read_header()
        self.add_attribute(FORMAT_VERSION, DEFAULT_VERSION)
        if hasattr(self.filename, 'name'):
            name = self.filename.name.replace(".gz", '').rsplit('.', 1)[0].split(os.sep)[-1]
        else:
            name = self.filename.replace(".gz", '').rsplit(".", 1)[0].split(os.sep)[-1]
        self.add_attribute(LIBRARY_NAME, name)
        self.add_attribute("MS:1003207|library creation software", "MS:1003253|DIA-NN")
        return result

    def create_index(self):
        with open_stream(self.filename, 'rb') as stream:
            header = stream.readline()
            header_cols = header.split(b'\t')
            column_key = header_cols.index(b'transition_group_id')
            offset = stream.tell()

            line = stream.readline()
            tokens = line.split(b'\t')
            key = tokens[column_key]
            index = {key: offset}
            n = 0
            self.index.add(
                number=n,
                offset=offset,
                name=key.decode("utf8"),
                analyte=None
            )
            n += 1

            # To hold previous values
            last_offset = None
            last_key = None
            while line:
                tokens = line.split(b'\t')
                if tokens[column_key] != key:
                    key = tokens[column_key]
                    self.index.add(
                        number=n,
                        offset=offset,
                        name=key.decode("utf8") ,
                        analyte=None
                    )
                    n += 1
                    last_offset = stream.tell()
                else:
                    offset = stream.tell()
                line = stream.readline()
        # self.index.add(
        #     number=n,
        #     offset=last_offset,
        #     name=tokens[column_key].decode('utf8'),
        #     analyte=None
        # )
        n += 1
        self.index.commit()
        return n

    def _parse_from_buffer(self, buffer: List[Dict[str, Any]], spectrum_index: Optional[int] = None) -> Spectrum:
        spec = self._new_spectrum()
        descr = buffer[0]

        spec.add_attribute(SPECTRUM_NAME, descr['transition_group_id'])
        spec.add_attribute(SPECTRUM_SELECTED_ION_MZ, float(descr['PrecursorMz']))

        if 'FileName' in descr:
            spec.add_attribute(SOURCE_FILE, urlify(descr['FileName']))

        if 'decoy' in descr and int(descr['decoy']):
            spec.add_attribute("MS:1003072|spectrum origin type", "MS:1003192|decoy spectrum")
        else:
            spec.add_attribute(*self._spectrum_origin_type())

        spec.add_attribute(*self._spectrum_aggregation_type())

        if 'IonMobility' in descr:
            spec.add_attribute("MS:1002476|ion mobility drift time", float(descr['IonMobility']))

        analyte = self._new_analyte('1')

        pf_seq = _rewrite_unimod_peptide_as_proforma(descr['FullUniModPeptideName'])
        peptide = proforma.ProForma.parse(pf_seq)

        if 'PeptideSequence' in descr:
            analyte.add_attribute(STRIPPED_PEPTIDE_TERM, descr['PeptideSequence'])
        peptide.charge_state = descr['PrecursorCharge']
        analyte.add_attribute(PROFORMA_ION, str(peptide))
        analyte.add_attribute("MS:1001117|theoretical mass", peptide.mass)
        spec.add_attribute(CHARGE_STATE, int(descr['PrecursorCharge']))

        protein_group_id = analyte.get_next_group_identifier()
        if "UniprotID" in descr:
            analyte.add_attribute(
                "MS:1000885|protein accession",
                descr['UniprotID'],
                group_identifier=protein_group_id
            )
        if "ProteinName" in  descr:
            if descr['ProteinName']:
                analyte.add_attribute(
                    "MS:1000886|protein name",
                    descr["ProteinName"],
                    group_identifier=protein_group_id
                )

        for key in self._custom_analyte_keys:
            if key in descr:
                analyte.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(descr[key])]
                ])

        spec.add_analyte(analyte)
        spec.peak_list = self._generate_peaks(buffer)
        spec.add_attribute("MS:1003059|number of peaks", len(spec.peak_list))

        for key in self._custom_spectrum_keys:
            if key in descr:
                spec.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(descr[key])]
                ])

        if spectrum_index is not None:
            spec.index = spectrum_index
        else:
            spec.index = -1
        spec.key = spec.index + 1
        return spec

    def _generate_peaks(self, batch: List[Dict[str, Any]]) -> List[Tuple[float, float, List[annotation.IonAnnotationBase], List]]:
        peaks = []
        for row in batch:
            mz = float(row['ProductMz'])
            intensity = float(row['LibraryIntensity'])

            series = row['FragmentType']
            ordinal = int(row['FragmentSeriesNumber'])
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

    def _batch_rows(self, iterator: Iterator[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
        group_key = None
        group = []
        for row in iterator:
            key = row['transition_group_id']
            if group_key is None:
                group_key = key
                group.append(row)
            elif group_key == key:
                group.append(row)
            else:
                yield group
                group = [row]
                group_key = key
        if group:
            yield group


DiaNNTSVSpectralLibrary = DIANNTSVSpectralLibrary
