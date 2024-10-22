from dataclasses import dataclass
import io
import re

from typing import Any, Dict, List, Tuple

from mzpaf.annotation import IonAnnotationBase, tokenize_signed_symbol_list

from mzspeclib.analyte import Analyte, Interpretation, InterpretationMember
from mzspeclib.attributes import AttributeManager, Attributed

from mzspeclib.annotation import AnnotationStringParser
from mzspeclib.spectrum import Spectrum

from .utils import open_stream, CaseInsensitiveDict, try_cast
from .base import DEFAULT_VERSION, FORMAT_VERSION, LIBRARY_NAME
from .msp import (
    DispatchingAttributeHandler,
    MSPSpectralLibrary as _MSPSpectralLibrary,
    analyte_terms,
    AttributeRulesEngine,
    AttributeRuleBase,
    AttributeHandler,
    MappingAttributeHandler,
    FunctionAttributeHandler,
    PeakParsingStrategy,
    FullNameParsingRule,
    StrippedPeptideParserRule,
    collision_energy_handler,
    nreps_handler,
    mz_diff_handler,
    protein_handler as msp_protein_handler,
    PeptideNotationUnpackingRule,
)

ANNOTATION_PATTERN = re.compile(
    r"""^(?P<is_auxiliary>&)?
(?:(?:(?P<series>[abcxyz]\.?)(?P<ordinal>\d+))|
   (:?Int/(?P<series_internal>[ARNDCEQGHKMFPSTWYVILJarndceqghkmfpstwyvilj]+))|
   (?P<precursor>p)|
   (:?I(?P<immonium>[ARNDCEQGHKMFPSTWYVIL]?)A?(?:(?P<immonium_modification>CAM)|[A-Z])?)|
   (?:_(?P<external_ion>[^\s,/]+))
)
(?P<neutral_losses>[+-]\d*(?:\.\d+)?)?
(?:\^(?P<charge>[+-]?\d+))?
(?:\[M(?P<adducts>(:?[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?
(?:(?P<isotope>\d*)i)?
(?:/(?P<mass_error>[+-]?\d+(?:\.\d+))(?P<mass_error_unit>ppm)?)?
""",
    re.X,
)


HEADER_MAX_DEPTH = 1000


class SPTXTAnnotationParser(AnnotationStringParser):
    nominal_loss_to_formula = {
        "18": ["H2O"],
        "17": ["NH3"],
        "28": ["CO"],
        "34": ["NH3", "NH3"],
        "35": ["H2O", "NH3"],
        "36": ["H2O", "H2O"],
        "44": ["CO2"],
        "43": ["HNCO"],
        "45": ["HCONH2"],
        "46": ["HCOOH"],
        "64": ["CH4OS"],
        "82": ["CH4OS", "H2O"],
        "91": ["C2H5NOS"],
        "92": ["C2H4O2S"],
    }

    def __init__(self, pattern):
        super().__init__(pattern)

    def _dispatch(
        self,
        annotation_string,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        data["sequence_ordinal"] = None
        return super()._dispatch(
            annotation_string,
            data,
            adducts,
            charge,
            isotope,
            neutral_losses,
            analyte_reference,
            mass_error,
            confidence,
            **kwargs,
        )

    def parse_annotation(self, annotation_string: str, **kwargs) -> List[IonAnnotationBase]:
        parts = annotation_string.split(",")
        patched_token = []
        for part in parts:
            if part.startswith("[") and part.endswith("]"):
                patched_token.append("&" + part[1:-1])
            else:
                patched_token.append(part)
        annotation_string = ','.join(patched_token)
        return super().parse_annotation(annotation_string, **kwargs)

    def _convert_losses(self, neutral_losses: str):
        tokens = []
        for i, loss in enumerate(tokenize_signed_symbol_list(neutral_losses)):
            if loss[0] == '-':
                sign = '-'
                loss = loss[1:]
            else:
                sign = '+'
            if loss in self.nominal_loss_to_formula:
                loss = sign.join(self.nominal_loss_to_formula[loss])
            tokens.append(sign)
            tokens.append(loss)
        return ''.join(tokens)


    def _parse_string(self, annotation_string: str, **kwargs) -> Tuple[re.Match, Dict[str, str]]:
        match, data = super()._parse_string(annotation_string, **kwargs)
        if data["isotope"] is not None and not data["isotope"]:
            data["isotope"] = "1"
        if data['neutral_losses']:
            try:
                data["neutral_losses"] = self._convert_losses(data['neutral_losses'])
            except KeyError as err:
                print(annotation_string)
                raise err from None
        return match, data


parse_annotation = SPTXTAnnotationParser(ANNOTATION_PATTERN)

@dataclass
class FullNameParserSPTXT(FullNameParsingRule):
    attribute = "FullName"

    name_pattern = re.compile(r"([A-Z\-\*])\.((?:[A-Z](?:\[[^\]]+\])?)+)\.([A-Z\-\*])/*([\d]*)")
    tokenizer = re.compile(r"([A-Z])(\[[^\]]+?\])?")

    def parse_value(self, value: str):
        match = self.name_pattern.match(value)
        if match is None:
            return None
        peptide = match.group(2)
        nterm = match.group(1)
        cterm = match.group(3)
        charge = match.group(4)

        peptide = ''.join([aa for aa, mod in self.tokenizer.findall(peptide)])
        return peptide, nterm, cterm, charge


def tpp_spectrum_to_usi(spectrum: str) -> str:
    tokens = spectrum.split(".")
    source_file = tokens[0]
    scan_number = int(tokens[1])
    scan_number2 = tokens[2]
    return (f"mzspec:USI000000:{source_file}:scan:{scan_number}", source_file, scan_number)


@FunctionAttributeHandler.wraps("RawSpectrum", "BestRawSpectrum")
def raw_spectrum_handler(key, value, container: Attributed):
    if key == "RawSpectrum":
        usi, source_file, scan_number = tpp_spectrum_to_usi(value)
        container.add_attribute("MS:1003203|constituent spectrum file", source_file)
        container.add_attribute("MS:1003299|contributing replicate spectrum USI", usi)
        return True
    elif key == "BestRawSpectrum":
        usi, source_file, scan_number = tpp_spectrum_to_usi(value)
        container.add_attribute("MS:1003203|constituent spectrum file", source_file)
        container.add_attribute("MS:1003299|contributing replicate spectrum USI", usi)
        return True
    else:
        return False


@FunctionAttributeHandler.wraps("RetentionTime")
def retention_time_handler(key: str, value: Any, container: Attributed):
    rt_rep, rt_min, rt_max = map(float, value.split(","))
    group_identifier = container.get_next_group_identifier()
    container.add_attribute(
        "MS:1000894|retention time",
        rt_rep,
        group_identifier,
    )
    container.add_attribute("UO:0000000|unit", "UO:0000010|second", group_identifier)
    return True


_HCD = [
    "MS:1000044|dissociation method",
    "MS:1000422|beam-type collision-induced dissociation",
]
_CID = [
    "MS:1000044|dissociation method",
    "MS:1002472|trap-type collision-induced dissociation",
]
_CAD = ["MS:1000044|dissociation method", "MS:1000133|collision-induced dissociation"]

INSTRUMENT_DISPATCH = CaseInsensitiveDict(
    {
        "it": [_CID],
        "hcd": [_HCD],
        "qtof": [_CAD],
        "qqq": [_CAD],
        "qit": [_CAD],
    }
)


@FunctionAttributeHandler.wraps("Inst")
def instrument_handler(key, value, container: Attributed):
    value = str(value)
    if "/" in value:
        _inst_id, rest = value.split("/", 1)
        inst_type, rest = rest.split(",", 1)
        for k, v in INSTRUMENT_DISPATCH[inst_type]:
            container.add_attribute(k, v)
        return True
    else:
        container.add_attribute(*_HCD)
        return True


INTERPRETATION_TERMS = CaseInsensitiveDict({"FracUnassigned": "MS:1003079|total unassigned intensity fraction"})


@FunctionAttributeHandler.wraps("DeltaCn", "DotConsensus", "Se", "Prob")
def search_engine_keys(key: str, value: Any, container: Attributed):
    container.add_attribute_group(
        [
            ["MS:1003275|other attribute name", key],
            ["MS:1003276|other attribute value", value],
        ]
    )
    return True


INTERPRETATION_MEMBER_TERMS = DispatchingAttributeHandler(
    [MappingAttributeHandler(CaseInsensitiveDict({"XCorr": "MS:1002252|Comet:xcorr"})), (search_engine_keys)]
)

OTHER_TERMS = CaseInsensitiveDict(
    {
        "Single": [
            "MS:1003065|spectrum aggregation type",
            "MS:1003066|singleton spectrum",
        ],
        "Consensus": [
            "MS:1003065|spectrum aggregation type",
            "MS:1003067|consensus spectrum",
        ],
        "Inst": instrument_handler,
        "Spec": {
            "Consensus": [
                [
                    "MS:1003065|spectrum aggregation type",
                    "MS:1003067|consensus spectrum",
                ]
            ],
            "Single": [
                [
                    "MS:1003065|spectrum aggregation type",
                    "MS:1003066|singleton spectrum",
                ]
            ],
        },
        "Scan": "MS:1003057|scan number",
        "Origfile": "MS:1003203|constituent spectrum file",
        "filename": "MS:1003203|constituent spectrum file",
        "file_name": "MS:1003203|constituent spectrum file",
        "Sample": "MS:1000002|sample name",
        "Filter": "MS:1000512|filter string",
        "FTResolution": "MS:1000028|detector resolution",
        "PrecursorIntensity": "MS:1003085|previous MS1 scan precursor intensity",
        "OrigMaxIntensity": "MS:1003086|precursor apex intensity",
        "Purity": "MS:1009013|isolation window precursor purity",
        "NumPeaks": "MS:1003059|number of peaks",
        "Run": "MS:1003203|constituent spectrum file",
        "CollisionEnergy": collision_energy_handler,
    }
)


SPTXT_SPECTRUM_ATTRIBUTE_MAP = CaseInsensitiveDict(
    {
        "TotalIonCurrent": "MS:1000285|total ion current",
        "RetentionTime": retention_time_handler,
    }
)


@FunctionAttributeHandler.wraps("Protein")
def protein_handler(key, value: str, container):
    if "/" in value:
        pid, rest = value.split("/", 1)
        if pid.isdigit():
            return msp_protein_handler(key, rest, container)
        else:
            return msp_protein_handler(key, value, container)
    else:
        return msp_protein_handler(key, value, container)


SPTXT_ANALYTE_ATTRIBUTE_MAP = CaseInsensitiveDict(
    {
        "MassDiff": mz_diff_handler,
        "Protein": protein_handler,
    }
)


SPTXT_SPECTRUM_ATTRIBUTE_HANDLER = DispatchingAttributeHandler()

SPTXT_SPECTRUM_ATTRIBUTE_HANDLER.add(nreps_handler)
SPTXT_SPECTRUM_ATTRIBUTE_HANDLER.add(MappingAttributeHandler(SPTXT_SPECTRUM_ATTRIBUTE_MAP))
SPTXT_SPECTRUM_ATTRIBUTE_HANDLER.add(raw_spectrum_handler)


@dataclass
class SkippedAttributeRule(AttributeRuleBase):
    attribute: str

    def is_applicable(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> bool:
        return self.attribute in attributes

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> List[str]:
        attributes.pop(self.attribute)
        return []


@dataclass
class LibIDAttributeRule(AttributeRuleBase):
    attribute: str

    def is_applicable(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> bool:
        return self.attribute in attributes

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> List[str]:
        spectrum.key = attributes.pop(self.attribute)
        return []


class SPTXTSpectralLibrary(_MSPSpectralLibrary):
    file_format = "sptxt"
    format_name = "sptxt"

    @classmethod
    def guess_from_header(cls, filename: str) -> bool:
        with open_stream(filename, "r") as stream:
            first_line = stream.readline()
            i = 1
            while first_line.startswith("###") and i < HEADER_MAX_DEPTH:
                first_line = stream.readline()
                i += 1
            if re.match("Name: ", first_line):
                return True
            if first_line.startswith("###"):
                return True
        return False

    def read_header(self) -> bool:
        with open_stream(self.filename, "rb") as stream:
            match, offset = self._parse_header_from_stream(stream)
            return match
        return False

    def _parse_header_from_stream(self, stream: io.RawIOBase) -> Tuple[bool, int]:
        nbytes = 0
        attributes = AttributeManager()
        attributes.add_attribute(FORMAT_VERSION, DEFAULT_VERSION)

        line = stream.readline()
        nbytes += len(line)
        i = 0
        header_buffer = []
        while not re.match(b"Name:", line) and i < HEADER_MAX_DEPTH:
            header_buffer.append(line)
            line = stream.readline()
            nbytes += len(line)
            i += 1

        self._parse_header_into_attributes(header_buffer, attributes)
        if not attributes.has_attribute(LIBRARY_NAME):
            name = self._infer_lib_name()
            if name:
                attributes.add_attribute(LIBRARY_NAME, name)

        self.attributes.clear()
        self.attributes._from_iterable(attributes)
        if re.match(b"Name:", line):
            header_buffer.append(line)

            return True, nbytes - len(line)
        return False, 0

    def _parse_header_into_attributes(self, lines: List[bytes], attributes: AttributeManager):
        for i, line in enumerate(lines):
            if line.startswith(b"###"):
                if i == 0:
                    name = line.rsplit(b'(', 1)[0].strip().decode('utf8').strip("# ")
                    attributes.add_attribute(LIBRARY_NAME, name)
                elif i == 1:
                    version = line.split(b"(", 1)[1].strip().decode("utf8").split(",")[0].split(" ")[1]
                    gid = attributes.get_next_group_identifier()
                    attributes.add_attribute("MS:1003207|library creation software", "MS:1001477|SpectraST", gid)
                    attributes.add_attribute("MS:1003200|software version", version, gid)
            else:
                continue

    def _make_peak_parsing_strategy(self) -> PeakParsingStrategy:
        return PeakParsingStrategy(parse_annotation)

    def _make_attribute_handlers(self) -> AttributeRulesEngine:
        return AttributeRulesEngine(
            MappingAttributeHandler(OTHER_TERMS),
            MappingAttributeHandler(analyte_terms),
            MappingAttributeHandler(INTERPRETATION_TERMS),
            INTERPRETATION_MEMBER_TERMS,
            SPTXT_SPECTRUM_ATTRIBUTE_HANDLER,
            MappingAttributeHandler(SPTXT_ANALYTE_ATTRIBUTE_MAP),
            early_rules=[
                FullNameParserSPTXT(),
                FullNameParserSPTXT("Fullname"),
                LibIDAttributeRule("LibID"),
                SkippedAttributeRule("BinaryFileOffset"),
                # Redundant with the Pep parsing rule
                SkippedAttributeRule("NTT"),
            ],
            late_rules=[StrippedPeptideParserRule(), PeptideNotationUnpackingRule(self.modification_parser)],
        )
