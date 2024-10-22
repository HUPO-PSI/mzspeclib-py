"""
Read NIST MSP formats.

There are an uncountable number of dialects of MSP. This library attempts to
cover a representative subset from the NIST and some found in the wild. If you
encounter an MSP file that does not parse properly, please report it.
"""

import re
import io
import os
import math
import logging
import itertools
import warnings

from dataclasses import dataclass, field
from fractions import Fraction

from typing import (
    Any,
    Callable,
    Collection,
    Dict,
    List,
    Mapping,
    Optional,
    Set,
    Tuple,
    Iterable,
    DefaultDict,
)

from pyteomics import proforma

from mzspeclib import annotation

from mzspeclib.analyte import (
    FIRST_ANALYTE_KEY,
    FIRST_INTERPRETATION_KEY,
    Analyte,
    Interpretation,
    InterpretationMember,
    ProteinDescription,
)

from mzspeclib import const
from mzspeclib.const import (
    CHARGE_STATE,
    CONSENSUS_SPECTRUM,
    MOLECULAR_FORMULA,
    MOLECULAR_MASS,
    PROTON,
    PROFORMA_ION,
    PROFORMA_SEQ,
    SCAN_NUMBER,
    SINGLETON_SPECTRUM,
    SOURCE_FILE,
    SPECTRUM_AGGREGATION_TYPE,
    STRIPPED_PEPTIDE_SEQ as STRIPPED_PEPTIDE_TERM,
    PEAK_ATTRIB,
    PEAK_OBSERVATION_FREQ,
    SELECTED_ION_MZ,
    PRECURSOR_MZ,
    THEORETICAL_MZ,
    CUSTOM_ATTRIBUTE_NAME,
    CUSTOM_ATTRIBUTE_VALUE,
    THEORETICAL_MASS,
    INTENSITY_OF_HIGH_UNASSIGNED_PEAK,
    TOP_20_UNASSIGNED_INTENSITY_FRACTION,
    TOTAL_UNASSIGNED_INTENSITY_FRACTION,
    NUM_UNASSIGNED_PEAKS_IN_TOP_20,
    NUM_PEAKS,
)
from mzspeclib.spectrum import Spectrum, SPECTRUM_NAME
from mzspeclib.attributes import Attribute, AttributeManager, AttributeSet, Attributed

from .base import (
    DEFAULT_VERSION,
    FORMAT_VERSION,
    _PlainTextSpectralLibraryBackendBase,
    LIBRARY_NAME,
    AttributeSetTypes,
    SpectralLibraryBackendBase,
    SpectralLibraryWriterBase,
)
from .utils import try_cast, open_stream, CaseInsensitiveDict


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# TODO: Name: could be Compound or SpectrumName
leader_terms = {
    "Name": SPECTRUM_NAME,
    "NAME": SPECTRUM_NAME,
    "Compound": SPECTRUM_NAME,
    "COMPOUND": SPECTRUM_NAME,
}


def _generate_numpeaks_keys():
    w1 = "num"
    w2 = "peaks"
    seps = (" ", "")
    cases = (str.lower, str.title, str.upper)
    w1_cases = [c(w1) for c in cases]
    w2_cases = [c(w2) for c in cases]
    return {
        (w1c + sep + w2c)
        for (sep, w1c, w2c) in itertools.product(seps, w1_cases, w2_cases)
    }


NUM_PEAKS_KEYS = _generate_numpeaks_keys()

LEADER_TERMS_PATTERN = re.compile(r"(Name|NAME|Compound|COMPOUND)\s*:")
LEADER_TERMS_LINE_PATTERN = re.compile(r"(?:Name|NAME|Compound|COMPOUND)\s*:\s+(.+)")

SPACE_SPLITTER = re.compile(r"\s+")

PEPTIDE_MODIFICATION_TERM = "MS:1001471|peptide modification details"


RawPeakLine = Tuple[float, float, str, str]
PeakLine = Tuple[float, float, List[annotation.IonAnnotationBase], List[Any]]

PeakAggregateParseFn = Callable[[str], Any]


class _AttributeHelperBase:
    def add_value(self, key: str, value: Any, container: Attributed):
        """
        Add a single attribute to an :class:`~.Attributed` container.

        Explicitly check for whether the attribute is already present with
        the same value and skip the operation to avoid duplication.
        """
        if container.has_attribute(key):
            existing = container.get_attribute(key)
            if existing == value:
                return
        container.add_attribute(key, value)

    def add_group(self, keys: List[str], values: List[Any], container: Attributed):
        """
        Add a collection of attributes to an :class:`~.Attributed` container
        that will be part of a single group.
        """
        group_id = container.get_next_group_identifier()
        for k, v in zip(keys, values):
            container.add_attribute(k, v, group_id)

    def deduplicate_attribute(self, container: Attributed, attribute_name: str) -> int:
        try:
            attrs = container.get_attribute(attribute_name, raw=True)
        except KeyError:
            return 0
        if isinstance(attrs, list):
            # Sometimes there are multiple entries that map to molecular mass, and one is the nominal
            # mass. Prefer listing the most precise mass first. Apply the same logic to all attributes.
            attrs.sort(key=lambda x: len(str(x.value)), reverse=True)
            container.remove_attribute(attrs[0].key)
            seen_values = set()
            for attr in attrs:
                key = (attr.value, attr.group_id)
                if key in seen_values:
                    continue
                seen_values.add(key)
                container.add_attribute(attr.key, attr.value, attr.group_id)
            return len(attrs) - len(seen_values)
        return 0


class AttributeHandler(_AttributeHelperBase):
    """
    A base type for managing conversion of specific keys or sets of keys and their values
    to controlled vocabulary attributes.
    """

    keys: Collection[str]

    def __init__(self, keys: Collection[str]):
        if isinstance(keys, str):
            keys = (keys,)
        self.keys = keys

    def __contains__(self, key):
        return key in self.keys

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        """
        Determine what, if any, and how to add this attribute to an
        :class:`~.Attributed` container.

        .. note::
            If this operation is successful, return :const:`True`, otherwise
            return :const:`False` to signal that another :class:`AttributeHandler`
            should try to handle this input, like "bubbling" in event handling
            systems.

        The :meth:`__call__` syntax is short-hand for this method.

        Arguments
        ---------
        key : :class:`str`
            The MSP comment key to map to an attribute
        value : :class:`Any`
            The interpreted value of the MSP comment. It may be a :class:`str` or
            already coerced into another value like a number.
        container : :class:`~.Attributed`
            The object to add the translated attribute to if appropriate.

        Returns
        -------
        :class:`bool`
        """
        raise NotImplementedError()

    def __call__(self, key: str, value: Any, container: Attributed) -> bool:
        """See :meth:`handle`"""
        return self.handle(key, value, container)

    def chain(self, handler: "AttributeHandler") -> "AttributeHandlerChain":
        """
        Combine this :class:`AttributeHandler` with another :class:`AttributeHandler`
        to be applied in sequence.

        The :meth:`__and__` syntax is short-hand for this method.

        Returns
        -------
        :class:`AttributeHandlerChain`
        """
        return AttributeHandlerChain([self, handler])

    def __and__(self, handler: "AttributeHandler") -> "AttributeHandlerChain":
        return self.chain(handler)


class MappingAttributeHandler(AttributeHandler):
    """
    An :class:`AttributeHandler` which uses a :class:`Mapping` to conveniently
    define the relationship between a comment key and the appropriate controlled
    vocabulary term.

    The behavior is handled differently depending upon the value in ``self.keys[key]``:
        - Scalar: ``self.keys[key]`` is the controlled vocabulary term used,
          and the comment value is minimally processed.
        - List of 2 Items and comment value is :const:`None`: ``self.keys[key][0]`` is the controlled vocabulary term
          and ``self.keys[key][1]`` is the value.
        - :class:`AttributeHandler`: All handling is delegated to ``self.keys[key](key, value, container)``
        - Mapping: The comment value is used to find the appropriate scalar or list to add.
        - List of Pairs: All pairs from the list are added with :meth:`add_group`.

    More complex transformations should be defined using an explicit rule
    like :class:`FunctionAttributeHandler` or some other sub-class if appropriate.
    """

    keys: Dict[str, Any]

    def __init__(self, keys: Dict[str, Any]):
        self.keys = keys

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        trans_key = self.keys[key]
        if value is None:
            if isinstance(trans_key, list):
                k, v = trans_key
                v = try_cast(v)
                self.add_value(k, v, container)
            else:
                return False
        elif isinstance(trans_key, str):
            self.add_value(trans_key, try_cast(value), container)
        elif isinstance(trans_key, dict):
            if value in trans_key:
                # If the mapping is a plain string, add it
                if isinstance(trans_key[value], str):
                    key, value = trans_key[value].split("=")
                    self.add_value(key, try_cast(value), container)
                # Or if it is a list, then there are multiple terms to add within a group
                elif isinstance(trans_key[value], list):
                    if len(trans_key[value]) == 1:
                        for item in trans_key[value]:
                            self.add_value(item[0], try_cast(item[1]), container)
                    else:
                        k, v = zip(*trans_key[value])
                        self.add_group(k, map(try_cast, v), container)
            else:
                return False
        elif isinstance(trans_key, AttributeHandler):
            return trans_key(key, value, container)
        return True


class AttributeHandlerChain:
    """
    An :class:`AttributeHandler` that sequentially applies a list
    of :class:`AttributeHandler`s until it finds one that works.

    This type exists to make handler composition easier.
    """

    chain: List[AttributeHandler]

    def __init__(self, chain: List[AttributeHandler]):
        self.chain = chain

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        result = True
        for handler in self.chain:
            result |= handler.handle(key, value, container)
            if not result:
                return result
        return result

    def __contains__(self, key):
        for handler in self.chain:
            if key in handler:
                return True
        return False

    def __iter__(self):
        return iter(self.chain)

    def __call__(self, key: str, value: Any, container: Attributed) -> bool:
        return self.handle(key, value, container)

    def chain(self, handler: "AttributeHandler") -> "AttributeHandlerChain":
        return self.__class__(self.chain + [handler])

    def add(self, handler: AttributeHandler):
        self.chain.append(handler)
        return handler


class FunctionAttributeHandler(AttributeHandler):
    func: Callable[[str, Any, Attributed], bool]

    def __init__(
        self, keys: Collection[str], func: Callable[[str, Any, Attributed], bool]
    ):
        super().__init__(keys)
        self.func = func

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        return self.func(key, value, container)

    @classmethod
    def wraps(cls, *keys):
        def wrapper(func):
            return cls(keys, func)

        return wrapper


class DispatchingAttributeHandler(AttributeHandlerChain):
    mapping: Dict[str, AttributeHandler]

    def __init__(self, chain: List[AttributeHandler] = None):
        if not chain:
            chain = []
        super().__init__(chain)
        self.mapping = CaseInsensitiveDict()
        for handler in self:
            for key in handler.keys:
                self.mapping[key] = handler

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        handler = self.mapping[key]
        return handler(key, value, container)

    def __contains__(self, key):
        return key in self.mapping

    def add(self, handler: AttributeHandler):
        super().add(handler)
        for key in handler.keys:
            self.mapping[key] = handler
        return handler


spectrum_terms = CaseInsensitiveDict(
    {
        "Charge": CHARGE_STATE,
        "precursor_charge": CHARGE_STATE,
        "precursorcharge": CHARGE_STATE,
        "Parent": SELECTED_ION_MZ,
        "PrecursorMonoisoMZ": PRECURSOR_MZ,
        "ObservedPrecursorMZ": PRECURSOR_MZ,
        "PrecursorMZ": PRECURSOR_MZ,
        "PRECURSORMZ": PRECURSOR_MZ,
        "precursor": PRECURSOR_MZ,
        "precursor_mass": PRECURSOR_MZ,
        "precursormass": PRECURSOR_MZ,
        "Mz_exact": PRECURSOR_MZ,
    }
)

analyte_terms = CaseInsensitiveDict(
    {
        # "MW": "MS:1000224|molecular mass",
        # "total exact mass": "MS:1000224|molecular mass",
        # "ExactMass": "MS:1000224|molecular mass",
        # "exact_mass": "MS:1000224|molecular mass",
        # "exact mass": "MS:1000224|molecular mass",
        "molecular formula": MOLECULAR_FORMULA,
        "Formula": MOLECULAR_FORMULA,
        "formula": MOLECULAR_FORMULA,
        "SMILES": "MS:1000868|SMILES formula",
        "InChIKey": "MS:1002894|InChIKey",
        "InChI": "MS:1003403|InChI",
        "Precursor_type": "MS:1002813|adduct ion formula",
        # "Theo_mz_diff": "MS:1003209|monoisotopic m/z deviation",
        "Scan": {
            "Mods": PEPTIDE_MODIFICATION_TERM,
            "Naa": "MS:1003043|number of residues",
        },
        "Pep": {
            "Tryptic": [
                ["MS:1003048|number of enzymatic termini", 2],
                ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"],
            ],
            "SemiTryptic": [
                ["MS:1003048|number of enzymatic termini", 1],
                ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"],
            ],
            "N-Semitryptic": [
                ["MS:1003048|number of enzymatic termini", 1],
                ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"],
            ],
            "C-Semitryptic": [
                ["MS:1003048|number of enzymatic termini", 1],
                ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"],
            ],
            "Tryptic/miss_good_confirmed": [
                ["MS:1003048|number of enzymatic termini", 2],
                ["MS:1003044|number of missed cleavages", "0"],
                ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"],
            ],
            "Tryptic/miss_bad_confirmed": [
                ["MS:1003048|number of enzymatic termini", 2],
                ["MS:1003044|number of missed cleavages", ">0"],
                ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"],
            ],
        },
        "MC": "MS:1003044|number of missed cleavages",
        "NMC": "MS:1003044|number of missed cleavages",
        "Mods": PEPTIDE_MODIFICATION_TERM,
        "Naa": "MS:1003043|number of residues",
        "Mz_av": "MS:1003054|theoretical average m/z",
    }
)


_HCD = [
    const.DISSOCIATION_METHOD,
    const.HCD,
]

# TODO: qtof -> CAD
instrument_dispatch = CaseInsensitiveDict(
    {
        "it": [
            [
                const.DISSOCIATION_METHOD,
                const.TRAP_CID,
            ]
        ],
        "hcd": [_HCD],
        "QExactive": [_HCD],
        "Elite": [_HCD],
    }
)


other_terms = CaseInsensitiveDict(
    {
        "Single": [
            const.SPECTRUM_AGGREGATION_TYPE,
            const.SINGLETON_SPECTRUM,
        ],
        "Consensus": [
            const.SPECTRUM_AGGREGATION_TYPE,
            const.CONSENSUS_SPECTRUM,
        ],
        "Inst": instrument_dispatch,
        "Instrument_type": instrument_dispatch,
        "Spec": {
            "Consensus": [
                [
                    const.SPECTRUM_AGGREGATION_TYPE,
                    const.CONSENSUS_SPECTRUM,
                ]
            ]
        },
        "Scan": SCAN_NUMBER,
        "Origfile": SOURCE_FILE,
        "filename": SOURCE_FILE,
        "file_name": SOURCE_FILE,
        "Sample": const.SAMPLE_NAME,
        "Filter": const.FILTER_STRING,
        "FTResolution": "MS:1000028|detector resolution",
        "ms1PrecursorAb": "MS:1003085|previous MSn-1 scan precursor intensity",
        "Precursor1MaxAb": const.PRECURSOR_APEX_INTENSITY,
        "Purity": const.ISOLATION_WINDOW_PRECURSOR_PURITY,
        "Num peaks": NUM_PEAKS,
        "Num Peaks": NUM_PEAKS,
        "num_peaks": NUM_PEAKS,
        "numpeaks": NUM_PEAKS,
        "Run": SOURCE_FILE,
        "Splash": "MS:1002599|splash key",
    }
)


@FunctionAttributeHandler.wraps("num_unassigned_peaks")
def unassigned_peaks_handler(key: str, value: str, container: Attributed) -> bool:
    is_top_20 = False
    if isinstance(value, str):
        if "/" in value:
            if "/20" in value:
                is_top_20 = True
            value = value.split("/")[0]
        value = int(value)
    assert isinstance(value, int)
    if is_top_20:
        container.add_attribute(
            NUM_UNASSIGNED_PEAKS_IN_TOP_20, value
        )
    else:
        container.add_attribute("MS:1003288|number of unassigned peaks", value)
    return True


interpretation_terms = CaseInsensitiveDict(
    {
        "Unassigned_all_20ppm": TOTAL_UNASSIGNED_INTENSITY_FRACTION,
        "Unassign_all": TOTAL_UNASSIGNED_INTENSITY_FRACTION,
        "top_20_num_unassigned_peaks_20ppm": unassigned_peaks_handler,
        "num_unassigned_peaks_20ppm": unassigned_peaks_handler,
        "num_unassigned_peaks": unassigned_peaks_handler,
        "max_unassigned_ab_20ppm": INTENSITY_OF_HIGH_UNASSIGNED_PEAK,
        "max_unassigned_ab": INTENSITY_OF_HIGH_UNASSIGNED_PEAK,
        "Unassigned_20ppm": TOP_20_UNASSIGNED_INTENSITY_FRACTION,
        "Unassigned": TOP_20_UNASSIGNED_INTENSITY_FRACTION,
    }
)


interpretation_member_terms = CaseInsensitiveDict(
    {
        "Q-value": const.Q_VALUE,
    }
)


species_map = {
    "human": [
        [const.TAXONOMY_NCBI_TAX_ID, "NCBITaxon:9606|Homo sapiens"],
        [const.TAXONOMY_SCIENTIFIC_NAME, "Homo sapiens"],
        [const.TAXONOMY_COMMON_NAME, "human"],
    ],
    "zebrafish": [
        [const.TAXONOMY_NCBI_TAX_ID, "NCBITaxon:7955|Danio rerio"],
        [const.TAXONOMY_SCIENTIFIC_NAME, "Danio rerio"],
        [const.TAXONOMY_COMMON_NAME, "zebra fish"],
    ],
    "chicken": [
        [const.TAXONOMY_NCBI_TAX_ID, "NCBITaxon:9031|Gallus gallus"],
        [const.TAXONOMY_SCIENTIFIC_NAME, "Gallus gallus"],
        [const.TAXONOMY_COMMON_NAME, "chicken"],
    ],
    "rat": [
        [const.TAXONOMY_NCBI_TAX_ID, "NCBITaxon:10116|Rattus norvegicus"],
        [const.TAXONOMY_SCIENTIFIC_NAME, "Rattus norvegicus"],
        [const.TAXONOMY_COMMON_NAME, "rat"],
    ],
}


MODIFICATION_NAME_MAP = {
    "CAM": "Carbamidomethyl",
    "Pyro_glu": "Glu->pyro-Glu",  # Resolves UNIMOD ambiguity
    "Pyro-glu": "Gln->pyro-Glu",
    "Oxidation": "Oxidation",
    "Phospho": "Phospho",
    "TMT6plex": "TMT6plex",
    "iTRAQ": "iTRAQ",
    "Acetyl": "Acetyl",
    "TMT": "TMT6plex",
}

TERMINAL_MODIFICATIONS = {"Acetyl", "TMT6plex"}


# TODO: ppm is unsigned, add mass calculation to determine true mass accuracy

annotation_pattern = re.compile(
    r"""^
(?:(?:(?P<series>[abyxcz]\.?)(?P<ordinal>\d+))|
   (:?(?P<series_internal>(?:Int/[ARNDCEQGHKMFPSTWYVILJarndceqghkmfpstwyvilj]+)|(?:[mM](?P<internal_start>\d+):(?P<internal_end>\d+))))|
   (?P<precursor>p)|
   (:?I(?P<immonium>[ARNDCEQGHKMFPSTWYVIL])(?:(?P<immonium_modification>CAM)|[A-Z])?)|
   (?P<reporter>(?P<reporter_label>TMT|ITRAQ|iTRAQ)(?P<reporter_mass>\d+[NC]?))|
   (?:_(?P<external_ion>[^\s,/]+))|
   (?P<unannotated>\?)
)
(?P<neutral_losses>((?:[+-]\d*[A-Z][A-Za-z0-9]*)|(?:[+-]iTRAQ|TMT)|(?:[+-]\d+?(?!i)))+)?
(?:(?:\[M(?P<adducts>(:?[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?
(?:\^(?P<charge>[+-]?\d+))?
(?:(?P<isotope>[+-]?\d*)i)?)+
(?:/(?:Dev=)?(?P<mass_error>[+-]?\d+(?:\.\d+))(?P<mass_error_unit>ppm)?)?
""",
    re.X,
)


class MSPAnnotationStringParser(annotation.AnnotationStringParser):
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
        data["sequence_internal"] = None
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

    def _dispatch_internal_peptide_fragment(
        self,
        data: Dict[str, Any],
        adducts: List,
        charge: int,
        isotope: int,
        neutral_losses: List,
        analyte_reference: Any,
        mass_error: Any,
        **kwargs,
    ):
        if data["internal_start"]:
            # The hybrid pattern detected an m<start>:<end> type string, not an Int/seq string
            return super(
                MSPAnnotationStringParser, self
            )._dispatch_internal_peptide_fragment(
                data,
                adducts,
                charge,
                isotope,
                neutral_losses,
                analyte_reference,
                mass_error,
                **kwargs,
            )

        spectrum = kwargs.get("spectrum")
        if spectrum is None:
            raise ValueError(
                "Cannot infer sequence coordinates from MSP internal fragmentation notation without"
                " a reference to the spectrum, please pass spectrum as a keyword argument"
            )
        sequence = self._get_peptide_sequence_for_analyte(spectrum, analyte_reference)
        subseq = data["series_internal"]
        if subseq.startswith("Int/"):
            subseq = subseq[4:]
        subseq = subseq.upper()
        try:
            start_index = sequence.index(subseq)
        except ValueError as err:
            raise ValueError(
                f"Cannot locate internal subsequence {subseq} in {sequence}"
            ) from err
        end_index = start_index + len(subseq)
        data["internal_start"] = start_index + 1
        data["internal_end"] = end_index
        return super(
            MSPAnnotationStringParser, self
        )._dispatch_internal_peptide_fragment(
            data,
            adducts,
            charge,
            isotope,
            neutral_losses,
            analyte_reference,
            mass_error,
            **kwargs,
        )

    def _coerce_isotope(self, data):
        value = data.get("isotope")
        if value is not None:
            if value == "":
                data["isotope"] = "1"
        return super()._coerce_isotope(data)

    def _coerce_neutral_losses(self, data: Dict[str, str]) -> List:
        return super()._coerce_neutral_losses(data)

    def _coerce_analyte_reference(self, data: Dict[str, str]) -> str:
        return None

    def _dispatch_immonium(
        self,
        data: Dict[str, Any],
        adducts: List,
        charge: int,
        isotope: int,
        neutral_losses: List,
        analyte_reference: Any,
        mass_error: Any,
        **kwargs,
    ):
        modification = data["immonium_modification"]
        if modification is not None:
            try:
                modification = MODIFICATION_NAME_MAP[modification]
                data["immonium_modification"] = modification
            except KeyError as err:
                print(f"Failed to convert immonium ion modification {modification}")
        return super(MSPAnnotationStringParser, self)._dispatch_immonium(
            data,
            adducts,
            charge,
            isotope,
            neutral_losses,
            analyte_reference,
            mass_error,
            **kwargs,
        )

    def _get_peptide_sequence_for_analyte(
        self, spectrum: Spectrum, analyte_reference: Optional[Any] = None
    ) -> str:
        if analyte_reference is None:
            if len(spectrum.analytes) == 0:
                return None
            else:
                analyte_reference = spectrum.analytes[FIRST_ANALYTE_KEY].id
        analyte = spectrum.analytes.get(analyte_reference)
        if analyte is None:
            return None
        return analyte.get_attribute("MS:1000888|stripped peptide sequence")


parse_annotation = MSPAnnotationStringParser(annotation_pattern)
MODIFICATION_LIST_PARSER = re.compile(
    r"(\d+),([ARNDCEQGHKMFPSTWYVIL_\-]),([A-Za-z0-9_\-]+)"
)


class ModificationParser:
    pattern: re.Pattern
    modification_map: Dict[str, str]
    unknown_modifications: Set[str]

    def __init__(self, pattern: str = None, modification_map: Dict[str, str] = None):
        if pattern is None:
            pattern = MODIFICATION_LIST_PARSER.pattern
        if modification_map is None:
            modification_map = dict(MODIFICATION_NAME_MAP)
        self.pattern = re.compile(pattern)
        self.modification_map = modification_map
        self.unknown_modifications = set()

    def __call__(self, text: str):
        return self.parse(text)

    def parse(self, text: str) -> List[Tuple[int, str, str]]:
        if not isinstance(text, str) or not text:
            return []
        i = 0
        n = len(text)

        mods = []
        while i < n and text[i].isdigit():
            i += 1

        for position, residue, mod in self.pattern.findall(text):
            position = int(position)
            if mod not in self.modification_map:
                warnings.warn(
                    f"{mod} is not found in the known MSP modification mapping. Using this name verbatim"
                )
                modification_name = mod
                self.modification_map[mod] = mod
                self.unknown_modifications.add(mod)
            else:
                modification_name = self.modification_map[mod]
            mods.append((position, residue, modification_name))
        return mods


def null_handler(key: str, value: str, container: Attributed) -> bool:
    return True


msp_spectrum_attribute_handler = DispatchingAttributeHandler()
msp_analyte_attribute_handler = DispatchingAttributeHandler()

msp_interpretation_attribute_handler = DispatchingAttributeHandler()
msp_interpretation_attribute_handler.add(MappingAttributeHandler(interpretation_terms))

msp_spectrum_attribute_handler.add(FunctionAttributeHandler("Nprot", null_handler))
msp_spectrum_attribute_handler.add(FunctionAttributeHandler("Peptype", null_handler))

msp_spectrum_attribute_handler.add(MappingAttributeHandler(spectrum_terms))


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("Spectrum_type")
def ms_level_handler(key: str, value: str, container: Attributed) -> bool:

    attr_name = const.MS_LEVEL
    if value is None:
        return False
    if isinstance(value, str):
        if value.lower().startswith("ms"):
            value = int(value[2:])
    container.add_attribute(attr_name, value)
    return True


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps(
    "ionmode", "ion_mode", "ionization mode", "IONMODE", "Ion_mode", "MS_mode"
)
def polarity_handler(key: str, value: str, container: Attributed) -> bool:
    polarity_term = const.SCAN_POLARITY
    positive = const.POSITIVE_SCAN
    negative = const.NEGATIVE_SCAN


    if isinstance(value, str):
        value = value.lower().strip("\"'")
        if value == "positive" or value == "+":
            value = positive
        elif value == "negative" or value == "-":
            value = negative
        else:
            raise ValueError(f"Can't infer scan polarity from {value}")
        container.add_attribute(polarity_term, value)
        return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("HCD")
def dissociation_method_handler(key: str, value: str, container: Attributed) -> bool:
    container.add_attribute(
        const.DISSOCIATION_METHOD,
        const.HCD,
    )
    found_match = False
    if value is not None:
        match = re.match(r"([\d\.]+)\s*ev", value, flags=re.IGNORECASE)
        if match is not None:
            found_match = True
            group_identifier = container.get_next_group_identifier()

            container.add_attribute(
                const.COLLISION_ENERGY,
                try_cast(match.group(1)),
                group_identifier,
            )
            container.add_attribute(
                const.UNIT, const.ELECTRONVOLT, group_identifier
            )
        match = re.match(r"([\d\.]+)\s*%", value)
        if match is not None:
            found_match = True
            group_identifier = container.get_next_group_identifier()
            container.add_attribute(
                const.COLLISION_ENERGY,
                try_cast(match.group(1)),
                group_identifier,
            )
            container.add_attribute(const.UNIT, const.PERCENT, group_identifier)
    return found_match


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps(
    "Collision_energy", "CE", "colenergy", "collisionenergy", "ionization energy"
)
def collision_energy_handler(key: str, value: str, container: Attributed) -> bool:
    if isinstance(value, str):
        if "NCE" in value:
            return normalized_collision_energy_handler(key, value, container)
        match = re.match(r"([\d\.]+)(:?eV|EV|ev|Ev)?", value.strip('"'))
        if match is not None:
            value = try_cast(match.group(1))
        else:
            warnings.warn(
                f"Failed to parse {value} for {key} in collision_energy_handler"
            )
    if value is not None:
        group_identifier = container.get_next_group_identifier()
        container.add_attribute(const.COLLISION_ENERGY, value, group_identifier)
        container.add_attribute(
            const.UNIT, const.ELECTRONVOLT, group_identifier
        )
        return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("NCE", "nce")
def normalized_collision_energy_handler(
    key: str, value: str, container: Attributed
) -> bool:
    if isinstance(value, str):
        match = re.match(r"([\d\.]+)", value)
        if match is not None:
            value = try_cast(match.group(1))
    if value is not None:
        group_identifier = container.get_next_group_identifier()
        container.add_attribute(
            "MS:1000138|normalized collision energy", value, group_identifier
        )
        container.add_attribute(
            const.UNIT, const.PERCENT, group_identifier
        )
        return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("RT", "rettime", "retentiontime", "rtinseconds")
def rt_handler(key, value, container) -> bool:
    if not isinstance(value, str):
        if key.lower() == "rtinseconds":
            group_identifier = container.get_next_group_identifier()
            container.add_attribute(
                const.RETENTION_TIME,
                try_cast(value),
                group_identifier,
            )
            container.add_attribute(
                const.UNIT, const.SECOND, group_identifier
            )
        else:
            warnings.warn(f"Unable to infer retention time unit from {key!r}")
            container.add_attribute(const.RETENTION_TIME, try_cast(value))
        return True
    else:
        match = re.match(r"([\d\.]+)\s*(\D*)", value)
        if match is not None:
            if match.group(2):
                container.add_attribute(
                    "ERROR", f"Need more RT parsing code to handle this value"
                )
                return False
            else:
                group_identifier = container.get_next_group_identifier()
                container.add_attribute(
                    const.RETENTION_TIME,
                    try_cast(match.group(1)),
                    group_identifier,
                )
                #### If the value is greater than 250, assume it must be seconds
                if float(match.group(1)) > 250 or key.lower() == "rtinseconds":
                    container.add_attribute(
                        const.UNIT, const.SECOND, group_identifier
                    )
                #### Although normally assume minutes
                else:
                    container.add_attribute(
                        const.UNIT, const.MINUTE, group_identifier
                    )
                return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("ms2IsolationWidth")
def isolation_width_handler(key, value, container) -> bool:
    if value is None:
        return False
    group_identifier = container.get_next_group_identifier()
    container.add_attribute(
        "MS:1000828|isolation window lower offset", (float(value) / 2), group_identifier
    )
    container.add_attribute(const.UNIT, const.MZ, group_identifier)
    group_identifier = container.get_next_group_identifier()
    container.add_attribute(
        "MS:1000829|isolation window upper offset", (float(value) / 2), group_identifier
    )
    container.add_attribute(const.UNIT, const.MZ, group_identifier)
    return True


@msp_interpretation_attribute_handler.add
@FunctionAttributeHandler.wraps("Mz_diff", "Theo_mz_diff")
def mz_diff_handler(key, value, container: Attributed) -> bool:
    if isinstance(value, float):
        # We must be dealing with a unit-less entry.
        group_identifier = container.get_next_group_identifier()
        container.add_attribute(const.DELTA_MZ, abs(value), group_identifier)
        container.add_attribute(const.UNIT, const.MZ, group_identifier)
    else:
        match = re.match(r"([\-\+e\d\.]+)\s*ppm", value, flags=re.IGNORECASE)
        if match is not None:
            group_identifier = container.get_next_group_identifier()
            container.add_attribute(
                const.DELTA_MZ, try_cast(match.group(1)), group_identifier
            )
            container.add_attribute(
                const.UNIT, const.PPM, group_identifier
            )
        else:
            match = re.match(r"([\-\+e\d\.]+)\s*", value)
            if match is not None:
                group_identifier = container.get_next_group_identifier()
                container.add_attribute(
                    const.DELTA_MZ, try_cast(match.group(1)), group_identifier
                )
                container.add_attribute(
                    const.UNIT, const.MZ, group_identifier
                )
            else:
                return False
    return True


@msp_interpretation_attribute_handler.add
@FunctionAttributeHandler.wraps("Dev_ppm")
def dev_ppm_handler(key, value, container) -> bool:
    if value is None:
        return False
    group_identifier = container.get_next_group_identifier()
    container.add_attribute(const.DELTA_MZ, try_cast(value), group_identifier)
    container.add_attribute(
        const.UNIT, const.PPM, group_identifier
    )
    return True


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("Nrep", "Nreps")
def nreps_handler(key, value, container):
    if value is None:
        return False
    if not isinstance(value, str):
        value = str(value)
    match = re.match(r"(\d+)/(\d+)", value)
    if match is not None:
        container.add_attribute(
            const.NUMBER_OF_REPLICATE_SPECTRA_USED, try_cast(match.group(1))
        )
        container.add_attribute(
            const.NUMBER_OF_REPLICATE_SPECTRA_AVAILABLE, try_cast(match.group(2))
        )
        return True
    else:
        match = re.match(r"(\d+)", value)
        if match is not None:
            container.add_attribute(
                const.NUMBER_OF_REPLICATE_SPECTRA_USED, try_cast(match.group(1))
            )
            container.add_attribute(
                const.NUMBER_OF_REPLICATE_SPECTRA_AVAILABLE,
                try_cast(match.group(1)),
            )
            return True
        else:
            return False


@msp_analyte_attribute_handler.add
@FunctionAttributeHandler.wraps("Organism")
def organism_handler(key, value, container):
    if value is None:
        return False
    value = value.strip('"')

    if value in species_map:
        group_identifier = container.get_next_group_identifier()
        for item in species_map[value]:
            container.add_attribute(item[0], try_cast(item[1]), group_identifier)
        return True
    return False


@msp_analyte_attribute_handler.add
@FunctionAttributeHandler.wraps("Protein")
def protein_handler(key, value, container: Attributed):
    if value is None:
        return False
    key = "MS:1000885|protein accession"
    match = re.match(r"\(pre=(.),post=(.)\)", value)
    group_identifier = None
    if match is not None:
        value = value[: match.start()]
        group_identifier = container.get_next_group_identifier()
        container.add_attribute(
            "MS:1001112|n-terminal flanking residue",
            match.group(1),
            group_identifier=group_identifier,
        )
        container.add_attribute(
            "MS:1001113|c-terminal flanking residue",
            match.group(2),
            group_identifier=group_identifier,
        )
    container.add_attribute(
        key,
        re.sub(r"\(pre=(.),post=(.)\)", "", value.strip('"').strip("'")),
        group_identifier=group_identifier,
    )
    return True


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("BasePeak")
def base_peak_handler(key, value, container: Attributed):
    if value is None:
        return False
    value = float(value)
    group_id = container.get_next_group_identifier()
    container.add_attribute("MS:1000505|base peak intensity", value, group_id)
    container.add_attribute(
        const.UNIT, "MS:1000131|number of detector counts", group_id
    )
    return True


class _UnknownTermTracker:
    counts: DefaultDict

    def add(self, key: str, value: Optional[str] = None):
        """Add an unknown attribute to the tracker"""
        raise NotImplementedError()

    def items(self):
        return self.counts.items()


class UnknownKeyValueTracker(_UnknownTermTracker):
    """
    A diagnostic tool for tracking attributes with values that the parser doesn't know how to interpret.

    This tracker holds both keys and values, and can grow quite large. For debugging purposes only.
    """

    def __init__(self) -> None:
        self.counts = DefaultDict(lambda: DefaultDict(int))

    def add(self, key: str, value: Optional[str] = None):
        self.counts[key][value] += 1


class UnknownKeyTracker(_UnknownTermTracker):
    """A diagnostic tool for tracking attributes that the parser doesn't know how to interpret."""

    def __init__(self) -> None:
        self.counts = DefaultDict(int)

    def add(self, key: str, value: Optional[str] = None):
        self.counts[key] += 1


protein_attributes_to_group = [
    "MS:1003048|number of enzymatic termini",
    "MS:1001045|cleavage agent name",
    "MS:1001112|n-terminal flanking residue",
    "MS:1001113|c-terminal flanking residue",
    "MS:1003044|number of missed cleavages",
    "MS:1000885|protein accession",
]


@dataclass
class AttributeRuleBase(_AttributeHelperBase):
    engine: Optional["AttributeRulesEngine"] = field(default=None, init=False)

    def is_applicable(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> bool:
        raise NotImplementedError()

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
        unknown_terms: Optional[List[str]] = None,
    ) -> List[str]:
        raise NotImplementedError()


class StrippedPeptideParserRule(AttributeRuleBase):
    def is_applicable(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> bool:
        return not analyte.has_attribute(
            STRIPPED_PEPTIDE_TERM
        ) and spectrum.has_attribute(SPECTRUM_NAME)

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
        unknown_terms: Optional[List[str]] = None,
    ) -> List[str]:
        name = spectrum.get_attribute(SPECTRUM_NAME)
        if name:
            match = re.match(r"([ARNDCEQGHKMFPSTWYVIL]+)/(\d+)", name)
            if match:
                self.add_value(STRIPPED_PEPTIDE_TERM, match.group(1), analyte)
                self.add_value(CHARGE_STATE, try_cast(match.group(2)), analyte)
        []


@dataclass
class FullNameParsingRule(AttributeRuleBase):
    attribute: str = "Fullname"

    def is_applicable(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> bool:
        return self.attribute in attributes

    def parse_value(self, value: str):
        match = re.match(
            r"([A-Z\-\*])\.([A-Z]+)\.([A-Z\-\*])/*([\d]*)",
            value,
        )
        if match is None:
            return None
        peptide = match.group(2)
        nterm = match.group(1)
        cterm = match.group(3)
        charge = match.group(4)
        return peptide, nterm, cterm, charge

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
        unknown_terms: Optional[List[str]] = None,
    ) -> List[str]:
        value = attributes[self.attribute]
        more_unknown_terms = []
        if value is not None:
            match = self.parse_value(value)
            if match is not None:
                (peptide, nterm, cterm, charge) = match
                analyte.add_attribute(STRIPPED_PEPTIDE_TERM, peptide)
                analyte.add_attribute("MS:1001112|n-terminal flanking residue", nterm)
                analyte.add_attribute("MS:1001113|c-terminal flanking residue", cterm)
                if charge:
                    analyte.add_attribute(CHARGE_STATE, try_cast(charge))
                attributes.pop(self.attribute)
            else:
                more_unknown_terms.append(self.attribute)
        else:
            more_unknown_terms.append(self.attribute)
        return more_unknown_terms


@dataclass
class PeptideNotationUnpackingRule(AttributeRuleBase):
    modification_parser: ModificationParser

    def is_applicable(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> bool:
        return bool(analyte.attributes)

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
        unknown_terms: Optional[List[str]] = None,
    ) -> List[str]:
        self.complete_analyte(analyte, spectrum)
        return []

    def complete_analyte(self, analyte: Analyte, spectrum: Spectrum):
        if analyte.has_attribute(STRIPPED_PEPTIDE_TERM):
            peptide = proforma.ProForma.parse(
                analyte.get_attribute(STRIPPED_PEPTIDE_TERM)
            )

            if analyte.has_attribute(PEPTIDE_MODIFICATION_TERM):
                modification_details = analyte.get_attribute(PEPTIDE_MODIFICATION_TERM)
                mods = self.modification_parser(modification_details)
                for position, residue, mod in mods:
                    if position == 0 and mod in TERMINAL_MODIFICATIONS:
                        peptide.n_term = [proforma.GenericModification(mod)]
                    else:
                        seqpos = list(peptide.sequence[position])
                        if not seqpos[1]:
                            seqpos[1] = [proforma.GenericModification(mod)]
                        peptide.sequence[position] = tuple(seqpos)
                        assert seqpos[0] == residue

            charge = None
            if analyte.charge:
                charge = analyte.charge
            elif spectrum.precursor_charge:
                charge = spectrum.precursor_charge

            if charge is not None:
                peptide.charge_state = charge
                analyte.add_attribute(PROFORMA_ION, str(peptide))
            else:
                analyte.add_attribute(PROFORMA_SEQ, str(peptide))

            if analyte.has_attribute(PEPTIDE_MODIFICATION_TERM):
                analyte.remove_attribute(PEPTIDE_MODIFICATION_TERM)
            analyte.add_attribute(THEORETICAL_MASS, peptide.mass)

        self.pack_protein_description(analyte)

        self.deduplicate_attribute(analyte, MOLECULAR_MASS)
        self.deduplicate_attribute(analyte, MOLECULAR_FORMULA)

    def pack_protein_description(self, analyte: Analyte):
        table = {}
        max_len = 0

        # Collect terms that describe a protein
        for term in protein_attributes_to_group:
            if not analyte.has_attribute(term):
                table[term] = []
                continue
            values = analyte.get_attribute(term, raw=True)
            if not isinstance(values, list):
                values = [values]
            max_len = max(max_len, len(values))
            table[term] = values

        # Ensure that all arrays are the same length
        for k, v in table.items():
            if len(v) < max_len:
                v.extend([None] * (max_len - len(v)))

        # Group together terms and remove the previous entries
        groups = []
        for i in range(max_len):
            group = []
            for k, v in table.items():
                if v[i] is None:
                    continue
                group.append((k, v[i].value))
                try:
                    analyte.remove_attribute(v[i].key, v[i].group_id)
                except KeyError:
                    # Sometimes, there are multiple entries for the same attribute,
                    # some with a group ID, some without. If an entry without a group
                    # ID is removed, every entry is removed, and subsequent removals
                    # do not find the entry.
                    pass
            groups.append(group)

        # Now add back the groups
        for group in groups:
            analyte.add_attribute_group(group)


@dataclass
class AnalyteMassDescriptorDiscriminationRule(AttributeRuleBase):

    SELECTED_ION_MZ_TERM = "MS:1000744|selected ion m/z"
    ADDUCT_ION_MASS_TERM = "MS:1003243|adduct ion mass"
    MOLECULAR_MASS_TERM = "MS:1000224|molecular mass"
    EXP_PREC_MONO_MZ_TERM = "MS:1003208|experimental precursor monoisotopic m/z"
    THEO_AVG_MZ_TERM = "MS:1003054|theoretical average m/z"

    precursor_ion_mz_terms = {
        "PrecursorMonoisoMZ": EXP_PREC_MONO_MZ_TERM,
        "ObservedPrecursorMZ": EXP_PREC_MONO_MZ_TERM,
        "PrecursorMZ": EXP_PREC_MONO_MZ_TERM,
        "PRECURSORMZ": EXP_PREC_MONO_MZ_TERM,
        "precursor": EXP_PREC_MONO_MZ_TERM,
        "precursor_mass": EXP_PREC_MONO_MZ_TERM,
        "precursormass": EXP_PREC_MONO_MZ_TERM,
        "Mz_exact": EXP_PREC_MONO_MZ_TERM,
        "Mz_av": THEO_AVG_MZ_TERM,
        "Parent": SELECTED_ION_MZ_TERM,
    }

    molecular_mass_terms = {
        "MW": ADDUCT_ION_MASS_TERM,
        "total exact mass": MOLECULAR_MASS_TERM,
        "ExactMass": MOLECULAR_MASS_TERM,
        "exact_mass": MOLECULAR_MASS_TERM,
    }

    @dataclass
    class _CommentAttrValue:
        comment: str
        attribute: str
        value: float
        error: Optional[float] = None

    def is_applicable(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ) -> bool:
        attributes = CaseInsensitiveDict(attributes)
        for key in itertools.chain.from_iterable(
            (self.molecular_mass_terms, self.precursor_ion_mz_terms)
        ):
            if key in attributes:
                return True
        return False

    def collect_mz_terms(
        self, attributes: Dict[str, Any]
    ) -> List["AnalyteMassDescriptorDiscriminationRule._CommentAttrValue"]:
        terms = []
        for k, v in self.precursor_ion_mz_terms.items():
            if k in attributes:
                terms.append(self._CommentAttrValue(k, v, try_cast(attributes[k])))
        return terms

    def collect_mass_terms(
        self, attributes: Dict[str, Any]
    ) -> List["AnalyteMassDescriptorDiscriminationRule._CommentAttrValue"]:
        terms = []
        for k, v in self.molecular_mass_terms.items():
            if k in attributes:
                terms.append(self._CommentAttrValue(k, v, try_cast(attributes[k])))
        return terms

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
        unknown_terms: Optional[List[str]] = None,
    ) -> List[str]:
        peptide = analyte.peptide

        adduct = PROTON
        molecular_mass_computed = None
        adduct_mass_computed = None

        molecular_mass_candidates: List[
            AnalyteMassDescriptorDiscriminationRule._CommentAttrValue
        ] = []
        adduct_mass_candidates: List[
            AnalyteMassDescriptorDiscriminationRule._CommentAttrValue
        ] = []

        explained_attributes = []
        # This entry has a peptide analyte
        if peptide is not None:
            peptide_mass = peptide.mass
            analyte_charge = spectrum.precursor_charge

            mass_terms = self.collect_mass_terms(attributes)

            if analyte_charge is None:
                analyte_charge = spectrum.precursor_charge

            if analyte_charge is not None:
                molecular_mass_computed = peptide_mass
                adduct_mass_computed = peptide_mass + adduct * analyte_charge

            categories = [
                (
                    self.MOLECULAR_MASS_TERM,
                    molecular_mass_computed,
                ),
                (
                    self.ADDUCT_ION_MASS_TERM,
                    adduct_mass_computed,
                ),
            ]

            for attr in mass_terms:
                best_category = None
                best_error = float("inf")
                for kind, theoretical in categories:
                    err = abs(attr.value - theoretical)
                    if err < best_error:
                        best_category = kind
                        best_error = err

                attr.error = best_error
                if best_error > 0.5 and math.isfinite(best_error):
                    logger.warning(
                        "Large error detected in best mass category classification: %s %r",
                        best_category,
                        best_error,
                    )
                if best_category == self.MOLECULAR_MASS_TERM:
                    molecular_mass_candidates.append(attr)
                elif best_category == self.ADDUCT_ION_MASS_TERM:
                    adduct_mass_candidates.append(attr)
                else:
                    logger.warning(
                        "No mass terms detected for spectrum %s", spectrum.key
                    )

            if adduct_mass_candidates:
                adduct_mass_candidates.sort(key=lambda x: x.error)
                if len(adduct_mass_candidates) > 1:
                    warnings.warn(
                        "Multiple adduct mass candidates found: %r",
                        adduct_mass_candidates,
                    )
                term = adduct_mass_candidates[0]
                explained_attributes.append(term.comment)
                analyte.add_attribute(self.ADDUCT_ION_MASS_TERM, term.value)
            if molecular_mass_candidates:
                molecular_mass_candidates.sort(key=lambda x: x.error)
                if len(molecular_mass_candidates) > 1:
                    warnings.warn(
                        "Multiple molecular mass candidates found: %r",
                        molecular_mass_candidates,
                    )
                term = adduct_mass_candidates[0]
                explained_attributes.append(term.comment)
                analyte.add_attribute(self.MOLECULAR_MASS_TERM, term.value)

        # This entry has no molecular analyte or we don't know how to get a mass from it
        else:
            mass_terms = self.collect_mass_terms(attributes)
            for term in mass_terms:
                analyte.add_attribute(term.attribute, term.value)
                explained_attributes.append(term.comment)

        if unknown_terms:
            for term in explained_attributes:
                unknown_terms.remove(term)


@dataclass
class AttributeRulesEngine:
    other_manager: AttributeHandler
    analyte_manager: AttributeHandler
    interpretation_manager: AttributeHandler
    interpretation_member_manager: AttributeHandler
    spectrum_attribute_handler: AttributeHandler
    analyte_attribute_handler: AttributeHandler

    early_rules: List[AttributeRuleBase] = field(default_factory=list)
    late_rules: List[AttributeRuleBase] = field(default_factory=list)

    def bind_rule(self, rule: AttributeRuleBase):
        rule.engine = self
        return rule

    def __post_init__(self):
        for rule in self.early_rules:
            self.bind_rule(rule)
        for rule in self.late_rules:
            self.bind_rule(rule)

    def process_attributes(
        self,
        attributes: Dict[str, Any],
        spectrum: Spectrum,
        analyte: Analyte,
        interpretation: Interpretation,
        interpretation_member: InterpretationMember,
    ):
        unknown_terms = []

        for early_rule in self.early_rules:
            if early_rule.is_applicable(
                attributes, spectrum, analyte, interpretation, interpretation_member
            ):
                early_rule.process_attributes(
                    attributes, spectrum, analyte, interpretation, interpretation_member
                )

        for attribute in attributes:
            # Skip a leader term that we already processed
            if attribute in leader_terms:
                continue
            if not attribute:
                continue
            if attribute in self.other_manager:
                if not self.other_manager(attribute, attributes[attribute], spectrum):
                    unknown_terms.append(attribute)

            elif attribute in self.analyte_manager:
                if not self.analyte_manager(attribute, attributes[attribute], analyte):
                    unknown_terms.append(attribute)

            elif attribute in self.interpretation_manager:
                if not self.interpretation_manager(
                    attribute, attributes[attribute], interpretation
                ):
                    unknown_terms.append(attribute)

            elif attribute in self.interpretation_member_manager:
                if not self.interpretation_member_manager(
                    attribute, attributes[attribute], interpretation_member
                ):
                    unknown_terms.append(attribute)

            elif attribute in self.spectrum_attribute_handler:
                if not self.spectrum_attribute_handler(
                    attribute, attributes[attribute], spectrum
                ):
                    unknown_terms.append(attribute)

            elif attribute in self.analyte_attribute_handler:
                if not self.analyte_attribute_handler(
                    attribute, attributes[attribute], analyte
                ):
                    unknown_terms.append(attribute)
            else:
                unknown_terms.append(attribute)

        for late_rule in self.late_rules:
            if late_rule.is_applicable(
                attributes, spectrum, analyte, interpretation, interpretation_member
            ):
                late_rule.process_attributes(
                    attributes,
                    spectrum,
                    analyte,
                    interpretation,
                    interpretation_member,
                    unknown_terms=unknown_terms,
                )

        return unknown_terms


def proportion_parser(aggregation: str) -> float:
    """Parse for fractions or percentages"""
    if "/" in aggregation:
        aggregation: Fraction = Fraction(aggregation)
        return float(aggregation)
    else:
        return float(aggregation)


@dataclass
class PeakAggregationParser:
    """
    Parse peak aggregation information.

    Subtypes may produce different attributes.

    Attributes
    ----------
    peak_attributes : list[(:class:`~.Attribute`, :class:`PeakAggregateParseFn`)]
        The attributes this parser expects
    """

    peak_attributes: List[Tuple[Attribute, PeakAggregateParseFn]] = field(
        default_factory=lambda: [
            (Attribute(PEAK_ATTRIB, PEAK_OBSERVATION_FREQ), proportion_parser)
        ]
    )

    def __call__(
        self, aggregation: str, wrap_errors: bool = True, **kwargs
    ) -> List[Tuple[Attribute, float]]:
        parsed = []

        if isinstance(aggregation, str):
            aggregation = SPACE_SPLITTER.split(aggregation)

        for (i, token), (k, parser) in zip(
            enumerate(aggregation), self.peak_attributes
        ):
            try:
                result = parser(token)
                parsed.append((k, result))
            except Exception as err:
                if not wrap_errors:
                    raise err from None
                else:
                    logger.error(f"Failed to parse aggregation at {i}")
                    parsed.append((k, token))
        return parsed


@dataclass
class PeakParsingStrategy:
    """
    A combination of peak annotation parsing and peak aggregation parsing
    strategies.

    Attributes
    ----------
    annotation_parser : :class:`MSPAnnotationStringParser` or :const:`None`
        The peak annotation parser for this format.

    aggregation_parser : :class:`PeakAggregationParser` or :const:`None`
        The peak aggregation parser for this format
    """

    annotation_parser: Optional[MSPAnnotationStringParser] = field(
        default=parse_annotation
    )
    aggregation_parser: Optional[PeakAggregationParser] = field(
        default_factory=PeakAggregationParser
    )

    def has_aggregation(self) -> bool:
        return self.aggregation_parser is not None

    def has_annotation(self) -> bool:
        return self.annotation_parser is not None

    def parse_aggregation(self, aggregation: str, wrap_errors: bool = True, **kwargs):
        return self.aggregation_parser(
            aggregation=aggregation, wrap_errors=wrap_errors, **kwargs
        )

    def parse_annotation(self, annotation: str, wrap_errors: bool = True, **kwargs):
        return self.annotation_parser(
            annotation_string=annotation, wrap_errors=wrap_errors, **kwargs
        )

    def process_peak_list(self, spectrum: Spectrum):
        aggregation_metrics_used = DefaultDict(int)
        for i, peak in enumerate(spectrum.peak_list):
            try:
                parsed_interpretation = self.parse_annotation(
                    peak[2], spectrum=spectrum
                )
            except ValueError as err:
                message = str(err)
                raise ValueError(
                    f"An error occurred while parsing the peak annotation for peak {i}: {message}"
                ) from err

            if (
                parsed_interpretation
                and isinstance(parsed_interpretation[0], annotation.InvalidAnnotation)
                and peak[3]
            ):
                logger.debug(
                    "Failed to parse interpretation string %r with assumed %d aggregation fields, trying to parse the combined string",
                    peak[2],
                    len(peak[3]),
                )

                peak[2] = " ".join([peak[2]] + peak[3])
                peak[3] = []
                try:
                    parsed_interpretation = self.parse_annotation(
                        peak[2], spectrum=spectrum
                    )
                except ValueError as err:
                    message = str(err)
                    raise ValueError(
                        f"An error occurred while parsing the peak annotation for peak {i}: {message}"
                    ) from err

            peak[2] = parsed_interpretation
            if peak[3]:
                aggregations = []
                for agg_type, agg in self.parse_aggregation(peak[3]):
                    aggregation_metrics_used[agg_type] += 1
                    aggregations.append(agg)

                peak[3] = aggregations

        aggregation_metrics = [
            metric for metric, count in aggregation_metrics_used.items() if count > 0
        ]

        if aggregation_metrics:
            spectrum.add_attribute_group(aggregation_metrics)

    def parse_peak_list(self, peak_lines: Iterable[str]) -> List[RawPeakLine]:
        peak_list = []
        for values in peak_lines:
            interpretations = ""
            aggregation = ""
            if len(values) == 1:
                mz = values
                intensity = "1"
            if len(values) == 2:
                mz, intensity = values
            elif len(values) == 3:
                mz, intensity, interpretations = values
            elif len(values) > 3:
                mz, intensity, interpretations = values[0:2]
            else:
                mz = "1"
                intensity = "1"

            interpretations = interpretations.strip('"')
            aggregation = None
            if interpretations.startswith("?"):
                parts = SPACE_SPLITTER.split(interpretations)
                if len(parts) > 1:
                    # Some msp files have a concept for ?i, but this requires a definition
                    interpretations = "?"
                    aggregation = parts[1:]
            else:
                if " " in interpretations:
                    parts = SPACE_SPLITTER.split(interpretations)
                    interpretations = parts[0]
                    aggregation = parts[1:]

            #### Add to the peak list
            peak_list.append(
                [float(mz), float(intensity), interpretations, aggregation]
            )
        return peak_list

    def __call__(self, peak_lines: Iterable[str]) -> list:
        raw_peaks = self.parse_peak_list(peak_lines)
        return raw_peaks


class MSPSpectralLibrary(_PlainTextSpectralLibraryBackendBase):
    """
    A reader for the plain text NIST MSP spectral library format.

    The MSP format is only roughly defined, and does places few
    constraints on the meanings of spectrum attributes. This parser
    attempts to cover a variety of different ways that MSPs found
    "in the wild" have denoted different spectrum properties, but
    is neither exhaustive nor nuanced enough to know from context
    exactly what those files' authors intended, making a best guess
    at when they correspond to in the controlled vocabulary mapping
    for :mod:`mzspeclib`

    Attributes
    ----------
    modification_parser : :class:`ModificationParser`
        A parser for peptide modifications
    unknown_attributes : :class:`_UnknownTermTracker`
        A tracker for unknown attributes. Used to tell how much information
        the reader is unable to map onto the controlled vocabulary.
    """

    file_format = "msp"
    format_name = "msp"

    modification_parser: ModificationParser
    unknown_attributes: _UnknownTermTracker

    def __init__(
        self, filename, index_type=None, read_metadata=True, create_index: bool = True
    ):
        super().__init__(filename, index_type, read_metadata, create_index=create_index)
        self.modification_parser = ModificationParser()
        self.unknown_attributes = UnknownKeyTracker()

    @classmethod
    def guess_from_header(cls, filename: str) -> bool:
        with open_stream(filename, "r") as stream:
            first_line = stream.readline()
            if LEADER_TERMS_PATTERN.match(first_line):
                return True
        return False

    def read_header(self) -> bool:
        filename = self.filename
        file_like_object = isinstance(filename, io.IOBase)

        stream = open_stream(filename, "rt")
        match, offset = self._parse_header_from_stream(stream)
        if file_like_object:
            if stream.seekable():
                stream.seek(0)
            stream.detach()
        else:
            stream.close()
        return match

    def _parse_header_from_stream(self, stream: io.IOBase) -> Tuple[bool, int]:
        first_line = stream.readline()
        if isinstance(first_line, bytes):
            first_line = first_line.decode("utf8")
        attributes = AttributeManager()
        attributes.add_attribute(FORMAT_VERSION, DEFAULT_VERSION)
        if isinstance(self.filename, (str, os.PathLike)):
            attributes.add_attribute(
                LIBRARY_NAME, self.filename.rsplit(".msp", 1)[0].split(os.sep)[-1]
            )
        elif hasattr(stream, "name"):
            attributes.add_attribute(
                LIBRARY_NAME, stream.name.rsplit(".msp", 1)[0].split(os.sep)[-1]
            )
        self.attributes.clear()
        self.attributes._from_iterable(attributes)
        if LEADER_TERMS_PATTERN.match(first_line):
            return True, 0
        return False, 0

    def create_index(self) -> int:
        """
        Populate the spectrum index

        Returns
        -------
        n_spectra: int
            The number of entries read
        """
        #### Check that the spectrum library filename isvalid
        filename = self.filename

        begin_loc = None
        file_like_object = False
        #### Determine the filesize
        try:
            file_size = os.path.getsize(filename)
        except TypeError:
            if isinstance(filename, io.IOBase):
                file_like_object = True
                if filename.seekable():
                    begin_loc = filename.tell()
                    filename.seek(0, os.SEEK_END)
                    file_size = filename.tell()
                    filename.seek(0)
        infile = open_stream(filename, "r")

        state = "header"
        spectrum_buffer = []
        n_spectra = 0
        start_index = 0
        file_offset = 0
        line_beginning_file_offset = 0
        spectrum_file_offset = 0
        spectrum_name = ""

        # Required for counting file_offset manually (LF vs CRLF)
        infile.readline()
        file_offset_line_ending = len(infile.newlines) - 1
        infile.seek(0)

        # if debug:
        #     eprint("INFO: Reading..", end='', flush=True)
        logger.debug(f"Reading {filename} ({file_size} bytes)...")
        while 1:
            line = infile.readline()
            if len(line) == 0:
                break

            line_beginning_file_offset = file_offset

            #### tell() is twice as slow as counting it myself
            # file_offset = infile.tell()
            file_offset += len(line) + file_offset_line_ending

            line = line.rstrip()
            # TODO: Name: could be Compound or SpectrumName
            if state == "header":
                if LEADER_TERMS_PATTERN.match(line):
                    state = "body"
                    spectrum_file_offset = line_beginning_file_offset
                else:
                    continue
            if state == "body":
                if len(line) == 0:
                    continue
                if LEADER_TERMS_PATTERN.match(line):
                    if len(spectrum_buffer) > 0:
                        self.index.add(
                            number=n_spectra + start_index,
                            offset=spectrum_file_offset,
                            name=spectrum_name,
                            analyte=None,
                        )
                        n_spectra += 1
                        spectrum_buffer = []
                        #### Commit every now and then
                        if n_spectra % 10000 == 0:
                            self.index.commit()
                            logger.info(
                                f"... Indexed  {file_offset} bytes, {n_spectra} spectra read"
                            )

                    spectrum_file_offset = line_beginning_file_offset
                    spectrum_name = LEADER_TERMS_LINE_PATTERN.match(line).group(1)

                spectrum_buffer.append(line)
        logger.debug(f"Processed {file_offset} bytes, {n_spectra} spectra read")
        self.index.add(
            number=n_spectra + start_index,
            offset=spectrum_file_offset,
            name=spectrum_name,
            analyte=None,
        )
        self.index.commit()
        n_spectra += 1

        #### Flush the index
        self.index.commit()

        if not file_like_object:
            infile.close()
        else:
            infile.detach()
            if begin_loc is not None:
                filename.seek(begin_loc)
        return n_spectra

    def _buffer_from_stream(self, infile: Iterable[str]) -> List[str]:
        state = "body"
        spectrum_buffer = []

        file_offset = 0
        line_beginning_file_offset = 0
        for line in infile:
            line_beginning_file_offset = file_offset
            file_offset += len(line)
            line = line.rstrip()
            if state == "body":
                if len(line) == 0:
                    continue
                if LEADER_TERMS_PATTERN.match(line):
                    if len(spectrum_buffer) > 0:
                        return spectrum_buffer
                spectrum_buffer.append(line)
        return spectrum_buffer

    def _make_peak_parsing_strategy(self) -> PeakParsingStrategy:
        return PeakParsingStrategy(parse_annotation)

    def _parse(self, buffer: Iterable[str], spectrum_index: int = None) -> Spectrum:
        #### Start in the header section of the entry
        in_header = True

        #### Reset all spectrum properties in case there is already data here
        attributes = {}
        peak_list = []
        peak_line_parser: PeakParsingStrategy = self._make_peak_parsing_strategy()

        #### Loop through each line in the buffered list
        for line in buffer:
            #### If in the the header portion of the entry
            if in_header:
                key = value = None
                #### Extract the key,value pair by splitting on the *first* colon with optional whitespace
                match = re.match(
                    r"\s*#", line
                )  # Assume lines starting with # are comments
                if match:
                    continue
                elif line.count(":") > 0:
                    key, value = re.split(r":\s*", line, 1)
                    attributes[key] = value
                elif line.count("=") > 0:
                    key, value = re.split(r"=\s*", line, 1)
                    attributes[key] = value
                elif line.count("\t") > 0:
                    warnings.warn(f"Line {line!r} looks like a peak annotation?")
                    in_header = False
                else:
                    key = line
                    attributes[key] = value

                #### If the key is "Num peaks" then we're done with the header and peak list follows
                if key in NUM_PEAKS_KEYS:
                    in_header = False

                #### The "Comment" key requires special parsing
                if key == "Comment" or key == "Comments":
                    #### Remove it from attributes
                    del attributes[key]
                    self._parse_comment(value, attributes)

            #### Else in the peaks section. Parse the peaks.
            else:
                #### Split into the expected three values
                values = re.split(r"\s+", line, maxsplit=2)
                # Sometimes MSP files have multiple peaks on the same line, delimited by ';',
                # so we must potentially parse more than one peak per line.
                if values[1].endswith(";"):
                    buffered_peaks = []
                    buffered_peaks.append(values[:2])
                    buffered_peaks[0][1] = buffered_peaks[0][1].strip(";")
                    rest = values[0]
                    for block in re.split(r";\s?", rest):
                        if block:
                            buffered_peaks.append(re.split(r"\s+", block, maxsplit=2))
                    peak_list.extend(buffered_peaks)
                else:
                    peak_list.append(values)

        peak_list = peak_line_parser(peak_list)

        #### Now convert the format attributes to standard ones
        spectrum = self._make_spectrum(peak_list, attributes)
        if spectrum_index is not None:
            spectrum.index = spectrum_index
        else:
            spectrum.index = -1
        spectrum.key = spectrum.index + 1
        peak_line_parser.process_peak_list(spectrum)
        return spectrum

    def _parse_annotation(self, annotation: str, wrap_errors: bool = True, **kwargs):
        return parse_annotation(
            annotation_string=annotation, wrap_errors=wrap_errors, **kwargs
        )

    def _parse_comment(self, value: str, attributes: Attributed):
        comment_items = re.split(" ", value)

        #### Any spaces within quotes are then de-split
        fixed_comment_items = []
        quote_counter = 0
        new_item = ""
        for item in comment_items:
            if new_item > "":
                new_item = new_item + " "
            new_item = new_item + item
            n_quotes = new_item.count('"')
            if n_quotes % 2 == 0:
                fixed_comment_items.append(new_item)
                new_item = ""

        #### Try to split each item on the first = character and store in attributes
        for item in fixed_comment_items:
            #### If there is an = then split on the first one and store key and value
            if item.count("=") > 0:
                comment_key, comment_value = item.split("=", 1)
                cleaned_key = comment_key.strip('"')
                if len(cleaned_key) != len(comment_key):
                    comment_value = comment_value.strip('"')
                attributes[cleaned_key] = try_cast(comment_value)

            #### Otherwise just store the key with a null value
            else:
                attributes[item] = None

    def _make_attribute_handlers(self) -> AttributeRulesEngine:
        """
        Create the attribute handling scopes that map this flavor of MSP's
        attributes onto controlled vocabulary terms in context.

        This method should be overridden in sub-classes to allow them
        to change the meanings of attributes, add new ones, or otherwise
        redirect how they are interpreted.

        See the :class:`AttributeHandler` type tree for more details about
        how the distributed predicates are resolved.

        Returns
        -------
        other_manager : :class:`AttributeHandler`
            The attribute handler for uncategorized attributes that will be added
            to a :class:`Spectrum`.
        analyte_manager : :class:`AttributeHandler`
            The attribute handler for attributes that will be added to a :class:`Analyte`
        interpretation_manager : :class:`AttributeHandler`
            The attribute handler for attributes that will be added to a :class:`Interpretation`
        interpretation_member_manager : :class:`AttributeHandler`
            The attribute handler for attributes that will be added to a :class:`InterpretationMember`
        spectrum_manager : :class:`AttributeHandler`
            The attribute handler for attributes that will be added to a :class:`Spectrum`
        analyte_fallback_manager : :class:`AttributeHandler`
            The attribute handler for attributes that will be tried for any attribute
            that fails to be categorized by all of the other managers to be added to the
            :class:`Analyte` before labeling the attribute as "unknown".
        """
        other_manager = MappingAttributeHandler(other_terms)
        analyte_manager = MappingAttributeHandler(analyte_terms)
        interpretation_manager = msp_interpretation_attribute_handler
        interpretation_member_manager = MappingAttributeHandler(
            interpretation_member_terms
        )

        return AttributeRulesEngine(
            other_manager,
            analyte_manager,
            interpretation_manager,
            interpretation_member_manager,
            msp_spectrum_attribute_handler,
            msp_analyte_attribute_handler,
            early_rules=[FullNameParsingRule()],
            late_rules=[
                StrippedPeptideParserRule(),
                PeptideNotationUnpackingRule(self.modification_parser),
                AnalyteMassDescriptorDiscriminationRule(),
            ],
        )

    def _make_spectrum(self, peak_list: List, attributes: Mapping[str, str]):
        spectrum = self._new_spectrum()
        interpretation = self._new_interpretation(FIRST_INTERPRETATION_KEY)
        interpretation_member = self._new_interpretation_member(1)
        analyte = self._new_analyte(FIRST_ANALYTE_KEY)

        # Assume MSP spectra always have an Analyte, we can remove it later.
        # It is less common to assume that they will have Interpretation
        # or InterpretationMember component, so we will add these only if
        # we find attributes for them
        spectrum.add_analyte(analyte)
        spectrum.peak_list = peak_list

        #### Add special terms that we want to start off with
        for term in leader_terms:
            if term in attributes:
                spectrum.add_attribute(leader_terms[term], try_cast(attributes[term]))
                break
        else:
            logger.error(
                "Did not find any MSP leader terms (%s) for Spectrum Index = %s",
                leader_terms,
                spectrum.index,
            )

        #### Translate the rest of the known attributes and collect unknown ones
        attribute_rules = self._make_attribute_handlers()
        unknown_terms = attribute_rules.process_attributes(
            attributes, spectrum, analyte, interpretation, interpretation_member
        )

        #### Handle the uninterpretable terms
        for attribute in unknown_terms:
            self.unknown_attributes.add(attribute, attributes[attribute])
            if attributes[attribute] is None:
                spectrum.add_attribute(
                    CUSTOM_ATTRIBUTE_NAME, try_cast(attribute)
                )
            else:
                group_identifier = spectrum.get_next_group_identifier()
                spectrum.add_attribute(
                    CUSTOM_ATTRIBUTE_NAME,
                    try_cast(attribute),
                    group_identifier,
                )
                spectrum.add_attribute(
                    CUSTOM_ATTRIBUTE_VALUE,
                    try_cast(attributes[attribute]),
                    group_identifier,
                )

        if interpretation_member:
            interpretation.add_member_interpretation(interpretation_member)

        if not self._is_analyte_defined(analyte):
            self._hoist_analyte_attributes_on_rejection(analyte, spectrum)
            analyte.clear()
            spectrum.remove_analyte(analyte.id)
        elif analyte:
            spectrum.add_analyte(analyte)
            interpretation.add_analyte(analyte)
            spectrum.add_interpretation(interpretation)
        return spectrum

    def get_spectrum(
        self, spectrum_number: int = None, spectrum_name: str = None
    ) -> Spectrum:
        # keep the two branches separate for the possibility that this is not possible with all
        # index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError("Provide only one of spectrum_number or spectrum_name")
            index_record = self.index.record_for(spectrum_number)
            offset = index_record.offset
        elif spectrum_name is not None:
            index_record = self.index.record_for(spectrum_name)
            spectrum_number = index_record.number
            offset = index_record.offset
        else:
            raise ValueError(
                "Must provide either spectrum_number or spectrum_name argument"
            )
        buffer = self._get_lines_for(offset)
        spectrum = self._parse(buffer, index_record.index)
        return spectrum

    def summarize_parsing_errors(self) -> Dict:
        errors = super().summarize_parsing_errors()
        errors["unknown attributes"] = DefaultDict(dict)
        for k, v in self.unknown_attributes.items():
            errors["unknown attributes"][k] = v
        return errors


def _parse_fraction(x: str) -> float:
    a, b = x.split("/")
    return int(a) / int(b)


class MSPSpectralLibraryWriter(SpectralLibraryWriterBase):
    file_format = "msp"

    analyte_keys = {
        MOLECULAR_FORMULA: "Formula",
        "MS:1003044|number of missed cleavages": "MC",
        "MS:1001471|peptide modification details": "Mods",
        "MS:1003043|number of residues": "Naa",
        PRECURSOR_MZ: "PrecursorMonoisoMZ",
        "MS:1003054|theoretical average m/z": "Mz_av",
        PROFORMA_SEQ: "ProForma",
        STRIPPED_PEPTIDE_TERM: "Peptide",
        CHARGE_STATE: "Charge",
    }

    for species_name, keys in species_map.items():
        analyte_keys[tuple(keys[0])] = ("Organism", species_name)

    modification_map = {v: k for k, v in MODIFICATION_NAME_MAP.items()}

    spectrum_keys = {
        CHARGE_STATE: "Charge",
        (
            SPECTRUM_AGGREGATION_TYPE,
            SINGLETON_SPECTRUM,
        ): "Single",
        (
            SPECTRUM_AGGREGATION_TYPE,
            CONSENSUS_SPECTRUM,
        ): "Consensus",
        SCAN_NUMBER: "Scan",
        SOURCE_FILE: "Origfile",
        const.SAMPLE_NAME: "Sample",
        const.FILTER_STRING: "Filter",
        const.PRECURSOR_APEX_INTENSITY: "Precursor1MaxAb",
        const.ISOLATION_WINDOW_PRECURSOR_PURITY: "Purity",
        const.BASE_PEAK_INTENSITY: "BasePeak",
        "MS:1002599|splash key": "Splash",
        "MS:1003289|intensity of highest unassigned peak": "max_unassigned_ab",
        TOP_20_UNASSIGNED_INTENSITY_FRACTION: "Unassigned",
        TOTAL_UNASSIGNED_INTENSITY_FRACTION: "Unassign_all",
        NUM_UNASSIGNED_PEAKS_IN_TOP_20: "top_20_num_unassigned_peaks",
        (
            "MS:1000044|dissociation method",
            "MS:1000422|beam-type collision-induced dissociation",
        ): "HCD",
        (
            "MS:1000044|dissociation method",
            "MS:1002472|trap-type collision-induced dissociation",
        ): "CID",
    }

    # TODO: add these
    interpretation_keys = {
        const.Q_VALUE: "Q-value",
    }

    def __init__(self, filename, **kwargs):
        super(MSPSpectralLibraryWriter, self).__init__(filename)
        self._coerce_handle(self.filename)

    def write_header(self, library: SpectralLibraryBackendBase):
        pass

    def write_attribute_set(
        self, attribute_set: AttributeSet, attribute_set_type: AttributeSetTypes
    ):
        pass

    def _format_value(self, value):
        if isinstance(value, str):
            if not (value.startswith('"') and value.endswith('"')):
                value = f'"{value}"'
        return str(value)

    def _proforma_to_mods(self, proforma_seq: str) -> str:
        if isinstance(proforma_seq, proforma.ProForma):
            parsed = proforma_seq
        else:
            parsed = proforma.ProForma.parse(proforma_seq)
        mods = [(i, tok) for i, tok in enumerate(parsed) if tok[1]]
        if mods:
            tokens = []
            for i, mod_site in mods:
                tokens.append(str(i))
                tokens.append(mod_site[0])
                mod = mod_site[1][0]
                mod_name = self.modification_map.get(mod.name, mod.name)
                tokens.append(mod_name)
            return f'Mods={len(mods)}({",".join(tokens)})'
        else:
            return "Mods=0"

    def _protein_to_comments(self, analyte: Analyte) -> List[str]:
        acc = []
        protein: ProteinDescription
        for protein in analyte.proteins:
            accession = None
            pre = None
            post = None
            if protein.accession:
                accession = protein.accession
            if protein.flanking_n_terminal_residue:
                pre = protein.flanking_n_terminal_residue
            if protein.flanking_c_terminal_residue:
                post = protein.flanking_c_terminal_residue
            if accession:
                token = accession
                if token.startswith('"'):
                    token = token.strip('"')
                if pre or post:
                    token += f"(pre={pre or '-'},post={post or '-'})"
                acc.append(f"Protein={self._format_value(token)}")
                if protein.number_of_enzymatic_termini == 2:
                    if protein.cleavage_agent == "MS:1001251|Trypsin":
                        acc.append("Pep=Tryptic")
                elif protein.number_of_enzymatic_termini == 1:
                    if protein.cleavage_agent == "MS:1001251|Trypsin":
                        acc.append("Pep=SemiTryptic")
                if protein.missed_cleavages is not None:
                    acc.append(f"MC={self._format_value(protein.missed_cleavages)}")
                break
        return acc

    def _build_comments(
        self, spectrum: Spectrum, attribute_container: Attributed, rule_map: Dict
    ) -> List[Tuple[str, str]]:
        accumulator = []

        for attr_name, msp_name in rule_map.items():
            if isinstance(attr_name, tuple):
                if attribute_container.has_attribute(attr_name[0]):
                    value = attribute_container.get_attribute(attr_name[0])
                    if value == attr_name[1]:
                        if isinstance(msp_name, str):
                            accumulator.append(msp_name)
                        elif isinstance(msp_name, (list, tuple)) and len(msp_name) == 2:
                            accumulator.append("=".join(msp_name))
                        else:
                            raise TypeError(
                                f"Can't infer conversion for {msp_name} given {attr_name}"
                            )
            elif attribute_container.has_attribute(attr_name):
                value = attribute_container.get_attribute(attr_name)
                if isinstance(value, list):
                    logger.warn(
                        "Spectrum %r contains multiple values for %r, only the first will be saved",
                        spectrum.name,
                        attr_name,
                    )
                    accumulator.append(f"{msp_name}={self._format_value(value[0])}")
                else:
                    accumulator.append(f"{msp_name}={self._format_value(value)}")
        return accumulator

    def build_spectrum_comments(self, spectrum: Spectrum) -> List[Tuple[str, str]]:
        accumulator = self._build_comments(spectrum, spectrum, self.spectrum_keys)
        if spectrum.analytes:
            analyte = spectrum.get_analyte("1")
            accumulator += self._build_comments(spectrum, analyte, self.analyte_keys)
            if analyte.peptide:
                accumulator.append(self._proforma_to_mods(analyte.peptide))
            accumulator += self._protein_to_comments(analyte)
        if spectrum.interpretations:
            interp = spectrum.get_interpretation("1")
            accumulator += self._build_comments(
                spectrum, interp, self.interpretation_keys
            )
        return accumulator

    def write_spectrum(self, spectrum: Spectrum):
        if len(spectrum.analytes) > 1:
            logger.warning(
                "Spectrum %r contains multiple analytes, MSP will only contain the first",
                spectrum.name,
            )
        analyte = spectrum.get_analyte("1")
        self.handle.write(f"Name: {spectrum.name}\n")
        self.handle.write(f"MW: {analyte.mass}\n")
        self.handle.write(
            f"Comment: {' '.join(self.build_spectrum_comments(spectrum))}\n"
        )
        self._write_peaks(spectrum)
        self.handle.write("\n")

    def _format_annotation(self, annot: annotation.IonAnnotationBase):
        parts = []
        if isinstance(annot, annotation.PeptideFragmentIonAnnotation):
            parts.append(f"{annot.series}{annot.position}")
        elif isinstance(annot, annotation.ImmoniumIonAnnotation):
            parts.append(
                f"I{annot.amino_acid}{annot.modification if annot.modification else ''}"
            )
        elif isinstance(annot, annotation.ReporterIonAnnotation):
            parts.append(annot.reporter_label)
        if not parts:
            return "?"
        if annot.neutral_losses:
            f = annotation.NeutralName.combine(annot.neutral_losses)
            if f[0] not in ("-", "+"):
                f = "+" + f
            parts.append(f)
        if annot.adducts:
            parts.append(f"[{annotation.combine_formula(annot.adducts)}]")
        if annot.charge > 1:
            parts.append(f"^{annot.charge}")
        if annot.isotope:
            if annot.isotope > 0:
                parts.append(f"+{annot.isotope}i")
            else:
                parts.append(f"{annot.isotope}i")
        if annot.mass_error:
            parts.append(f"/{annot.mass_error}")
        return "".join(parts)

    def _write_peaks(self, spectrum: Spectrum):
        self.handle.write(f"Num peaks: {len(spectrum.peak_list)}\n")
        for peak in spectrum.peak_list:
            annot = self._format_annotation(peak[2][0]) if peak[2] else "?"
            self.handle.write(f'{peak[0]:0.4f}\t{peak[1]:0.4f}\t"{annot}"\n')

    def close(self):
        if not self.handle.closed:
            self.handle.close()
