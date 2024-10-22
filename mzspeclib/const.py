"""A collection of constants"""

PROTON = 1.00727646677

FIRST_ANALYTE_KEY = "1"
FIRST_INTERPRETATION_KEY = "1"


LIBRARY_CREATION_SW = "MS:1003207|library creation software"
SW_VERSION = "MS:1003200|software version"

ATTRIBUTE_SET_NAME = "MS:1003212|library attribute set name"
PEAK_ATTRIBUTE = "MS:1003254|peak attribute"

ANALYTE_MIXTURE = "MS:1003163|analyte mixture members"
CHARGE_STATE = "MS:1000041|charge state"
PROFORMA_ION = "MS:1003270|proforma peptidoform ion notation"
PROFORMA_SEQ = "MS:1003169|proforma peptidoform sequence"
STRIPPED_PEPTIDE_SEQ = "MS:1000888|stripped peptide sequence"
MOLECULAR_FORMULA = "MS:1000866|molecular formula"
RETENTION_TIME = "MS:1000894|retention time"

MOLECULAR_MASS = "MS:1000224|molecular mass"
ADDUCT_ION_MASS = "MS:1003243|adduct ion mass"
THEORETICAL_MASS = "MS:1001117|theoretical mass"

FORMAT_VERSION = "MS:1003186|library format version"
LIBRARY_NAME = "MS:1003188|library name"
LIBRARY_VERSION = "MS:1003190|library version"
LIBRARY_IDENTIFIER = "MS:1003187|library identifier"
LIBRARY_DESCRIPTION = "MS:1003189|library description"
LIBRARY_URI = "MS:1003191|library URI"

SAMPLE_NAME = "MS:1000002|sample name"
SOURCE_FILE = "MS:1003203|constituent spectrum file"
SCAN_NUMBER = "MS:1003057|scan number"
MS_LEVEL = "MS:1000511|ms level"
FILTER_STRING = "MS:1000512|filter string"

SPECTRUM_NAME = "MS:1003061|library spectrum name"
LIBRARY_SPECTRUM_KEY = "MS:1003237|library spectrum key"
LIBRARY_SPECTRUM_INDEX = "MS:1003062|library spectrum index"

SELECTED_ION_MZ = "MS:1000744|selected ion m/z"
PRECURSOR_MZ = "MS:1003208|experimental precursor monoisotopic m/z"
THEORETICAL_MZ = "MS:1003053|theoretical monoisotopic m/z"
DELTA_MZ = "MS:1001975|delta m/z"
ISOLATION_WINDOW_PRECURSOR_PURITY = "MS:1009013|isolation window precursor purity"
PRECURSOR_APEX_INTENSITY = "MS:1003086|precursor apex intensity"
BASE_PEAK_INTENSITY = "MS:1000505|base peak intensity"

SPECTRUM_ORIGIN_TYPE = "MS:1003072|spectrum origin type"
SELECTED_FRAGMENT_THEORETICAL_OBSERVED_INTENSITY_SPECTRUM = "MS:1003424|selected fragment theoretical m/z observed intensity spectrum"

SPECTRUM_AGGREGATION_TYPE = "MS:1003065|spectrum aggregation type"
SINGLETON_SPECTRUM = "MS:1003066|singleton spectrum"
CONSENSUS_SPECTRUM = "MS:1003067|consensus spectrum"
NUMBER_OF_REPLICATE_SPECTRA_USED = "MS:1003070|number of replicate spectra used"
NUMBER_OF_REPLICATE_SPECTRA_AVAILABLE = "MS:1003069|number of replicate spectra available"
DECOY_SPECTRUM = "MS:1003192|decoy spectrum"
DECOY_PEPTIDE_SPECTRUM = "MS:1003195|unnatural peptidoform decoy spectrum"

SCAN_POLARITY = "MS:1000465|scan polarity"
POSITIVE_SCAN = "MS:1000130|positive scan"
NEGATIVE_SCAN = "MS:1000129|negative scan"

CUSTOM_ATTRIBUTE_NAME = "MS:1003275|other attribute name"
CUSTOM_ATTRIBUTE_VALUE = "MS:1003276|other attribute value"

PEAK_OBSERVATION_FREQ = "MS:1003279|observation frequency of peak"
PEAK_ATTRIB = "MS:1003254|peak attribute"

DISSOCIATION_METHOD = "MS:1000044|dissociation method"
HCD = "MS:1000422|beam-type collision-induced dissociation"
TRAP_CID = "MS:1002472|trap-type collision-induced dissociation"
COLLISION_ENERGY = "MS:1000045|collision energy"

UNIT = "UO:0000000|unit"
ELECTRONVOLT = "UO:0000266|electronvolt"
PERCENT = "UO:0000187|percent"
SECOND = "UO:0000010|second"
MINUTE = "UO:0000031|minute"
MZ = "MS:1000040|m/z"
PPM = "UO:0000169|parts per million"


INTENSITY_OF_HIGH_UNASSIGNED_PEAK = "MS:1003289|intensity of highest unassigned peak"
TOP_20_UNASSIGNED_INTENSITY_FRACTION = "MS:1003080|top 20 peak unassigned intensity fraction"
TOTAL_UNASSIGNED_INTENSITY_FRACTION = "MS:1003079|total unassigned intensity fraction"
NUM_UNASSIGNED_PEAKS_IN_TOP_20 = "MS:1003290|number of unassigned peaks among top 20 peaks"
NUM_PEAKS = "MS:1003059|number of peaks"

Q_VALUE = "MS:1002354|PSM-level q-value"

PROTEIN_ACCESSION = "MS:1000885|protein accession"
PROTEIN_NAME = "MS:1000886|protein name"
PROTEIN_DESCRIPTION = "MS:1001088|protein description"

TAXONOMY_NCBI_TAX_ID = "MS:1001467|taxonomy: NCBI TaxID"
TAXONOMY_SCIENTIFIC_NAME = "MS:1001469|taxonomy: scientific name"
TAXONOMY_COMMON_NAME = "MS:1001468|taxonomy: common name"