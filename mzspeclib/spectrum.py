from __future__ import print_function

import textwrap

from typing import Any, Dict,  List, Optional, TYPE_CHECKING

from mzspeclib.attributes import (
    AttributeManager, AttributeManagedProperty, AttributeListManagedProperty,
    AttributeProxy as _AttributeProxy, AttributeFacet
)
from mzspeclib.analyte import Analyte, InterpretationCollection, Interpretation
from .const import (SPECTRUM_NAME, LIBRARY_SPECTRUM_INDEX, LIBRARY_SPECTRUM_KEY, PRECURSOR_MZ, CHARGE_STATE,
                    SPECTRUM_AGGREGATION_TYPE, NUMBER_OF_REPLICATE_SPECTRA_USED, NUMBER_OF_REPLICATE_SPECTRA_AVAILABLE,
                    THEORETICAL_MZ, PEAK_ATTRIBUTE)

if TYPE_CHECKING:
    from mzspeclib.spectrum_library import SpectrumLibrary

#A class that holds data for each spectrum that is read from the SpectralLibrary class


class SpectrumAggregation(_AttributeProxy):
    aggregation_type = AttributeManagedProperty(SPECTRUM_AGGREGATION_TYPE)
    replicates_used = AttributeManagedProperty[int](NUMBER_OF_REPLICATE_SPECTRA_USED)
    replicates_available = AttributeManagedProperty[int](NUMBER_OF_REPLICATE_SPECTRA_AVAILABLE)


class Spectrum(AttributeManager):
    """
    A mass spectrum stored in a spectral library.

    The :class:`Spectrum` type is a :class:`~.AttributeManager`, interacting
    with the spectrum's attribute list itself. It also encloses :class:`~.Interpretation`
    and :class:`~.Analyte` members which have their own attributes.

    Attributes
    ----------
    peak_list : List[Tuple[float, float, Optional[List[:class:`mzpaf.PeakAnnotation`]], Optional[List[Any]]]]
        The list of peaks which may or may not be annotated and may or may not have
        aggregation metrics.
    interpretations : :class:`~.InterpretationCollection`
        The set of all interpretations associated with this spectrum, and their analytes.
    """

    peak_list: List
    analytes: Dict[str, Analyte]
    interpretations: InterpretationCollection
    _source: Optional['SpectrumLibrary']

    #### Constructor
    def __init__(self, attributes=None, peak_list=None, analytes=None,
                 interpretations=None):
        """

        Parameters
        ----------
        attributes : list
            A list of attribute [key, value (, group)] sets to initialize to.
        peak_list : list
            A list of tuples representing (annotated) peaks
        analytes : dict[str, :class:`~.Analyte`]
            A mapping from identifier to :class:`~.Analyte` unique within this
            :class:`Spectrum`.
        interpretations : :class:`~.InterpretationCollection`
            A mapping from identifier to :class:`~.Interpretation` unique within
            this :class:`Spectrum`.
        """
        if peak_list is None:
            peak_list = []
        if analytes is None:
            analytes = dict()
        if interpretations is None:
            interpretations = InterpretationCollection()
        else:
            interpretations = InterpretationCollection(interpretations)
        super(Spectrum, self).__init__(attributes)
        self.peak_list = peak_list
        self.analytes = analytes
        self.interpretations = interpretations

    name = AttributeManagedProperty[str](SPECTRUM_NAME)
    key = AttributeManagedProperty[int](LIBRARY_SPECTRUM_KEY)
    index = AttributeManagedProperty[int](LIBRARY_SPECTRUM_INDEX)

    precursor_mz = AttributeListManagedProperty[float](
        [PRECURSOR_MZ,
         THEORETICAL_MZ])

    _precursor_charge = AttributeManagedProperty[int](CHARGE_STATE)

    spectrum_aggregation = AttributeFacet[SpectrumAggregation](SpectrumAggregation)
    peak_aggregations = AttributeManagedProperty(PEAK_ATTRIBUTE, multiple=True)

    @property
    def precursor_charge(self) -> int:
        """Obtain the spectrum's precursor ion charge or analyte charge"""
        if self._precursor_charge:
            return self._precursor_charge
        for analyte in self.analytes.values():
            if analyte.charge:
                return analyte.charge

    def add_analyte(self, analyte: Analyte):
        """Add an :class:`~.Analyte` to the spectrum"""
        self.analytes[str(analyte.id)] = analyte

    def get_analyte(self, analyte_id: str) -> Analyte:
        """Get an :class:`~.Analyte` by ID from the spectrum"""
        return self.analytes[str(analyte_id)]

    def remove_analyte(self, analyte_id: str):
        """Remove an :class:`~.Analyte` by ID from the spectrum"""
        analyte_id = str(analyte_id)
        del self.analytes[analyte_id]
        interpretation: Interpretation
        for interpretation in self.interpretations.values():
            if interpretation.has_analyte(analyte_id):
                interpretation.remove_analyte(analyte_id)

    def add_interpretation(self, interpretation: Interpretation):
        self.interpretations.add_interpretation(interpretation)

    def get_interpretation(self, interpretation_id: str) -> Interpretation:
        return self.interpretations.get_interpretation(interpretation_id)

    def __eq__(self, other):
        result = super(Spectrum, self).__eq__(other)
        if result:
            result = self.peak_list == other.peak_list
        if result:
            result = self.analytes == other.analytes
        return result

    def __repr__(self):  # pragma: no cover
        template = f"{self.__class__.__name__}("
        lines = list(map(str, self.attributes))
        if not lines:
            template += "[], "
        else:
            template += "[\n%s], " % textwrap.indent(',\n'.join(lines), ' ' * 2)
        lines = list(map(str, self.peak_list))
        if not lines:
            template += "peak_list=[])"
        else:
            template += "peak_list=[\n%s])" % textwrap.indent(
                ',\n'.join(lines), ' ' * 2)
        return template

    def __str__(self):  # pragma: no cover
        return self.write("text")

    def write(self, format="text", **kwargs):  # pragma: no cover
        """
        Write out the spectrum in any of the supported formats

        Parameters
        ----------
        format : str
            The name of the format to write in
        **kwargs
            Passed to implementation
        """
        #### Set a buffer to fill with string data
        buffer = ''

        #### Make the format string lower case to facilitate comparisons
        format = format.lower()

        #### If the format is text
        if format == "text":
            from mzspeclib.backends.text import format_spectrum
            return format_spectrum(self, **kwargs)

        #### If the format is TSV
        elif format == "tsv" or format == "csv":

            #### Set the appropriate delimiter for format
            delimiter = "\t"
            if format == "csv":
                delimiter = ","

            #### Create the header line and columns line
            buffer += "# Spectrum attributes\n"
            buffer += delimiter.join("cv_param_group", "accession", "name", "value_accession", "value\n")

            for attribute in self.attributes:
                if len(attribute) == 2:
                    key,value = attribute
                    cv_param_group = ''
                elif len(attribute) == 3:
                    key,value,cv_param_group = attribute
                    if format == "csv" and ',' in str(value):
                        value = '"' + value + '"'
                else:
                    print("ERROR: Unsupported number of items in attribute")
                    print(attribute)
                    raise ValueError(
                        f"Unsupported number of items in attribute: {attribute}")
                components = key.split('|',1)
                if len(components) == 2:
                    accession,name = components
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    raise ValueError(f"Unsupported number of items in components: {components}")
                components = str(value).split('|',1)
                if len(components) == 2:
                    value_accession,value = components
                    value = str(value)
                    if format == "csv" and ',' in value:
                        value = '"' + value + '"'
                elif len(components) == 1:
                    value = str(value)
                    if format == "csv" and ',' in value:
                        value = '"' + value + '"'
                    value_accession = ''
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    raise ValueError(
                        f"Unsupported number of items in components: {components}")

                #### Create the data line
                buffer += delimiter.join(map(str, [cv_param_group,accession,name,value_accession,value]))+"\n"

            #### Create the header line and columns line
            buffer += "# Peak list\n"
            buffer += "mz\tintensity\tinterpretation\n"

            #### Write out the peak list
            for peak in self.peak_list:
                mz,intensity,interpretation = peak
                if format == "csv" and ',' in str(interpretation):
                    interpretation = '"' + interpretation + '"'
                buffer += delimiter.join(map(str, [mz,intensity,interpretation]))+"\n"

        #### If the format is JSON
        elif format == "json":
            from mzspeclib.backends.json import format_spectrum
            return format_spectrum(self, **kwargs)

        #### Otherwise we don't know this format
        else:
            raise ValueError(f"ERROR: Unrecogized format '{format}'")

        return(buffer)
