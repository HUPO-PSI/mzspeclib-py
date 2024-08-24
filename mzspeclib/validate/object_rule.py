
import logging

from typing import TYPE_CHECKING, List, Tuple

from mzspeclib.attributes import Attributed

from mzspeclib.spectrum import Spectrum
from mzspeclib.annotation import IonAnnotationBase, InvalidAnnotation

from mzspeclib.validate.level import RequirementLevel

if TYPE_CHECKING:
    from .validator import ValidatorBase
    from mzspeclib.spectrum_library import SpectrumLibrary

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class ValidationWarning(UserWarning):
    """
    Indicates that something was parsed that did not halt the parser but
    which violates the expectations of the parser.

    The parser will make a best-effort attempt to interpret the value
    correctly but when validating this will count as a violation.
    """


class ScopedObjectRuleBase:
    id: str
    path: str
    requirement_level: RequirementLevel

    def __init__(self, id, path: str, requirement_level: RequirementLevel=RequirementLevel.should) -> None:
        self.id = id
        self.path = path
        self.requirement_level = requirement_level

    def __eq__(self, other):
        if other is None:
            return False
        props_eq = self.id == other.id and self.path == other.path and self.requirement_level == other.requirement_level
        if not props_eq:
            return False
        return isinstance(other, self.__class__) or isinstance(self, other.__class__)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.id)

    def __call__(self, obj: Attributed, path: str, identifier_path: Tuple, validator_context: "ValidatorBase") -> bool:
        return self.validate(obj, path, identifier_path, validator_context)

    def validate(self, obj: Attributed, path: str, identifier_path: Tuple, validator_context: "ValidatorBase") -> bool:
        raise NotImplementedError()


class LibraryFormatVersionFirstRule(ScopedObjectRuleBase):
    def __init__(self, requirement_level: RequirementLevel = RequirementLevel.must):
        super().__init__("Library_format_version_first_rule", "/Library", requirement_level=requirement_level)

    def validate(self, obj: "SpectrumLibrary", path: str, identifier_path: Tuple, validator_context: "ValidatorBase") -> bool:
        errors = []
        attr = next(iter(obj.attributes))
        is_format_version_key = attr.key == "MS:1003186|library format version"
        if not is_format_version_key:
            validator_context.add_warning(
                obj,
                path,
                identifier_path,
                self,
                attr,
                self.requirement_level,
                f"The first attribute of the library is not 'MS:1003186|library format version': {attr}",
            )
        if errors:
            return False
        return True


class SpectrumPeakAnnotationRule(ScopedObjectRuleBase):
    def __init__(self, requirement_level: RequirementLevel=RequirementLevel.should):
        super().__init__(
            "Spectrum_peak_annotation_rule", "/Library/Spectrum", requirement_level=requirement_level)

    def validate(self, obj: Spectrum, path: str, identifier_path: Tuple, validator_context: "ValidatorBase") -> bool:
        errors = []
        for i, peak in enumerate(obj.peak_list):
            annotations = peak[2]
            annotations: List[IonAnnotationBase]
            for annot in annotations:
                if isinstance(annot, InvalidAnnotation):
                    errors.append(annot)
                    validator_context.add_warning(
                        obj,
                        path,
                        identifier_path,
                        self,
                        annot,
                        self.requirement_level,
                        f"Invalid peak annotation for peak {i} of Spectrum={obj.key} {annot.content}"
                    )
                else:
                    analyte_id = annot.analyte_reference
                    if analyte_id is not None:
                        if analyte_id not in obj.analytes:
                            validator_context.add_warning(
                                obj,
                                path,
                                identifier_path,
                                self,
                                annot,
                                self.requirement_level,
                                (f"Analyte of peak annotation {annot} for peak {i} of "
                                 f"Spectrum={obj.key} {annot.analyte_reference} is missing")
                            )


        if errors:
            return False
        return True
