import itertools
import logging
import warnings
import re
import numbers

from dataclasses import dataclass, field
from typing import Any, Callable, Deque, Dict, Iterator, List, Optional, Sequence, Tuple, Union

from psims.controlled_vocabulary.entity import Entity, ListOfType

from mzspeclib.attributes import Attribute, Attributed

from mzspeclib.spectrum import Spectrum
from mzspeclib.analyte import Analyte, Interpretation
from mzspeclib.spectrum_library import SpectrumLibrary

from mzspeclib.ontology import _VocabularyResolverMixin


from mzspeclib.validate.level import RequirementLevel
from mzspeclib.validate.semantic_rule import ScopedSemanticRule, load_rule_set
from mzspeclib.validate.object_rule import ScopedObjectRuleBase, SpectrumPeakAnnotationRule, ValidationWarning, LibraryFormatVersionFirstRule
from mzspeclib.defaults import DEFAULT_UNITS

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _walk_children(term: Entity):
    queue = Deque([term])
    while queue:
        term = queue.popleft()
        yield term
        queue.extend(term.children)


_is_curie = re.compile(r"(([A-Z]+):(.+))(?:\|.+)?")

def is_curie(value: str) -> bool:
    if isinstance(value, str):
        return _is_curie.match(value)
    return False


@dataclass
class ValidationContext:
    attributes_visited: Dict[Tuple[str, str], bool] = field(default_factory=dict)
    rule_states: Dict[str, Any] = field(default_factory=dict)

    def clear_attributes(self):
        self.attributes_visited.clear()

    def record_attribute(self, attribute: Union[Tuple[str, str], Attribute], result: bool):
        if isinstance(attribute, Attribute):
            attribute = (attribute.key, attribute.group_id)
        self.attributes_visited[attribute] = result

    def visited_attribute(self, attribute: Union[Tuple[str, str], Attribute]) -> bool:
        if isinstance(attribute, Attribute):
            attribute = (attribute.key, attribute.group_id)
        return attribute in self.attributes_visited



def _warning_iterator(iterator: Iterator[Spectrum]) -> Iterator[Spectrum]:
    # coerce to an actual iterator in case we were passed only an iterable
    iterator = iter(iterator)
    while True:
        try:
            with warnings.catch_warnings(record=True) as w:
                value = next(iterator)
            vw = [a for a in w if issubclass(a.category, ValidationWarning)]
            yield value, vw
        except StopIteration:
            break
        except:
            raise


def _is_of_type(attrib, relation) -> Tuple[bool, Optional[str]]:
    if isinstance(relation.value_type.type_definition, type):
        if relation.value_type.type_definition is float:
            tp = numbers.Number
        else:
            tp = relation.value_type.type_definition
        return (isinstance(attrib.value, tp), None)
    else:
        return _try_convert(attrib.value, relation.value_type.type_definition)


def _try_convert(value, converter) -> Tuple[bool, Optional[str]]:
    try:
        converter(value)
        return (True, None)
    except (ValueError, TypeError) as err:
        return (False, str(err))


class ValidatorBase(_VocabularyResolverMixin):
    name: str = None
    error_log: List
    current_context: ValidationContext

    def __init__(self, error_log=None, current_context=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.error_log = error_log or []
        self.current_context = current_context or ValidationContext()

    def reset_context(self):
        self.current_context.clear_attributes()

    def add_warning(
        self,
        obj: Attributed,
        path: str,
        identifier_path: Tuple,
        attrib: Any,
        value: Any,
        requirement_level: RequirementLevel,
        message: str,
    ):
        if hasattr(obj, "key"):
            key = obj.key
        elif hasattr(obj, "id"):
            key = obj.id
        else:
            key = ""
        warning = f"{attrib.id if hasattr(attrib, 'id') else attrib} failed to validate {path}:{identifier_path or key} ({requirement_level.name.upper()}): {message} ({self.name})"
        log_level = logging.WARN
        if requirement_level == RequirementLevel.may:
            log_level = logging.DEBUG
        logger.log(log_level, warning, stacklevel=2)
        self.error_log.append(ValidationError(path, identifier_path, attrib, value, requirement_level, warning, self.name))

    def validate_spectrum(
        self,
        spectrum: Spectrum,
        path: str,
        library: SpectrumLibrary,
        parsing_warnings: Optional[List[warnings.WarningMessage]] = None,
    ) -> bool:
        path = f"{path}/Spectrum"
        identifier_path = (spectrum.key,)
        result = self.apply_rules(spectrum, path, identifier_path)
        result &= self.check_attributes(spectrum, path, identifier_path)
        self.reset_context()

        if parsing_warnings:
            result = False
            for parsing_warning in parsing_warnings:
                logger.warn(str(parsing_warning.message))
                self.add_warning(
                    spectrum,
                    path,
                    identifier_path,
                    '',
                    RequirementLevel.must,
                    parsing_warning.message
                )

        for _key, analyte in spectrum.analytes.items():
            result &= self.validate_analyte(analyte, path, spectrum, library)

        for _key, interp in spectrum.interpretations.items():
            result &= self.validate_interpretation(interp, path, spectrum, library)
        return result

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary) -> bool:
        path = f"{path}/Analyte"
        identifier_path = (spectrum.key, analyte.id)
        result = self.apply_rules(analyte, path, identifier_path)
        result &= self.check_attributes(analyte, path, identifier_path)
        self.reset_context()
        return result

    def validate_interpretation(
        self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary
    ) -> bool:
        path = f"{path}/Interpretation"
        identifier_path = (spectrum.key, interpretation.id)
        result = self.apply_rules(interpretation, path, identifier_path)
        result &= self.check_attributes(interpretation, path, identifier_path)
        self.reset_context()
        return result

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        raise NotImplementedError()

    def check_attributes(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        """Stub implementation for any attribute rule checking"""
        return True

    def validate_library(self, library: SpectrumLibrary, spectrum_iterator: Optional[Iterator[Spectrum]]=None):
        path = "/Library"
        result = self.apply_rules(library, path, (library.identifier, ))
        result &= self.check_attributes(library, path, (library.identifier, ))
        self.reset_context()

        if spectrum_iterator is None:
            spectrum_iterator = library
        for spectrum, warns in _warning_iterator(spectrum_iterator):
            result &= self.validate_spectrum(spectrum, path, library, parsing_warnings=warns)
        return result

    def chain(self, validator: 'ValidatorBase') -> 'ValidatorBase':
        """
        Combine this validator with another validator, applying both rulesets.

        See Also
        --------
        ValidatorChain
        """
        return ValidatorChain([self, validator])

    def __or__(self, other: 'ValidatorBase') -> 'ValidatorBase':
        return self.chain(other)

    def walk_terms_for(self, curie: str) -> Iterator[str]:
        term = self.find_term_for(curie)
        for entity in _walk_children(term):
            yield f"{entity.id}|{entity.name}"


INTENSITY_UNITS = ("MS:1000131", "MS:1000132", "MS:1000905", " MS:1000814", "UO:0000269")

def _is_intensity_units(units: List):
    for unit_rel in units:
        if unit_rel.accession in INTENSITY_UNITS:
            return True
    return False


class ControlledVocabularyAttributeValidator(ValidatorBase):
    name: str = "ControlledVocabularyAttributeValidator"

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        return True

    def _check_value_types(self, obj: Attributed, path: str, identifier_path: Tuple, attrib, term: Entity) -> bool:
        valid: bool = True
        value_types = term.get("has_value_type")
        if not value_types:
            value_parsed = isinstance(attrib.value, Attribute)
            if value_parsed or is_curie(attrib.value):
                if value_parsed:
                    value_key = attrib.value.key.split("|")[0]
                else:
                    value_key = attrib.value.split("|")[0]
                if is_curie(value_key):
                    value_term = self.find_term_for(value_key)
                else:
                    value_term = None
                if not value_term or not value_term.is_of_type(term):
                    self.add_warning(
                        obj,
                        path,
                        identifier_path,
                        attrib.key,
                        attrib.value,
                        RequirementLevel.must,
                        f"The value type of {attrib.key} must be a term derived from {attrib.key}, but found {attrib.value}",
                    )
                    valid = False
                    return valid
            else:
                self.add_warning(
                    obj,
                    path,
                    identifier_path,
                    attrib.key,
                    attrib.value,
                    RequirementLevel.must,
                    f"The value type of {attrib.key} must be a term derived from {attrib.key}",
                )
                valid = False
                return valid
        else:
            type_message_errors = []
            for rel in value_types:
                if isinstance(rel.value_type, ListOfType):
                    hit = False
                    for tp in rel.value_type.type_definition.entity.has_value_type:
                        if isinstance(attrib.value, Sequence) and all(
                            isinstance(v, tp.value_type.type_definition) for v in attrib.value
                        ):
                            hit = True
                            break
                    if hit:
                        break
                type_match, message = _is_of_type(attrib, rel)
                if message is not None:
                    type_message_errors.append(message)
                if type_match:
                    break
            else:
                type_message_errors = ", ".join(type_message_errors)
                self.add_warning(
                    obj,
                    path,
                    identifier_path,
                    attrib.key,
                    attrib.value,
                    RequirementLevel.must,
                    f"The value type of {attrib.key} must be a value of type {', '.join([rel.value_type.id for rel in value_types])}, but got {type(attrib.value)} {type_message_errors}",
                )
                valid = False
        return valid

    def _check_units(self, obj: Attributed, path: str, identifier_path: Tuple, attrib, term: Entity) -> bool:
        valid: bool = True
        units = term.get("has_units")
        if units:
            if not isinstance(units, list):
                units = [units]

            is_intensity_measure = _is_intensity_units(units)

            if attrib.group_id is not None:
                try:
                    unit_attrib = obj.get_attribute("UO:0000000|unit", group_identifier=attrib.group_id, raw=True)
                except KeyError:
                    unit_attrib = None
                    if len(units) == 1:
                        logger.warning(f"{attrib.key}'s unit is missing, defaulting to {units[0]}")
                        return valid
            else:
                unit_attrib = None
                if len(units) == 1:
                    logger.debug(f"{attrib.key}'s unit is missing, defaulting to {units[0]}")
                    return valid

            if not unit_attrib and is_intensity_measure:
                try:
                    unit_attrib = obj.get_attribute("MS:1000043|intensity unit", raw=True)
                    if isinstance(unit_attrib, list):
                        unit_attrib = unit_attrib[0]
                except KeyError:
                    pass

            if unit_attrib:
                unit_acc, unit_name = unit_attrib.value.split("|", 1)
                for unit in units:
                    if unit_acc == unit.accession or unit_name == unit.comment:
                        break
                else:
                    self.add_warning(
                        obj,
                        path,
                        identifier_path,
                        attrib.key,
                        attrib.value,
                        RequirementLevel.must,
                        f"The attribute {attrib.key} must have a unit {', '.join([rel.accession + '|' + rel.comment for rel in units])}, but got {unit_acc}|{unit_name}",
                    )
                    valid = False
                    return valid
            else:
                if not term.id in DEFAULT_UNITS:
                    self.add_warning(
                        obj,
                        path,
                        identifier_path,
                        attrib.key,
                        attrib.value,
                        RequirementLevel.must,
                        f"The attribute {attrib.key} must have a unit {', '.join([rel.accession + '|' + rel.comment for rel in units])}, but none were found",
                    )
                    valid = False
                    return valid
        return valid

    def check_attributes(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        valid: bool = True
        for attrib in obj.attributes:
            if attrib.key == "MS:1003276|other attribute value" or attrib.key == "MS:1003275|other attribute name":
                continue
            if self.current_context.visited_attribute(attrib):
                continue
            try:
                acc, name = attrib.key.split("|", 1)
                term = self.find_term_for(acc)
                try:
                    term_by_name = self.find_term_by_name(name)
                except KeyError:
                    term_by_name = None
                if term != term_by_name and term_by_name is not None:
                    self.add_warning(
                        obj,
                        path,
                        identifier_path,
                        attrib.key,
                        attrib.value,
                        RequirementLevel.must,
                        f'{attrib.key}\'s accession {acc} resolved to "{term.id}|{term.name}", '
                        f'but it\'s name "{name}" resolved to "{term_by_name.id}|{term_by_name.name}"',
                    )
                    valid = False
                elif term_by_name is None:
                    self.add_warning(
                        obj,
                        path,
                        identifier_path,
                        attrib.key,
                        attrib.value,
                        RequirementLevel.must,
                        f'{attrib.key}\'s accession {acc} resolved to "{term.id}|{term.name}", '
                        f'but it\'s name "{name}" did not resolve',
                    )
                    valid = False
            except KeyError:
                logger.warn(f"Could not resolve term for {attrib.key} at {path} {identifier_path}")
                valid = False
                self.add_warning(
                    obj,
                    path,
                    identifier_path,
                    attrib,
                    None,
                    RequirementLevel.must,
                    f"Could not resolve term for {attrib.key} at {path} {identifier_path}",
                )
                continue
            except ValueError:
                logger.error(f"Could not parse {attrib.key} at {path} {identifier_path}", exc_info=True)
                valid = False
                self.add_warning(
                    obj,
                    path,
                    identifier_path,
                    attrib,
                    None,
                    RequirementLevel.must,
                    f"Could not parse {attrib.key} at {path} {identifier_path}",
                )

            valid &= self._check_value_types(obj, path, identifier_path, attrib, term)
            valid &= self._check_units(obj, path, identifier_path, attrib, term)

        return valid


@dataclass(frozen=True)
class ValidationError:
    path: str
    identifier_path: Tuple
    attribute: Any
    value: Any
    requirement_level: RequirementLevel
    message: str
    source: str = None

    def __hash__(self):
        return hash((self.path, self.identifier_path, self.attribute, self.message))


class RuleValidator(ValidatorBase):
    name: str
    semantic_rules: List[ScopedSemanticRule]
    object_rules: List[ScopedObjectRuleBase]

    def __init__(self, name, semantic_rules: Optional[List[ScopedSemanticRule]] = None, object_rules: Optional[List[ScopedObjectRuleBase]] = None, error_log: Optional[List] = None):
        super().__init__(error_log)
        self.name = name
        self.semantic_rules = semantic_rules or []
        self.object_rules = object_rules or []

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        result = True
        for rule in itertools.chain(self.semantic_rules, self.object_rules):
            if rule.path == path:
                v = rule(obj, path, identifier_path, self)
                level = logging.DEBUG
                if not v and rule.requirement_level > RequirementLevel.may:
                    # level = logging.WARN
                    result &= v
                logger.log(level, f"Applied {rule.id} to {path}:{identifier_path} {v}/{result}")
        return result

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name!r}, {self.semantic_rules}, {self.object_rules})"


class ValidatorChain(ValidatorBase):
    validators: List[ValidatorBase]

    def __init__(self, validators: List[ValidatorBase], *args, **kwargs):
        self.validators = list(validators)
        super().__init__(*args, **kwargs)

    @property
    def error_log(self):
        log = []
        seen = set()
        for validator in self.validators:
            for message in validator.error_log:
                if message not in seen:
                    log.append(message)
                    seen.add(message)
        return log

    @error_log.setter
    def error_log(self, _value):
        pass

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary, parsing_warnings: Optional[List[warnings.WarningMessage]] = None):
        result = True
        for validator in self.validators:
            result &= validator.validate_spectrum(spectrum, path, library, parsing_warnings)
        return result

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        result = True
        for validator in self.validators:
            result &= validator.validate_analyte(analyte, path, spectrum, library)
        return result

    def validate_interpretation(self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        result = True
        for validator in self.validators:
            result &= validator.validate_interpretation(interpretation, path, spectrum, library)
        return result

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        result = True
        for validator in self.validators:
            result &= validator.apply_rules(obj, path, identifier_path)
        return result

    def check_attributes(self, obj: Attributed, path: str, identifer_path: Tuple) -> bool:
        result = True
        for validator in self.validators:
            result &= validator.check_attributes(obj, path, identifer_path)
        return result

    def reset_context(self):
        for validator in self.validators:
            validator.reset_context()

    def chain(self, validator: ValidatorBase) -> ValidatorBase:
        self.validators.append(validator)
        return self


predicates = {
    "single_spectrum": lambda spec: spec.get_attribute("MS:1003065|spectrum aggregation type") == "MS:1003066|singleton spectrum",
    "consensus_spectrum": lambda spec: spec.get_attribute("MS:1003065|spectrum aggregation type") == "MS:1003066|singleton spectrum",
}


object_rules = {
    "peak_annotations": [SpectrumPeakAnnotationRule()],
    "base": [LibraryFormatVersionFirstRule()]
}


def get_validator_for(name: str) -> RuleValidator:
    """Load a :class:`RuleValidator` with semantic rules for the given rule profile"""
    rules = load_rule_set(name)
    validator = RuleValidator(name, rules)
    return validator


def get_object_validator_for(name: str) -> RuleValidator:
    """Prepare a :class:`RuleValidator` configured with the object rules for the given rule profile"""
    rules = object_rules[name]
    validator = RuleValidator(name, object_rules=rules)
    return validator


def load_default_validator() -> ValidatorChain:
    """Load the core :class:`RuleValidator` and add all the basic object rules for library validation"""
    chain = get_validator_for("base")
    chain |= ControlledVocabularyAttributeValidator()
    chain |= get_object_validator_for("base")
    chain |= get_object_validator_for("peak_annotations")
    return chain
