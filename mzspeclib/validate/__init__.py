"""
The components of the semantic validator for :title-reference:`MzSpecLib`

To obtain the base validator, use :func:`~.load_default_validator`, and add
successive validation rule sets using :func:`~.get_validator_for` and add
them to the base rules using the :meth:`~.ValidatorBase.chain` or the `|=`
operator.

.. code-block:: python

    from mzspeclib.validate import *

    chain = load_default_validator()
    chain.validate_library(library)

    by_level: DefaultDict[RequirementLevel, List[ValidationError]] = DefaultDict(list)
    for message in chain.error_log:
        by_level[message.requirement_level].append(message)

    for level, bucket in sorted(by_level.items()):
        log_level = logging.WARN
        if level == RequirementLevel.may:
            log_level = logging.DEBUG
        logger.log(log_level, f"Found {len(bucket)} violations for {level.name.upper()} rules")
        for err in bucket:
            logger.log(log_level, f"... {err.message}")
"""

from .level import CombinationLogic, RequirementLevel
from .semantic_rule import (RuleSet, ScopedSemanticRule, AttributeSemanticPredicate, AttributeSemanticRule)
from .object_rule import (LibraryFormatVersionFirstRule, SpectrumPeakAnnotationRule, ScopedObjectRuleBase)
from .validator import (
    ValidationError, ValidationWarning, ControlledVocabularyAttributeValidator, ValidatorBase,
    RuleValidator, ValidatorChain, load_default_validator, get_object_validator_for, get_validator_for
)