import logging
import sys

from typing import List, DefaultDict

from mzspeclib import SpectrumLibrary
from mzspeclib.validate import *

logging.basicConfig(level=logging.DEBUG, stream=sys.stderr)
logger = logging.getLogger()

library = SpectrumLibrary(filename="examples/fetal_brain_tiny.mzlib.txt")

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