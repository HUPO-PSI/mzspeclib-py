from typing import Iterable, Any, Union, TypeVar

T = TypeVar("T")

def ensure_iter(obj: Union[Any, T, Iterable[T]]) -> Iterable[Union[Any, T]]:
    if isinstance(obj, Iterable):
        if not isinstance(obj, (str, bytes)):
            return obj
    return (obj, )


def flatten(obj: Iterable[Union[Iterable[T], T]]) -> Iterable[T]:
    acc = []
    for val in obj:
        if isinstance(val, Iterable) and not isinstance(val, (str, bytes)):
            acc.extend(flatten(val))
        else:
            acc.append(val)
    return acc


class ValidationWarning(UserWarning):
    """
    Indicates that something was parsed that did not halt the parser but
    which violates the expectations of the parser.

    The parser will make a best-effort attempt to interpret the value
    correctly but when validating this will count as a violation.
    """