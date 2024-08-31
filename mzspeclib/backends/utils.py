"""Helper module for I/O operations"""

import os
import io
import gzip
import logging

from collections import deque
from urllib import parse as urlparse
from typing import Any, Dict, Iterable, Mapping, Optional, Union

DEFAULT_BUFFER_SIZE = int(2e6)
GZIP_MAGIC = b"\037\213"

GzipFile = gzip.GzipFile
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

try:
    # Fast random acces with Gzip compatibility
    import idzip

    idzip.compressor.IdzipWriter.enforce_extension = False

    GzipFile = idzip.IdzipFile
except ImportError:
    pass


class _LineBuffer(object):
    """
    An implementation detail that treats a stream/iterator over line strings as LIFO
    queue that can have lines pushed back onto it.
    """

    lines: deque
    stream: io.IOBase
    last_line: str
    encoding: Optional[str]
    _stream_is_file_like: bool

    def __init__(
        self,
        stream: io.IOBase,
        lines: Iterable = None,
        last_line: str = None,
        encoding: Optional[str] = None,
    ):
        if lines is None:
            lines = []
        self.lines = deque(lines)
        self.stream = stream
        self.last_line = last_line
        self.encoding = encoding
        self._stream_is_file_like = hasattr(self.stream, "readline")

    def readline(self) -> Union[bytes, str]:
        if self.lines:
            line = self.lines.popleft()
        else:
            line = (
                self.stream.readline()
                if self._stream_is_file_like
                else next(self.stream)
            )
        self.last_line = line
        if self.encoding:
            return line.decode(self.encoding)
        return line

    def push_line(self, line=None):
        if line is None:
            line = self.last_line
            self.last_line = None
        if line is None:
            raise ValueError(
                "Cannot push empty value after the backtrack line is consumed"
            )
        self.lines.appendleft(line)

    def __iter__(self):
        while self.lines:
            line = self.lines.popleft()
            self.last_line = line
            if self.encoding:
                yield line.decode(self.encoding)
            else:
                yield line
        for line in self.stream:
            self.last_line = line
            if self.encoding:
                yield line.decode(self.encoding)
            else:
                yield line

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


def try_cast(value: Any) -> Union[str, int, float, Any]:
    """
    Given a value, if it is a string, attempt to convert it to a numeric type,
    or else return it as is.
    """
    if value is None:
        return value
    if not isinstance(value, str):
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def test_gzipped(f) -> bool:
    """
    Checks the first two bytes of the
    passed file for gzip magic numbers

    Parameters
    ----------
    f : file-like or path-like
        The file to test

    Returns
    -------
    bool
    """
    if isinstance(f, os.PathLike):
        f = io.open(f, "rb")
    try:
        current = f.tell()
        assert current >= 0
    except OSError:
        return False
    f.seek(0)
    magic = f.read(2)
    f.seek(current)
    return magic == GZIP_MAGIC


def starts_with_gz_magic(bytestring):
    """
    Test whether or not a byte string starts with
    the GZIP magic bytes.

    Parameters
    ----------
    bytestring : bytes
        The bytes to test.

    Returns
    -------
    bool
    """
    return bytestring.startswith(GZIP_MAGIC)


class _NotClosingWrapper:
    stream: io.BufferedIOBase

    def __init__(self, stream) -> None:
        self.stream = stream

    def __getattr__(self, attrib: str):
        attr = getattr(self.stream, attrib)
        return attr

    def close(self):
        logger.debug("Resetting stream handle %r", self.stream)
        if self.stream.seekable():
            self.stream.seek(0)

    def __iter__(self):
        return iter(self.stream)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        logger.debug("Resetting stream handle %r", self.stream)
        if self.stream.seekable():
            self.stream.seek(0)
        return


def open_stream(
    f: Union[io.IOBase, os.PathLike],
    mode="rt",
    buffer_size: Optional[int] = None,
    encoding: Optional[str] = "utf8",
    newline=None,
    closing=False,
):
    """
    Select the file reading type for the given path or stream.

    Detects whether the file is gzip encoded.
    """
    if isinstance(f, PreBufferedStreamReader):
        if closing:
            return f
        else:
            return _NotClosingWrapper(f)
    if buffer_size is None:
        buffer_size = DEFAULT_BUFFER_SIZE
    is_stream = False
    offset = None
    if "r" in mode:
        if not hasattr(f, "read"):
            f = io.open(f, "rb")
        else:
            is_stream = True
            offset = f.tell()
        if not isinstance(f, io.BufferedReader) and not isinstance(f, io.TextIOWrapper):
            if not closing and is_stream:
                f = _NotClosingWrapper(f)
            buffered_reader = io.BufferedReader(f, buffer_size)
        else:
            if not closing and is_stream:
                f = _NotClosingWrapper(f)
            buffered_reader = f
        if test_gzipped(buffered_reader):
            handle = GzipFile(fileobj=buffered_reader, mode="rb")
        else:
            handle = buffered_reader
    else:
        raise NotImplementedError(
            "Haven't implemented automatic output stream determination"
        )
    try:
        fmode = f.mode
        if isinstance(fmode, int):
            # gzip.GzipFile uses ints in `mode`
            fmode = "b"
    except AttributeError:
        fmode = "b"
    if "b" not in mode and "b" in fmode:
        handle = io.TextIOWrapper(handle, encoding=encoding, newline=newline)
    if is_stream and handle.seekable():
        if offset:
            handle.seek(offset)
    return handle


class CaseInsensitiveDict(Dict[str, Any]):
    """A case sensitive version of a dictionary with string keys."""

    def __init__(self, base=None, **kwargs):
        if base is not None:
            self.update(base)
        if kwargs:
            self.update(kwargs)

    def __getitem__(self, key: str):
        return super().__getitem__(key.lower())

    def __setitem__(self, key: str, value):
        super().__setitem__(key.lower(), value)

    def __delitem__(self, key: str):
        super().__delitem__(key.lower())

    def __contains__(self, __o: str) -> bool:
        if hasattr(__o, "lower"):
            return super().__contains__(__o.lower())
        else:
            return super().__contains__(__o)

    def get(self, key: str, default=None):
        return super().get(key.lower(), default)

    def update(self, value: Mapping[str, Any]):
        super().update({k.lower(): v for k, v in value.items()})


def urlify(path: str) -> str:
    """Convert a path into a URL if it is not already one."""
    parsed = urlparse.urlparse(path)
    if parsed.scheme == "":
        parsed = parsed._replace(scheme="file")
    return urlparse.urlunparse(parsed)


import io
from typing import Literal

DEFAULT_BUFFER_SIZE = 2**18


class PreBufferedStreamReader(io.IOBase):
    """
    A file-like object that can wrap a non-rewindable (e.g. un-seekable)
    file stream like :obj:`sys.stdin` and support seeking over the first **N**
    bytes, but blocking seek operations once the stream has been read beyond
    that point.

    This type exists to adapt :obj:`sys.stdin` or socket files to work with the
    behavior of readers like :class:`~ms_deisotope.data_source.mzml.MzMLLoader`
    that rely on early random access to extract metadata even when full random
    access indexing isn't enabled.

    Attributes
    ----------
    stream : :class:`io.IOBase`
        The object to actually read from. Must support :meth:`io.IOBase.read`,
        :meth:`io.IOBase.tell`, and :meth:`close`. It will also proxy
        :meth:`io.IOBase.fileno` if it is available.
    """

    stream: io.IOBase
    buffer: io.BytesIO

    _cursor: int
    _total_cursor: int
    _buffer_size: int

    def __init__(self, stream, buffer_size=DEFAULT_BUFFER_SIZE):
        self.stream = stream
        dat = self.stream.read(buffer_size)
        self.buffer = io.BytesIO(dat)
        self._cursor = 0
        self._total_cursor = len(dat)
        self._buffer_size = buffer_size

    def read(self, n=-1):
        if n == -1:
            buffer = bytearray()
            buffer.extend(self.buffer.read())
            buffer.extend(self.stream.read())
            self._cursor = self._buffer_size
            self._total_cursor = len(buffer) - self._cursor
            return bytes(buffer)
        buffer = bytearray(n)
        view = memoryview(buffer)
        k = self.buffer.readinto(view)
        self._cursor += k
        j = 0
        if k < n:
            j = self.stream.readinto(view[k:])
            self._total_cursor += j
            k += j
        return bytes(view[:k])

    def fileno(self):
        return self.stream.fileno()

    def seek(self, offset: int, whence: Literal[0, 1, 2] = io.SEEK_SET):
        if whence == io.SEEK_SET and offset > self._buffer_size:
            raise io.UnsupportedOperation("Cannot seek beyond the pre-loaded buffer")
        elif whence == io.SEEK_CUR or whence == io.SEEK_END:
            raise io.UnsupportedOperation("Cannot seek not from the start of the file")
        if self._total_cursor > self._buffer_size:
            raise io.UnsupportedOperation("Already read beyond the pre-loaded buffer")
        else:
            self.buffer.seek(offset, whence)
            self._cursor = offset

    def truncate(self, size=None):
        raise io.UnsupportedOperation("Read-only")

    def close(self):
        self.buffer.close()
        return self.stream.close()

    @property
    def closed(self):
        return self.stream.closed

    def seekable(self):
        return True

    def tell(self):
        if self._cursor < self._buffer_size:
            return self._cursor
        return self._total_cursor

    def writable(self):
        return False

    def __enter__(self):
        self.stream.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
