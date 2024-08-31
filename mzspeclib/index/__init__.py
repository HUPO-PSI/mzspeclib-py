"""
Spectral Library Indexing
=========================

:title-reference:`MzSpecLib` recommends that each file format be indexed in
a way that is appropriate to the format and the application, and
recognizes that most formats will not have any kind of built-in
index.

:mod:`mzspeclib` provides data structures for both in-memory and on-disk
index storage.

The :class:`~.IndexBase` type is an abstract base class for writing
indices. It isn't necessary unless you're writing a new backend or
index format.

The :class:`~.MemoryIndex` type holds all offset information and metadata
in memory, which can make it fast to access but problematic if several
large libraries are open at once. The index is not saved, requiring a
potentially costly scan of the entire library file when opening a file.

The :class:`~.SQLIndex` type holds its information in a SQLite3 database on
disk, and executes queries to bring only the required information into memory
when requested. There is some I/O overhead involved in each look-up, but
the index information persists between runs which can greatly improve start-up
time.

"""

from .memory import MemoryIndex
from .sql import SQLIndex
from .base import IndexBase

__all__ = [
    'IndexBase',
    'MemoryIndex',
    'SQLIndex'
]