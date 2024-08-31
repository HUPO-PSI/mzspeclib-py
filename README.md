# PSI Spectral Library format - Python implementation

[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/HUPO-PSI/mzspeclib-py/test.yml?style=for-the-badge)](https://github.com/HUPO-PSI/mzspeclib-py/actions/workflows/test.yml)
[![Read the Docs](https://img.shields.io/readthedocs/mzspeclib-py?style=for-the-badge)](https://mzspeclib-py.readthedocs.io/en/latest/)

This is a Python reference implementation of the [`mzSpecLib`](https://www.psidev.info/mzspeclib)
spectral library format. It provides readers and writers for the Text and JSON serialization
of `mzSpecLib`, as well as readers for the following spectral library formats found in the wild:

- NIST MSP (also minimal writing facility)
- SPTXT
- BiblioSpec
- EncyclopeDIA
- DIA-NN
- Spectronaut

Once installed, it can be used programmatically by importing the `mzspeclib` Python library, or using
the `mzspeclib` command line tool to read, write, and manipulate spectral libraries.

## Development

For development, run:
```sh
pip install --editable .
```

Test with
```
pytest
````

## Usage

### Python library

`mzspeclib-py` provides a Python API with the name `mzspeclib`. The top-level type, `SpectrumLibrary` is
the main entry point into the library:
<!-- [[[cog
    import cog
    cog.outl("```python")
    cog.outl(open("docs/examples/readme.py").read())
    cog.outl("```")
]]] -->
```python
from mzspeclib import SpectrumLibrary

# Open a spectrum library
lib = SpectrumLibrary(filename="examples/fetal_brain_tiny.mzlib.txt")
print(lib)

# Get the number of spectra in the library
n_spectra = len(lib)
print(n_spectra)

# Get a specific spectrum from the library
spec = lib.get_spectrum(spectrum_number=3)
print(f"Key={spec.key}; Name={spec.name}; Num Peaks={len(spec.peak_list)}")
print(spec.get_interpretation('1'))
```
<!-- [[[end]]] -->

### Command Line Tool

It also provides a command line tool, also called `mzspeclib` that lets the user convert
supported spectrum library formats into one of the PSI MzSpecLib formats.

All of the commands provide a limited automatic file format detection, but you can specify the
input format if needed.

<!-- [[[cog
    import cog
    import subprocess

    buf = subprocess.check_output(["mzspeclib", "--help"])
    cog.outl("```bash")
    cog.outl("$ mzspeclib --help")
    cog.outl(buf.decode('utf8'))
    cog.outl("```")
]]] -->
```bash
$ mzspeclib --help
Usage: mzspeclib [OPTIONS] COMMAND [ARGS]...

  A collection of utilities for inspecting and manipulating spectral
  libraries.

Options:
  -d, --debug-logging  Enable debug logging
  -l, --log-file PATH  Write log messages to this file as well as STDERR
  -h, --help           Show this message and exit.

Commands:
  convert   Convert a spectral library from one format to another
  describe  Produce a minimal textual description of a spectral library
  index     Build an on-disk index for a spectral library
  validate  Semantically validate a spectral library

```
<!-- [[[end]]] -->

#### File Conversion

`mzspeclib convert` can read a variety of different text and binary file formats and write
them out in PSI MzSpecLib's Text and JSON encodings, as well as a dialect of MSP.

<!-- [[[cog
    import cog
    import subprocess

    buf = subprocess.check_output(["mzspeclib", "convert", "--help"])
    cog.outl("```bash")
    cog.outl("$ mzspeclib convert --help")
    cog.outl(buf.decode('utf8'))
    cog.outl("```")
]]] -->
```bash
$ mzspeclib convert --help
Usage: mzspeclib convert [OPTIONS] INPATH OUTPATH

  Convert a spectral library from one format to another. If `outpath` is `-`,
  instead of writing to file, data will instead be sent to STDOUT.

Options:
  -i, --input-format [bibliospec|blib|dia-nn.tsv|dlib|encyclopedia|json|msp|mzlb.json|mzlb.txt|mzlib.json|mzlib.txt|spectronaut.tsv|sptxt|text]
                                  The file format of the input file. If
                                  omitted, will attempt to infer
                                  automatically.
  -f, --format [text|json|msp]
  -k, --library-attribute <TEXT TEXT>...
                                  Specify an attribute to add to the library
                                  metadata section. May be repeated.
  -K, --header-file PATH          Specify a file to read name-value pairs
                                  from. May be JSON or TAB-separated
  -h, --help                      Show this message and exit.

```
<!-- [[[end]]] -->

#### Semantic Validation

`mzspeclib` includes a reference implementation of a semantic validator tool for spectrum libraries, testing
whether or not a library is compliant with the specification and a subset of community recommended practices.

<!-- [[[cog
    import cog
    import subprocess

    buf = subprocess.check_output(["mzspeclib", "validate", "--help"])
    cog.outl("```bash")
    cog.outl("$ mzspeclib validate --help")
    cog.outl(buf.decode('utf8'))
    cog.outl("```")
]]] -->
```bash
$ mzspeclib validate --help
Usage: mzspeclib validate [OPTIONS] INPATH

  Semantically and structurally validate a spectral library.

Options:
  -p, --profile [consensus|single|silver|peptide|gold]
  -i, --input-format [bibliospec|blib|dia-nn.tsv|dlib|encyclopedia|json|msp|mzlb.json|mzlb.txt|mzlib.json|mzlib.txt|spectronaut.tsv|sptxt|text]
                                  The file format of the input file. If
                                  omitted, will attempt to infer
                                  automatically.
  -h, --help                      Show this message and exit.

```
<!-- [[[end]]] -->