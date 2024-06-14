# PSI Spectral Library format - Python implementation

This is a Python reference implementation of the [`mzSpecLib`](https://www.psidev.info/mzspeclib)
spectral library format. It provides readers and writers for the Text and JSON serialization
of `mzSpecLib`, as well as readers for the following spectral library formats found in the wild:

- NIST MSP (also minimal writing facility)
- BiblioSpec
- EncyclopeDIA
- DIA-NN
- Spectronaut

Once installed, it can be used programmatically by importing the `mzlib` Python library, or using
the `mzspeclib` command line tool to read, write, and manipulate spectral libraries.

## Development

For development, run:
```sh
pip install --editable .
```

Or from the root of this repository:
```sh
pip install --editable implementations/python
```

Test with
```
pytest
````
