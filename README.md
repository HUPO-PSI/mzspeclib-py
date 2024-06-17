# PSI Spectral Library format - Python implementation

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/HUPO-PSI/mzspeclib-py/test.yml?style=for-the-badge)
![Read the Docs](https://img.shields.io/readthedocs/mzspeclib-py?style=for-the-badge)

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

Test with
```
pytest
````
