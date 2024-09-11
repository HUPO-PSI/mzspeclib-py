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

# Loop over a library and count the number of unique peptide analytes
unique_analytes = set()
n_spectra = 0
for spec in lib:
    for analyte in spec.analytes.values():
        unique_analytes.add(str(analyte.peptide))
    n_spectra += 1

print(f"\n{len(unique_analytes)} unique analytes over {n_spectra} spectra")
print(unique_analytes)