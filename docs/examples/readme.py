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