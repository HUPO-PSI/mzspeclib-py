"""A little Streamlit viewer for mzSpecLib"""

import pandas as pd
import matplotlib

matplotlib.use("agg")
from matplotlib import pyplot as plt

import mzspeclib
from mzspeclib import SpectrumLibrary, Spectrum

import streamlit as st

from tables import attributes_to_table, draw_spectrum, peaks_to_table, render_spectrum

st.write("# Upload a spectrum library file!")

uploaded_file = st.file_uploader("Upload a spectrum library file", None)

spectrum: Spectrum = None
lib: SpectrumLibrary | None = None

if uploaded_file is None:
    pass
else:
    lib = SpectrumLibrary(filename=uploaded_file)

    st.write("## <mzSpecLib>")
    attributes_to_table(lib.attributes)

    entries = pd.DataFrame(
        [
            {"name": entry.name, "key": entry.number, "index": entry.index}
            for entry in lib.index
        ]
    )

    st.markdown("""
    Select a spectrum by clicking on the row margin
""")
    entries_table = st.dataframe(
        entries,
        selection_mode='single-row',
        on_select="rerun",
        hide_index=True,
    )

    if entries_table.selection['rows']:
        spectrum = lib[entries_table.selection["rows"][0]]
        render_spectrum(spectrum)