import pandas as pd

import matplotlib

matplotlib.use("agg")
from matplotlib import pyplot as plt

from mzspeclib import SpectrumLibrary, Spectrum, draw
from mzspeclib.attributes import Attributed

import streamlit as st

spectrum: Spectrum = None
lib: SpectrumLibrary | None = None


def attributes_to_table(attributed: Attributed):
    attr_table = []
    for attr in attributed:
        attr_table.append({"group": attr.group_id or "", "name": attr.key, "value": str(attr.value)})
    attr_table = pd.DataFrame(attr_table).set_index("group")
    st.dataframe(attr_table, use_container_width=True)


def peaks_to_table(spectrum: Spectrum):
    rows = []
    for mz, intens, annots, aggs in spectrum.peak_list:
        row = {
            "m/z": mz,
            "intensity": intens,
            "annotations": ",".join(map(str, annots or [])),
            "aggregations": "\t".join(map(str, aggs or [])),
        }
        rows.append(row)
    st.write("### <Peaks>")
    st.dataframe(pd.DataFrame(rows), use_container_width=True)


def draw_spectrum(spectrum: Spectrum, **kwargs):
    fig, ax = plt.subplots(1, 1)
    draw.draw_spectrum(spectrum, ax=ax, **kwargs)
    return fig


def render_spectrum(spectrum: Spectrum):
    fig = draw_spectrum(spectrum)
    st.pyplot(fig, True, use_container_width=True)

    st.write(f"### <Spectrum {spectrum.key}>")
    attributes_to_table(spectrum.attributes)
    for analyte in spectrum.analytes.values():
        st.write(f"### <Analyte {analyte.id}>")
        attributes_to_table(analyte)
    for interpretation in spectrum.interpretations.values():
        if interpretation:
            st.write(f"### <Interpretation {interpretation.id}>")
            attributes_to_table(interpretation.attributes)
            for member_interp in interpretation.member_interpretations.values():
                st.write(f"### <InterpretationMember {member_interp.id}>")
                attributes_to_table(member_interp.attributes)
    peaks_to_table(spectrum)