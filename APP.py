"""
🧬 DNA Sequence Analyzer & Mutation Detector
=============================================
Main Streamlit application for interactive DNA sequence analysis.
Upload FASTA files to analyze sequences, detect mutations, and
visualize results through an interactive dashboard.
"""

import streamlit as st
import pandas as pd
import os

from Sequence import parse_fasta, calculate_gc_content, calculate_base_composition, analyze_sequence
from mutation import *
from visualizer import *

# ─── Page Configuration ───────────────────────────────────────────────────────

st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA Sequence Analyzer & Mutation Detector")
st.markdown("Upload a FASTA file to analyze sequences, detect mutations, and visualize results.")

# ─── Sidebar Settings ─────────────────────────────────────────────────────────

st.sidebar.header("⚙️ Settings")
input_method = st.sidebar.radio(
    "Choose input method:",
    ["📁 Upload FASTA file", "📋 Paste FASTA sequence", "📂 Use sample data"],
    index=2
)

uploaded_file = None
pasted_fasta = ""

if input_method == "📁 Upload FASTA file":
    uploaded_file = st.sidebar.file_uploader("Upload FASTA file", type=["fasta", "fa", "txt"])
elif input_method == "📋 Paste FASTA sequence":
    pasted_fasta = st.sidebar.text_area(
        "Paste your FASTA sequence(s) below:",
        height=200,
        placeholder=">Sequence_ID\nATGCGTACGTAGCTAGC...\n\n>Another_Sequence\nGCTAGCTAGCATGCATG..."
    )
    st.sidebar.caption("💡 Tip: Copy directly from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/)")

# ─── Data Loading ─────────────────────────────────────────────────────────────

records = []
if uploaded_file:
    with open("temp.fasta", "wb") as f:
        f.write(uploaded_file.read())
    records = parse_fasta("temp.fasta")
    st.sidebar.success(f"✅ Loaded {len(records)} sequences from file")
elif pasted_fasta.strip():
    with open("temp_pasted.fasta", "w") as f:
        f.write(pasted_fasta)
    records = parse_fasta("temp_pasted.fasta")
    st.sidebar.success(f"✅ Loaded {len(records)} sequences from pasted input")
elif input_method == "📂 Use sample data":
    sample_path = "sample.fasta"
    if os.path.exists(sample_path):
        records = parse_fasta(sample_path)
        st.sidebar.info(f"📂 Using sample data: {len(records)} sequences")

# ─── Main Analysis ────────────────────────────────────────────────────────────

if records:
    results = [analyze_sequence(r) for r in records]

    # ─── Tab Layout ───────────────────────────────────────────────────────────

    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Overview", "🔬 Sequence Detail", "🧪 Mutation Analysis", "🔗 ORF Finder"
    ])

    # ─── Tab 1: Overview ──────────────────────────────────────────────────────

    with tab1:
        st.subheader("Sequence Summary")

        summary_data = [{
            "ID": r["id"],
            "Length (bp)": r["length"],
            "GC Content (%)": r["gc_content"],
            "ORFs Found": len(r["orfs"]),
            "Protein Length": len(r["protein"])
        } for r in results]

        st.dataframe(pd.DataFrame(summary_data), use_container_width=True)
        st.plotly_chart(plot_gc_comparisons(results), use_container_width=True)

    # ─── Tab 2: Sequence Detail ───────────────────────────────────────────────

    with tab2:
        selected_id = st.selectbox("Select Sequence", [r["id"] for r in results])
        selected = next(r for r in results if r["id"] == selected_id)

        col1, col2, col3 = st.columns(3)
        col1.metric("Length", f"{selected['length']} bp")
        col2.metric("GC Content", f"{selected['gc_content']}%")
        col3.metric("ORFs Found", len(selected["orfs"]))

        st.plotly_chart(
            plot_base_composition(selected["base_composition"], selected["id"]),
            use_container_width=True
        )

        with st.expander("🔤 Raw Sequence"):
            st.code(selected["sequence"], language=None)
        with st.expander("↔️ Reverse Complement"):
            st.code(selected["reverse_complement"], language=None)
        with st.expander("🧫 Translated Protein"):
            st.code(selected["protein"], language=None)

    # ─── Tab 3: Mutation Analysis ─────────────────────────────────────────────

    with tab3:
        st.subheader("Compare Two Sequences for Mutations")

        if len(results) >= 2:
            col1, col2 = st.columns(2)
            seq1_id = col1.selectbox("Reference Sequence", [r["id"] for r in results], index=0)
            seq2_id = col2.selectbox("Query Sequence", [r["id"] for r in results], index=1)

            seq1 = next(r for r in results if r["id"] == seq1_id)["sequence"]
            seq2 = next(r for r in results if r["id"] == seq2_id)["sequence"]

            mutations = detect_point_mutations(seq1, seq2)
            rate = calculate_mutation_rate(seq1, seq2, mutations)
            summary = summarize_mutations(mutations)

            st.metric("Mutation Rate", f"{rate}%")
            st.write("**Mutation Summary:**", summary)

            if mutations:
                mut_fig = plot_mutation_positions(mutations, max(len(seq1), len(seq2)))
                if mut_fig:
                    st.plotly_chart(mut_fig, use_container_width=True)

                mut_df = pd.DataFrame(mutations)
                st.dataframe(mut_df, use_container_width=True)
        else:
            st.info("Need at least 2 sequences to compare mutations.")

    # ─── Tab 4: ORF Finder ────────────────────────────────────────────────────

    with tab4:
        st.subheader("Open Reading Frames (ORFs)")
        selected_id_orf = st.selectbox(
            "Select Sequence for ORF Analysis",
            [r["id"] for r in results],
            key="orf_select"
        )
        selected_orf = next(r for r in results if r["id"] == selected_id_orf)

        orfs = selected_orf["orfs"]
        if orfs:
            fig = plot_orf_map(orfs, selected_orf["length"], selected_orf["id"])
            if fig:
                st.plotly_chart(fig, use_container_width=True)
                st.subheader(f"Found {len(orfs)} ORFs")
                orf_df = pd.DataFrame([{
                    "Frame": o["frame"],
                    "Start": o["start"],
                    "End": o["end"],
                    "Length (bp)": o["length"]
                } for o in orfs])
                st.dataframe(orf_df, use_container_width=True)
        else:
            st.info("No ORFs found in this sequence (minimum length: 30 bp)")

else:
    st.info("👆 Upload a FASTA file or enable 'Use sample data' to get started.")
