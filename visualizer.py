"""
Visualization Module
====================
Interactive Plotly chart functions for visualizing DNA sequence analysis
results including base composition, GC content comparison, mutation
positions, and ORF maps.
"""

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd


def plot_base_composition(base_composition: dict, seq_id: str):
    """Create a color-coded bar chart of nucleotide base composition.

    Each base (A, T, G, C) is shown with a distinct color and percentage labels.

    Args:
        base_composition: Dict mapping bases to count/percentage dicts.
        seq_id: Sequence identifier for the chart title.

    Returns:
        Plotly Figure object.
    """
    bases = list(base_composition.keys())
    counts = [base_composition[b]["count"] for b in bases]
    percentages = [base_composition[b]["percentage"] for b in bases]

    colors = {"A": "#2ECC71", "T": "#E74C3C", "G": "#3498DB", "C": "#F39C12"}
    bar_colors = [colors[b] for b in bases]

    fig = go.Figure(data=[
        go.Bar(
            x=bases,
            y=counts,
            marker_color=bar_colors,
            text=[f"{p}%" for p in percentages],
            textposition="auto"
        )
    ])

    fig.update_layout(
        title=f"Base Composition — {seq_id}",
        xaxis_title="Nucleotide Base",
        yaxis_title="Count",
        template="plotly_dark"
    )
    return fig


def plot_gc_comparisons(results: list[dict]):
    """Create a bar chart comparing GC content across multiple sequences.

    Includes a dashed red line at 50% as a reference threshold.

    Args:
        results: List of analysis result dicts with 'id' and 'gc_content'.

    Returns:
        Plotly Figure object.
    """
    ids = [r["id"] for r in results]
    gc_values = [r["gc_content"] for r in results]

    fig = px.bar(
        x=ids,
        y=gc_values,
        labels={"x": "Sequence ID", "y": "GC Content (%)"},
        title="GC Content Comparison Across Sequences",
        color=gc_values,
        color_continuous_scale="viridis"
    )
    fig.add_hline(
        y=50, line_dash="dash",
        line_color="red", annotation_text="50% threshold"
    )
    fig.update_layout(template="plotly_dark")
    return fig


def plot_mutation_positions(mutations: list[dict], seq_length: int):
    """Create a scatter plot of mutation positions along the sequence.

    Mutations are color-coded by type (transition vs transversion).

    Args:
        mutations: List of mutation dicts from detect_point_mutations().
        seq_length: Length of the reference sequence for x-axis range.

    Returns:
        Plotly Figure object, or None if no mutations.
    """
    if not mutations:
        return None

    df = pd.DataFrame(mutations)
    fig = px.scatter(
        df,
        x="position",
        y="type",
        color="type",
        hover_data=["original", "mutated"],
        title="Mutation Position Map",
        labels={"position": "Position in Sequence", "type": "Mutation Type"}
    )
    fig.update_traces(marker=dict(size=12))
    fig.update_layout(template="plotly_dark", xaxis_range=[0, seq_length])
    return fig


def plot_orf_map(orfs: list[dict], seq_length: int, seq_id: str):
    """Create a horizontal bar chart showing ORF positions along the sequence.

    Each ORF is displayed as a colored bar positioned at its start location.

    Args:
        orfs: List of ORF dicts with start, end, length, frame keys.
        seq_length: Total sequence length for x-axis range.
        seq_id: Sequence identifier for the chart title.

    Returns:
        Plotly Figure object, or None if no ORFs.
    """
    if not orfs:
        return None

    fig = go.Figure()
    colors = ["#E74C3C", "#3498DB", "#2ECC71"]

    for i, orf in enumerate(orfs):
        fig.add_trace(go.Bar(
            x=[orf["length"]],
            y=[f"ORF {i+1} (Frame {orf['frame']})"],
            base=[orf["start"]],
            orientation="h",
            marker_color=colors[orf["frame"] - 1],
            name=f"ORF {i+1}",
            hovertemplate=(
                f"Start: {orf['start']}<br>"
                f"End: {orf['end']}<br>"
                f"Length: {orf['length']} bp"
            )
        ))

    fig.update_layout(
        title=f"ORF Map — {seq_id}",
        xaxis_title="Position in Sequence",
        xaxis_range=[0, seq_length],
        template="plotly_dark",
        barmode="overlay"
    )
    return fig