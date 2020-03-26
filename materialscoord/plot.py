"""Define materialscoord plotting functions."""

from typing import Dict, Optional, Tuple
from pandas import DataFrame

import seaborn as sns
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_benchmark_scores(
    scores: DataFrame,
    structure_mapping: Optional[Dict[str, str]] = None,
    nn_method_mapping: Optional[Dict[str, str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    vmax: Optional[float] = None,
    vmin: Optional[float] = None,
    round_dp: Optional[int] = 1,
    cbar_label: str = "Benchmark score",
) -> plt:
    """
    Plot MaterialsCoord benchmark scores as a heatmap.

    Args:
        scores: The scores dataframe from `Benchmark.scores`.
        structure_mapping: A dictionary to remap the structure labels. I.e., to rename
            "NaCl_rocksalt_100633" to "NaCl (rock-salt)", structure_mapping should be::

                {"NaCl_rocksalt_100633": "NaCl (rock-salt)"}
        nn_method_mapping: A dictionary to remap the near neighbor labels. I.e., to
            rename "BrunnerNN_reciprocal" to "BrunnerNN", nn_method_remapping should
            be::

                {"BrunnerNN_reciprocal": "BrunnerNN"}
        figsize: The figure size as a tuple of (width, height). If not provided this
            will be automatically estimated based on the number of structures and
            near neighbor methods.
        vmax: The maximum value used for scaling the colorbar of the heatmap.
        vmin: The minimum value used for scaling the colorbar of the heatmap.
        round_dp: Number of places to round the results to. If None, no rounding will
            occur. Note that seaborn may automatically round the numbers anyway.
        cbar_label: The y-axis label for the colorbar.

    Returns:
        A matplotlib pyplot object.
    """
    if round_dp is not None:
        scores = scores.round(round_dp)

    if not figsize:
        width = 4 + 1.5 * len(scores.columns)
        height = 2 + 0.45 * len(scores.index)
        figsize = (width, height)

    vmax = vmax if vmax is not None else scores.iloc[:-1].max().max()
    vmin = vmin if vmin is not None else 0

    structure_mapping = structure_mapping or {}
    nn_method_mapping = nn_method_mapping or {}

    scores = scores.rename(columns=nn_method_mapping)
    scores = scores.rename(structure_mapping)

    sns.set(
        font=["Helvetica", "Arial"],
        font_scale=1.5,
        rc={
            "figure.figsize": figsize,
            "axes.edgecolor": "black",
            "axes.linewidth": 1.3,
            'pdf.fonttype': 42,
        },
    )

    # want the colorbar to have a fixed width regardless of the number of structures or
    # near neighbor methods
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=0.3, pad=0.5)

    ax = sns.heatmap(
        scores,
        annot=True,
        cmap="Blues",
        vmax=vmax,
        vmin=vmin,
        fmt="g",
        ax=ax,
        cbar_ax=cax,
    )
    ax.set_xticklabels(scores.columns, rotation=60)

    # draw line above Total row and axis boundaries
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    ax.axhline(y=len(scores) - 1, color="k", linewidth=1.3)

    cbar = ax.collections[0].colorbar
    cbar.set_label(cbar_label, labelpad=40, rotation=270)
    cbar.outline.set(edgecolor="k", linewidth=1.3)

    return plt
