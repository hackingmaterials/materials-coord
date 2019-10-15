from typing import Dict, Optional, Tuple
from pandas import DataFrame

import seaborn as sns
import matplotlib.pyplot as plt


def plot_benchmark_scores(
    scores: DataFrame,
    structure_mapping: Optional[Dict[str, str]] = None,
    nn_method_mapping: Optional[Dict[str, str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    vmax: Optional[float] = None,
    round_dp: Optional[int] = 1,
) -> plt:
    """
    Plot MaterialsCoord benchmark scores as a heatmap.

    Args:
        scores: The scores dataframe from `Benchmark.scores`.
        structure_mapping: A dictionary to remap the structure labels. I.e., to rename
            "NaCl_rocksalt_100633" to "NaCl (rock-salt)", structure_mapping should be::

                {"NaCl_rocksalt_100633": "NaCl (rock-salt)"}
        nn_method_mapping: A dictionary to remap the near neighbor labels. I.e., to rename
            "BrunnerNN_reciprocal" to "BrunnerNN", nn_method_remapping should be::

                {"BrunnerNN_reciprocal": "BrunnerNN"}
        figsize: The figure size as a tuple of (width, height). If not provided this
            will be automatically estimated based on the number of structures and
            near neighbor methods.
        vmax: The maximum value used for scaling the colorbar of the heatmap.
        round_dp: Number of places to round the results to. If None, no rounding will
            occur. Note that seaborn may automatically round the numbers anyway.

    Returns:
        A matplotlib pyplot object.
    """
    if round_dp:
        scores.round(round_dp)

    if not figsize:
        width = 4 + 1.5 * len(scores.columns)
        height = 2 + 0.45 * len(scores.index)
        figsize = (width, height)

    if not vmax:
        vmax = scores.iloc[:-1].values.max()

    structure_mapping = structure_mapping or {}
    nn_method_mapping = nn_method_mapping or {}

    structure_labels = []
    for structure_label in scores.index:
        structure_labels.append(structure_mapping.get(structure_label, structure_label))

    nn_labels = []
    for nn_label in scores.columns:
        nn_labels.append(nn_method_mapping.get(nn_label, nn_label))

    sns.set(font="Helvetica")
    sns.set(font_scale=1.5)
    sns.set_style({"axes.edgecolor": "black", "axes.linewidth": 1.3})

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(scores, annot=True, cmap="Blues", vmax=vmax)

    ax.set_yticklabels(structure_labels)
    ax.set_xticklabels(nn_labels, rotation=60)

    # draw line above Total row and axis boundaries
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    ax.axhline(y=len(structure_labels) - 1, color="k", linewidth=1.3)

    cbar = ax.collections[0].colorbar
    cbar.set_label("Benchmark score", labelpad=40, rotation=270)
    cbar.outline.set(**{"edgecolor": "k", "linewidth": 1.3})

    return plt
