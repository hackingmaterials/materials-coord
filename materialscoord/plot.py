from typing import Dict, Optional
from pandas import DataFrame

import seaborn as sns
import matplotlib.pyplot as plt


def plot_benchmark_scores(scores: DataFrame,
                          structure_mapping: Optional[Dict[str, str]] = None,
                          nn_method_mapping: Optional[Dict[str, str]] = None,
                          figsize=(17, 10)
                          ) -> plt:
    """
    Plot MaterialsCoord benchmark scores.

    Args:
        scores: The scores dataframe from `Benchmark.scores`.
        structure_mapping: A dictionary to remap the structure labels. I.e., to rename
            "NaCl_rocksalt_100633" to "NaCl (rock-salt)", structure_mapping should be::

                {"NaCl_rocksalt_100633": "NaCl (rock-salt)"}
        nn_method_mapping: A dictionary to remap the near neighbor labels. I.e., to rename
            "BrunnerNN_reciprocal" to "BrunnerNN", nn_method_remapping should be::

                {"BrunnerNN_reciprocal": "BrunnerNN"}

    Returns:
        A matplotlib pyplot object.
    """
    structure_mapping = structure_mapping or {}
    nn_method_mapping = nn_method_mapping or {}

    structure_labels = []
    for structure_label in scores.index:
        structure_labels.append(structure_mapping.get(structure_label, structure_label))

    nn_labels = []
    for nn_label in scores.columns:
        nn_labels.append(nn_method_mapping.get(nn_label, nn_label))

    sns.set(font='Helvetica')
    sns.set(font_scale=1.5)
    sns.set_style({'axes.edgecolor': 'black', 'axes.linewidth': 1.3})

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(scores, annot=True, cmap="Blues", vmax=10)

    ax.set_yticklabels(structure_labels)
    ax.set_xticklabels(nn_labels, rotation=60)

    # draw line above Total row and axis boundaries
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    ax.axhline(y=len(structure_labels)-1, color='k', linewidth=1.3)

    cbar = ax.collections[0].colorbar
    cbar.set_label('Benchmark score', labelpad=40, rotation=270)
    cbar.outline.set(**{"edgecolor": 'k', "linewidth": 1.3})

    return plt
