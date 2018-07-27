import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt
import seaborn as sns

def heatmap(df, inline=False, figsize=(20, 12), fontsize=16, labels=True):

    #cmap = sns.cubehelix_palette(light=0.95, as_cmap=True)
    cmap = "Viridis_r"
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    ax0 = plt.subplot(gs[0])

    if labels:
        sns.heatmap(df, cmap=cmap, ax=ax0, cbar=True, vmax=10,
                    cbar_kws=dict(use_gridspec=False, location="left"),
                    annot=True, annot_kws={'size': fontsize - 2})
    else:
        sns.heatmap(df, cmap=cmap, ax=ax0, cbar=False)

    # Apply visual corrections and modifications.
    #ax0.xaxis.tick_bottom()
    ax0.yaxis.set_ticks_position('right')
    ax0.set_xticklabels(ax0.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=fontsize)
    ax0.set_yticklabels(ax0.yaxis.get_majorticklabels(), fontsize=fontsize, rotation=0)
    ax0.set_yticklabels(ax0.yaxis.get_majorticklabels(), rotation=0, fontsize=fontsize)
    ax0.xaxis.set_ticks_position('none')
    ax0.yaxis.set_ticks_position('none')
    ax0.patch.set_visible(False)

    for text in ax0.texts:
        t = float(text.get_text())
        if 0.95 <= t < 1:
            text.set_text('<1')
        elif -1 < t <= -0.95:
            text.set_text('>-1')
        elif t == 1:
            text.set_text('1')
        elif t == -1:
            text.set_text('-1')
        elif -0.05 < t < 0.05:
            text.set_text('')
        else:
            text.set_text(round(t, 1))

    if inline:
        plt.show()
    else:
        return ax0