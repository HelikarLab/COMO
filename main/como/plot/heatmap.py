from __future__ import annotations

import math
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FixedLocator


def condition_vs_pathway(
    data: pd.DataFrame,
    save_filepath: Path | None = None,
    *,
    copy_df: bool = False,
    exclude_zero_flux_pathways: bool = False,
) -> plt.Figure:
    """Build a heatmap of fluxes through pathways across conditions.

    :param data: Index values are conditions, column names are pathways
    :param save_filepath: If provided, the resulting figure will be saved to this location
    :param copy_df: Should the incoming dataframe be copied to prevent modifications to data?
    :param exclude_zero_flux_pathways: Should pathways that have 0 flux across all rows be excluded?

    :return: The resulting `matpotlib.pyplt.Figure` object
    """
    plot_df: pd.DataFrame = data.copy() if copy_df else data
    plot_df = plot_df.astype(np.float32)
    fig = plt.figure(figsize=(100, 40), dpi=175)

    if exclude_zero_flux_pathways:
        # Select pathways that have at least one non-zero value
        plot_df = plot_df.loc[:, plot_df.where(plot_df != 0).any(axis=0)]

    # Identify the second largest (pos or neg) value
    # This is needed in order to set the upper & lower bounds for the graph, excluding +/- 1000-flux values
    plot_df[plot_df > 1000] = 1001
    plot_df[plot_df < -1000] = -1001
    second_largest_positive = plot_df[plot_df > 0].stack().drop_duplicates().nlargest(2).iloc[-1]
    second_largest_negative = plot_df[plot_df < 0].stack().drop_duplicates().nsmallest(2).iloc[-1]
    vmax = max(abs(second_largest_negative), second_largest_positive)

    # Convert tick marks to reasonable values:
    #   max tick < 100: round to 10s place
    #   max tick < 1_000: round to 100s place
    #   max tick < 10_000: round to 1_000s place
    base = 10 if vmax < 100 else 100 if vmax < 10_000 else 1000
    vmax_root = math.ceil(vmax / base) * base

    # Create 5 evenly spaced ticks along the legend
    ticks = np.linspace(-vmax_root, vmax_root, 5)

    # Generate legend gradient
    norm = mcolors.TwoSlopeNorm(vmin=-vmax_root, vcenter=0, vmax=vmax_root)

    # If a value falls outside of `vmax_root`, set it to the following colors
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_over("#660033")
    cmap.set_under("#000099")

    ax: plt.Axes = sns.heatmap(
        data=plot_df,
        linewidths=1.0,
        linecolor="#686868",
        center=0,
        yticklabels=True,
        xticklabels=True,
        norm=norm,
        cmap=cmap,
        cbar_kws={"extend": "both", "label": f"Flux ratio (clipped at ±{vmax:.0f})"},
    )

    plt.title("Metabolic Model Flux Sum through Pathways", fontsize=100)

    plt.xlabel("Pathway", fontsize=80)
    ax.tick_params(axis="x", which="major", labelsize=55, labelrotation=90)

    plt.ylabel("Condition", fontsize=85)
    ax.tick_params(axis="y", which="major", labelsize=55, labelrotation=0)

    cbar = ax.collections[0].colorbar
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_major_locator(FixedLocator(ticks))
    cbar.update_ticks()

    # Add extended triangles that are detached from the colorbar (prevents very large pos/neg values from blowing out the legend)
    cbar.ax.text(0.5, 1.06, "> +1000", ha="center", va="bottom", transform=cbar.ax.transAxes, fontsize=40)
    cbar.ax.text(0.5, -0.06, "< -1000", ha="center", va="top", transform=cbar.ax.transAxes, fontsize=40)
    cbar.ax.tick_params(labelsize=40)
    cbar.set_label("Flux", rotation=270, labelpad=40)

    fig.tight_layout(h_pad=0.85)

    if save_filepath:
        save_filepath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_filepath, transparent=False, bbox_inches="tight")

    return fig
