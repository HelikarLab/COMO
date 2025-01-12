from pathlib import Path
from loguru import logger
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


__all__ = ["z_score_distribution"]


def z_score_distribution(
    df: pd.DataFrame,
    title: str,
    output_filepath: Path,
):
    if output_filepath.suffix not in {".png", ".pdf", ".svg"}:
        logger.warning(
            f"Expected .png, .pdf, or .svg suffix for output_png_filepath, got {output_filepath.suffix}. Defaulting to .pdf"
        )
        output_filepath = output_filepath.with_suffix(".pdf")
    logger.trace(f"Graphing z-score distribution")
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    output_filepath.unlink(missing_ok=True)

    plt.figure(figsize=(10, 6))

    if len(df["source"].unique()) == 1:
        ax = sns.histplot(df, x="zscore", bins=100, kde=True)
        sns.rugplot(df, x="zscore", ax=ax)
    else:
        sns.histplot(df, x="zscore", hue="source", bins=100, kde=True, element="step")
        plt.legend(loc="upper right", frameon=False, title=None)

    plt.title(title)
    plt.xlabel("Z-score")
    plt.ylabel("Frequency")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_filepath)
    plt.close()
    logger.success(f"Saved z-score distribution graph to '{output_filepath}'")
