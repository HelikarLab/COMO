from pathlib import Path
from loguru import logger
import pandas as pd
import plotly.express as px

__all__ = ["z_score_distribution"]


def z_score_distribution(
    df: pd.DataFrame,
    title: str,
    output_png_filepath: Path,
):
    if not output_png_filepath.suffix == ".png":
        logger.warning(
            f"Expected .png suffix for output_png_filepath, got {output_png_filepath.suffix}. Defaulting to .png"
        )
        output_png_filepath = output_png_filepath.with_suffix(".png")

    fig = px.histogram(
        df,
        x="zscore",
        color="source",
        nbins=100,
        marginal="rug",
        title=title,
    )

    fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font={"family": "sans-serif", "size": 12})

    # Simplified plot for many sources (optional)
    if len(df["source"].unique()) > 10:
        fig.update_layout(showlegend=False)

    fig.write_image(output_png_filepath)
