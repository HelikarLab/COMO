import colorsys
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from project import Config
from scipy import stats
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity


@dataclass
class ZResult:
    zfpkm: pd.DataFrame
    x_range: pd.DataFrame  # npt.NDArray[np.float64]
    density: pd.DataFrame  # npt.NDArray[np.float64]
    mu: npt.NDArray[np.float64]
    std_dev: npt.NDArray[np.float64]
    max_fpkm_peak: npt.NDArray[np.float64]


def z_score_calc(abundance: pd.DataFrame, min_thresh: int) -> ZResult:
    values = abundance.values.copy() + 1
    log_abundance_filt = np.log2(values[values > min_thresh]).reshape((len(abundance), len(abundance.columns)))
    log_abundance = np.log2(values)

    # np.zeros((1000, len(abundance.columns)), dtype=np.float64),
    z_result = ZResult(
        zfpkm=pd.DataFrame(
            data=np.nan * np.ones_like(values), index=abundance.index, columns=abundance.columns, dtype=np.float64
        ),
        x_range=pd.DataFrame(
            data=np.zeros((1000, len(abundance.columns))), columns=abundance.columns, dtype=np.float64
        ),
        density=pd.DataFrame(
            data=np.zeros((1000, len(abundance.columns))), columns=abundance.columns, dtype=np.float64
        ),
        mu=np.zeros(len(abundance.columns)),
        std_dev=np.zeros(len(abundance.columns)),
        max_fpkm_peak=np.zeros(len(abundance.columns)),
    )

    for i, col in enumerate(abundance.columns):
        kde: KernelDensity = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(
            log_abundance_filt[:, i].reshape(-1, 1)
        )  # type: ignore
        x_range = np.linspace(log_abundance[:, i].min(), log_abundance[:, i].max(), 1000)
        density = np.exp(kde.score_samples(x_range.reshape(-1, 1)))
        peaks, _ = find_peaks(density, height=0.02, distance=1.0)
        peak_positions = x_range[peaks]

        mu = peak_positions.max()
        max_fpkm_peak = density[peaks[np.argmax(peak_positions)]]

        # Select rows from `log_abundance` that are greater than 0 and less than mu in the column `i`
        u = log_abundance[:, i][(log_abundance[:, i] > 0) & (log_abundance[:, i] < mu)].mean()
        std_dev = (mu - u) * np.sqrt(np.pi / 2)
        zfpkm_values = (log_abundance[:, i] - mu) / std_dev

        z_result.zfpkm[col] = zfpkm_values
        z_result.x_range[col] = x_range
        z_result.density[col] = density
        z_result.mu[i] = mu
        z_result.std_dev[i] = std_dev
        z_result.max_fpkm_peak[i] = max_fpkm_peak

    return z_result


def lighten_color(red: int, green: int, blue: int, factor: float = 0.5) -> str:
    # Convert RGB to HLS
    hue, lightness, saturation = colorsys.rgb_to_hls(red / 255.0, green / 255.0, blue / 255.0)

    # Increase lightness
    lightness = min(1.0, lightness + (1 - lightness) * factor)

    # Covnert back to RGB values
    r, g, b = colorsys.hsv_to_rgb(hue, saturation, lightness)
    return f"rgb({int(r * 255)},{int(g * 255)},{int(b * 255)})"


# Plotting function
def plot_gaussian_fit(z_results: ZResult, facet_titles: bool = True, x_min: int = -4) -> go.Figure:
    zfpkm = z_results.zfpkm
    x_range = z_results.x_range
    density = z_results.density
    mu = z_results.mu
    std_dev = z_results.std_dev
    max_fpkm = z_results.max_fpkm_peak

    fig = make_subplots(
        rows=len(zfpkm.columns), cols=1, subplot_titles=zfpkm.columns if facet_titles else [None] * len(zfpkm.columns)
    )
    for i, col in enumerate(zfpkm.columns):
        fitted = stats.norm.pdf(x_range[col], loc=mu[i], scale=std_dev[i])
        scale_fit = fitted * (max_fpkm[i] / fitted.max())

        # Select random RGB values for each plot
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)
        color, lighten = f"rgb({r},{g},{b})", lighten_color(r, g, b)
        fig.add_trace(
            go.Scatter(x=x_range[col], y=density[col], name="Abundance Density", line=dict(color=color)),
            row=i + 1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(x=x_range[col], y=scale_fit, name="Fitted Gaussian", line=dict(dash="dash", color=lighten)),
            row=i + 1,
            col=1,
        )

        fig.update_xaxes(title_text="log2(abundance)", range=[x_min, x_range[col].max()], row=i + 1, col=1)
        fig.update_yaxes(title_text="[scaled] density", row=i + 1, col=1)

    fig.update_layout(height=len(zfpkm.columns) * 250, width=800, title_text="Gaussian Fit per Sample")
    return fig


# Main function for protein abundance transformation
def protein_transform_main(
    abundance_df: Union[pd.DataFrame, str, Path], out_dir: Union[str, Path], group_name: str, config: Config
) -> None:
    output_figure_directory = config.get_output_path(out_dir) / "figures"
    output_figure_directory.mkdir(parents=True, exist_ok=True)

    abundance_df: pd.DataFrame = (
        pd.read_csv(abundance_df) if isinstance(abundance_df, (str, Path)) else abundance_df.fillna(0)
    )
    abundance_df = abundance_df[np.isfinite(abundance_df).all(axis=1)]  # Remove +/- infinity values
    z_transform: ZResult = z_score_calc(abundance_df, min_thresh=0)

    fig = plot_gaussian_fit(z_results=z_transform, facet_titles=True, x_min=-4)
    fig.write_image(config.get_output_path(out_dir) / "gaussian_fit.png")
    fig.write_html(config.get_output_path(out_dir) / "gaussian_fit.html")
    print(f"Wrote image to {config.get_output_path(out_dir) / 'gaussian_fit.png'}")

    z_transformed_abundances = z_transform.zfpkm
    z_transformed_abundances[abundance_df == 0] = -4
    z_transformed_abundances.to_csv(
        config.get_output_path(out_dir) / f"protein_zscore_Matrix_{group_name}.csv", index=False
    )


def process_proteomics_data(context_name: str, config: Config) -> None:
    """Process proteomics data for a given context"""
    # Get paths using Config methods
    context_path = config.get_context_path(context_name)
    proteomics_path = context_path / "proteomics"
    proteomics_path.mkdir(parents=True, exist_ok=True)

    # Initialize R interface with config
    r_script = config.code_dir / "rscripts" / "protein_transform.R"
    protein_transform = rpy2_api.Rpy2(r_file_path=r_script, config=config)

    # Process abundance matrix
    abundance_matrix = config.get_matrix_path(context_name, f"protein_abundance_matrix_{context_name}.csv")
    if not abundance_matrix.exists():
        logger.error(f"Abundance matrix not found at {abundance_matrix}")
        return

    # Call R function with proper paths
    protein_transform.call_function("protein_transform_main", str(abundance_matrix), str(proteomics_path), context_name)

    logger.success(f"Proteomics data processed for {context_name}")
