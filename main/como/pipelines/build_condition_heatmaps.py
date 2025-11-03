"""Create heatmaps of conditions vs pathway flux.

This pipeline will generate heatmaps showing the flux through various pathways
"""

from __future__ import annotations

import concurrent.futures
from collections import defaultdict
from functools import partial
from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from pandas._libs.missing import NAType

from como.plot.heatmap import condition_vs_pathway


def find_possible_model_filepaths(search_dir: Path) -> list[Path]:
    """Find potential files that could be constraint-based metabolic models.

    Args:
        search_dir: The directory to search for models.

    Returns:
        Potential filepaths that could be loaded as a `cobra.Model` object
    """
    return [f for f in search_dir.rglob("*") if f.suffix in {".mat", ".json", ".sbml", ".xml"}]


def get_cobra_model_if_valid(filepath: Path) -> cobra.Model | None:
    """Evaluate if a given filepath can be read as a `cobra.Model`.

    Args:
        filepath: The filepath to read

    Returns:
        a `cobra.Model` object if the file can be read, otherwise None
    """
    if filepath.suffix == ".json":
        return cobra.io.load_json_model(filepath)
    elif filepath.suffix == ".mat":
        return cobra.io.load_matlab_model(filepath)
    elif filepath.suffix in {".yaml", ".yml"}:
        return cobra.io.load_yaml_model(filepath)
    elif filepath.suffix in {".xml", ".sbml"}:
        return cobra.io.read_sbml_model(filepath)
    return None


def get_model_flux(model: cobra.Model, objective: str = "biomass_maintenance", solver: str = "gurobi") -> pd.Series:
    """Get the flux through a CBMM.

    Args:
        model: A `cobra.Model` object
        objective: The objective function to optimize
        solver: The solver to use

    Returns:
        A pandas Series of reaction fluxes indexed by reaction ID
    """
    model.objective = objective
    model.solver = solver
    solution = model.optimize()
    return solution.fluxes


def get_many_model_flux(
    models: list[cobra.Model],
    objective: str = "biomass_maintenance",
    solver: str = "gurobi",
    cores: int = 4,
    process_pool: concurrent.futures.ProcessPoolExecutor | None = None,
    colnames: list[str] | None = None,
    na_value: NAType | int | float = NAType,
) -> pd.DataFrame:
    """Get the flux through many CBMMs.

    Args:
        models: A list of `cobra.Model` objects
        objective: The objective function to optimize
        solver: The solver to use
        cores: The number of CPU cores to use
        process_pool: An existing process pool to use
        colnames: Column names to use for the resulting dataframe
        na_value: Value to use for missing values in the dataframe

    Returns:
        A pandas DataFrame of reaction fluxes indexed by condition (row) and reaction ID (column)

    Raises:
        ValueError: If `colnames` is provided and its length does not match the number of models
    """
    if colnames and len(colnames) != len(models):
        raise ValueError("Length of colnames must match length of models")

    pool = process_pool or concurrent.futures.ProcessPoolExecutor(max_workers=cores)
    shutdown = not process_pool  # if the user provided a pool, do not shut it down

    func = partial(get_model_flux, objective=objective, solver=solver)
    series: list[pd.Series] = list(pool.map(func, models))
    for i, series_obj in enumerate(series):
        series_obj.name = colnames[i] if colnames else f"model_{i}"
    df: pd.DataFrame = pd.concat(list(series), axis="columns")

    if shutdown:
        pool.shutdown(wait=True)

    if na_value != NAType:  # no need to replace values that are already pd.NA
        df = df.fillna(na_value)

    df = df.T
    return df


def group_reactions_by_pathway(models: cobra.Model | list[cobra.Model], flux_df: pd.DataFrame) -> pd.DataFrame:
    """Group reactions by their subsystem/pathway and sum the fluxes.

    Args:
        models: A cobra.Model or list of cobra.Models
        flux_df: A dataframe of reaction fluxes, indexed by condition and with reaction IDs as columns

    Returns:
        A dataframe of pathway fluxes, indexed by condition and with pathways as columns
    """
    pathways_by_reaction: dict[str, set[str]] = defaultdict(set)
    models = [models] if isinstance(models, cobra.Model) else models

    for model in models:
        for reaction in model.reactions:
            reaction: cobra.Reaction
            pathways_by_reaction[reaction.subsystem].add(reaction.id)
    pathways_by_reaction.pop("", None)  # remove the empty pathway; faster than checking every reaction's subsystem

    # pathway_flux: pd.DataFrame = pd.DataFrame(index=flux_df.index, columns=list(pathways_by_reaction.keys()))
    # for condition in flux_df.index:
    #     for pathway, reactions in pathways_by_reaction.items():
    #         pathway_flux.loc[condition, pathway] = flux_df.loc[condition, list(reactions)].sum()
    pathway_fluxes: dict[str, pd.Series[npt.NDArray[float]]] = {}
    for pathway, reactions in pathways_by_reaction.items():
        reactions_in_df = list(reactions.intersection(flux_df.columns))
        if reactions_in_df:
            pathway_fluxes[pathway] = flux_df[reactions_in_df].sum(axis=1)
    return pd.DataFrame(pathway_fluxes)


def build_condition_vs_pathway_heatmap(
    data: pd.DataFrame | list[cobra.Model] | Path | list[Path],
    save_filepath: Path | None = None,
    objective: str = "biomass_maintenance",
    solver: str = "gurobi",
    process_pool: concurrent.futures.ProcessPoolExecutor | None = None,
    cores: int = 4,
    condition_names: list[str] | None = None,
    na_value: NAType | int | float = NAType,
    *,
    search_path: bool = False,
    copy_df_when_building_plot: bool = False,
    exclude_zero_flux_pathways: bool = False,
) -> plt.Figure:
    """Create a heatmap of conditions vs flux through pathways.

    If `data` is a pd.DataFrame:
        - The index names wile used as conditions and placed on the Y-axis
        - The column names will be used as pathways and placed on the X-axis. The columns should indicate pathways.

    If `data` is a Path and `search_path` is True:
        - Models will be recursively discovered under the given path
        - Models will be simulated with the given objective and  solver
        - A dataframe will be built from the resulting series based on the above rules

    If `data` is a list of Paths:
        - Models will be read and simulated for each path
        - A pd.DataFrame will be built from the resulting series based on the above rules

    If `data` is a list of cobra.Models:
        - Models will be simulated with the given objective and solver
        - A pd.DataFrame will be built from the resulting series based on the above rules

    Args:
        data: The data to use for the heatmap
        search_path: Whether to search the given path for models
        save_filepath: The filepath to save the heatmap to
        objective: The objective function to optimize
        solver: The solver to use
        process_pool: An existing process pool to use
        cores: The number of CPU cores to use
        condition_names: Column names to use for the resulting dataframe if `data` is a Path or list of Paths
        na_value: Value to use for missing values in the flux dataframe
        copy_df_when_building_plot: Whether to copy the dataframe when building the plot.
            This can be useful if the dataframe is going to be reused later.
        exclude_zero_flux_pathways: Whether to exclude pathways that have zero flux across all conditions

    Returns:
        A matplotlib Figure object containing the heatmap

    Raises:
        ValueError: If `search_path` is True and `data` is not a Path
    """
    if not isinstance(data, Path) and search_path:
        raise ValueError("If search_path is True, data must be a Path")

    flux_df: pd.DataFrame
    if isinstance(data, pd.DataFrame):
        return condition_vs_pathway(
            data,
            save_filepath=save_filepath,
            copy_df=copy_df_when_building_plot,
            exclude_zero_flux_pathways=exclude_zero_flux_pathways,
        )
    elif isinstance(data, list) and isinstance(data[0], cobra.Model):
        models = data
        flux_df = get_many_model_flux(
            models=data,
            objective=objective,
            solver=solver,
            cores=cores,
            process_pool=process_pool,
            colnames=condition_names,
            na_value=na_value,
        )
    elif isinstance(data, Path):
        if search_path:
            possible_model_fps: list[Path] = find_possible_model_filepaths(data)
            models = []
            for fp in possible_model_fps:
                if isinstance(model := get_cobra_model_if_valid(fp), cobra.Model):
                    models.append(model)
            flux_df = get_many_model_flux(
                models,
                objective=objective,
                solver=solver,
                cores=cores,
                process_pool=process_pool,
                colnames=condition_names,
                na_value=na_value,
            )
        else:
            models = get_cobra_model_if_valid(data)
            flux_df = pd.DataFrame(get_model_flux(models, objective=objective, solver=solver))
    elif isinstance(data, list) and isinstance(data[0], Path):
        models = [get_cobra_model_if_valid(fp) for fp in data]
        flux_df = get_many_model_flux(
            models,
            objective=objective,
            solver=solver,
            cores=cores,
            process_pool=process_pool,
            colnames=condition_names,
            na_value=na_value,
        )

    flux_df = group_reactions_by_pathway(models=models, flux_df=flux_df)
    return condition_vs_pathway(data=flux_df, save_filepath=save_filepath)


def _main():
    models = [
        Path("/home/joshl/projects/ImmunoMetabolism/data/model_build/A/A_pDCs/A_pDCs_model_imat.json"),
        Path("/home/joshl/projects/ImmunoMetabolism/data/model_build/B/B_pDCs/B_pDCs_model_imat.json"),
        Path("/home/joshl/projects/ImmunoMetabolism/data/model_build/C/C_pDCs/C_pDCs_model_imat.json"),
        Path("/home/joshl/projects/ImmunoMetabolism/data/model_build/D/D_pDCs/D_pDCs_model_imat.json"),
        Path("/home/joshl/projects/ImmunoMetabolism/data/model_build/E/E_pDCs/E_pDCs_model_imat.json"),
    ]
    save_path = Path(
        f"/home/joshl/projects/ImmunoMetabolism/results/figures/{models[0].stem.removeprefix('A_').removesuffix('_model_imat')}_heatmap.png"
    )

    fig = build_condition_vs_pathway_heatmap(
        data=models,
        cores=5,
        condition_names=["Age Group A", "Age Group B", "Age Group C", "Age Group D", "Age Group E"],
        save_filepath=save_path,
    )
    fig.show()


if __name__ == "__main__":
    _main()
