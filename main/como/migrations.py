import pandas as pd


def gene_info_migrations(df: pd.DataFrame) -> pd.DataFrame:
    """Migrate gene info DataFrame to the latest version.

    Args:
        df: The input DataFrame containing gene information.

    Returns:
        The migrated DataFrame with updated column names.
    """
    return df.rename(columns={"hgnc_symbol": "gene_symbol"}) if "hgnc_symbol" in df.columns else df
