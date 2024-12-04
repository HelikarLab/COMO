import pandas as pd


def gene_info_migrations(df: pd.DataFrame) -> pd.DataFrame:
    """Migrate gene info DataFrame to the latest version."""
    return df.rename(columns={"hgnc_symbol": "gene_symbol"}) if "hgnc_symbol" in df.columns else df
