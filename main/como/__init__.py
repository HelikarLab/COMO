import pandas as pd

from como.project import Config
from como.utils import stringlist_to_list

__all__ = ["stringlist_to_list", "Config"]
__version__ = "1.11.0"


def return_placeholder_data() -> pd.DataFrame:
    return pd.DataFrame(data=0, index=pd.Index(data=[0], name="entrez_gene_id"), columns=["expressed", "top"])
