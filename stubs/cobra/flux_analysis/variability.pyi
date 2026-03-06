import pandas as pd
from ..core import Configuration as Configuration, get_solution as get_solution
from ..util import ProcessPool as ProcessPool
from .deletion import single_gene_deletion as single_gene_deletion, single_reaction_deletion as single_reaction_deletion
from .helpers import normalize_cutoff as normalize_cutoff
from .loopless import loopless_fva_iter as loopless_fva_iter
from .parsimonious import add_pfba as add_pfba
from _typeshed import Incomplete
from cobra import Gene as Gene, Model as Model, Reaction as Reaction

logger: Incomplete
configuration: Incomplete

def flux_variability_analysis(model: Model, reaction_list: list[Reaction | str] | None = None, loopless: bool = False, fraction_of_optimum: float = 1.0, pfba_factor: float | None = None, processes: int | None = None) -> pd.DataFrame: ...
def find_blocked_reactions(model: Model, reaction_list: list[Reaction | str] | None = None, zero_cutoff: float | None = None, open_exchanges: bool = False, processes: int | None = None) -> list['Reaction']: ...
def find_essential_genes(model: Model, threshold: float | None = None, processes: int | None = None) -> set['Gene']: ...
def find_essential_reactions(model: Model, threshold: float | None = None, processes: int | None = None) -> set['Reaction']: ...
