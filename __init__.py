"""BioServices
import bioservices
u = bioservices.uniprot.UniProt()
u.search("ZAP70")
see online documentation for details on pypi or github.
"""
import sys

import pkg_resources

from easydev import CustomConfig

__version__ = "1.6.0"
try:
    version = pkg_resources.require("bioservices")[0].version
    __version__ = version
except:
    version = __version__

import colorlog

logger = colorlog.getLogger("bioservices")


# Initialise the config directory if not already done
configuration = CustomConfig("bioservices", verbose=False)
bspath = configuration.user_config_dir

# Add bioservices.uniprot to sys.modules to prevent cycles in our imports
# import bioservices.uniprot
# bioservices.uniprot  # Stop flake8 error

from . import settings
from .settings import *

from . import services
from .services import *

from . import biodbnet
from .biodbnet import *

from . import biocarta
from .biocarta import *

from . import biomodels
from .biomodels import *

from . import cog
from .cog import *

from . import chebi
from .chebi import *

from . import ena
from .ena import ENA

from . import eva
from .eva import *

from . import ensembl
from .ensembl import *

# moved to attic in bioservices v1.6
# from . import geneprof
# from .geneprof import *

from . import kegg
from .kegg import *

from . import hgnc
from .hgnc import *

from . import intact
from .intact import *

from . import pubchem
from .pubchem import *

from . import pfam
from .pfam import *

from . import rhea
from .rhea import *

from . import xmltools
from .xmltools import *

from . import wikipathway
from .wikipathway import *

from . import omnipath
from .omnipath import *

from . import pdb
from .pdb import *

from . import pdbe
from .pdbe import *

from . import pride
from .pride import *

from . import uniprot
from .uniprot import *

from . import unichem
from .unichem import *

from . import rnaseq_ebi
from .rnaseq_ebi import RNASEQ_EBI

from . import reactome
from .reactome import *

from . import quickgo
from .quickgo import *


from . import chembl
from .chembl import *

# from . import picr
# from .picr import *

from . import psicquic
from .psicquic import *

from . import ncbiblast
from .ncbiblast import *

# from . import readseq
# from .readseq import *

from . import biogrid
from .biogrid import *

from . import arrayexpress
from .arrayexpress import *

from . import biomart
from .biomart import *

from . import eutils
from .eutils import *

from . import pathwaycommons
from .pathwaycommons import *

from . import muscle
from .muscle import *

from . import dbfetch
from .dbfetch import *

from . import bigg
from .bigg import BiGG

from . import omicsdi
from .omicsdi import OmicsDI

# sub packages inside bioservices.

# import mapping
from . import apps

# import dev
