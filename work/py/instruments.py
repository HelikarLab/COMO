#!/usr/bin/python3
import os
import time
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from bioservices import BioDBNet

class RObject:
    def __init__(self):
        """
        This class is used to access the R-objects located at rscripts/fitAffy.R and rscripts/fitAligent.R
        This is done beacuse lines under the "if main. . ." were originally placed outside of this. This is not ideal
            https://www.freecodecamp.org/news/if-name-main-python-example/
        """
        self._home_dir = os.path.expanduser("~")
        pandas2ri.activate()

        # Import requried R libraries
        self._affy_library = importr("affy")
        self._limma_library = importr("limma")

        self._affy_R_file = open(f"{self._home_dir}/work/py/rscripts/fitAffy.R", "r").read()
        self._agilent_R_file = open(f"{self._home_dir}/work/py/rscripts/fitAgilent.R", "r").read()


    @property
    def affyio(self):
        return SignatureTranslatedAnonymousPackage(self._affy_R_file, "affyio")

    @property
    def agilent(self):
        return SignatureTranslatedAnonymousPackage(self._agilent_R_file, "agilentio")


class AffyIO:
    def __init__(self):
        """
        This class is created so GSEpipeline is able to grab the "affyio" rpy2 object
        This allows us to reuse this section of the code, while allowing us to keep code relevant to this file in a "main" if-statement
        """
        # Get the current user
        home_dir = os.path.expanduser("~")
        # read and translate R functions for handling affy
        # enable R to py conversions
        pandas2ri.activate()

        # import R libraries
        self.affy_library = importr("affy")
        self.limma_library = importr("limma")

        self.affy_R_file = open(f"{home_dir}/work/py/rscripts/fitAffy.R", "r").read()

    @property
    def affyio(self):
        affyio_object = SignatureTranslatedAnonymousPackage(self.affy_R_file, "affyio")
        return affyio_object


# setup agilent df
def agilent_raw(datadir, gsms):
    files = os.listdir(datadir)
    txts = []
    gzs = []
    keys = []
    if gsms:
        for gsm in gsms:
            for file in files:
                if gsm in file:
                    gzs.append(file)
                    txts.append(file[0:-3])
                    keys.append(gsm)
    else:
        for file in files:
            gzs.append(file)
            txts.append(file[0:-3])
            keys.append(file.split("_")[0])
    cols = dict(zip(txts, keys))

    targets = pd.DataFrame(gzs, columns=["FileName"], index=txts)
    df_agilent = RObject().agilent
    df_agilent = RObject().agilent.readagilent(datadir, targets)  # Error shows because this is an R object, PyCharm doesn't know this
    df_agilent = ro.conversion.rpy2py(df_agilent)
    df_temp = pd.read_csv(os.path.join(datadir, "ydf_temp.csv"), header=0)
    df_agilent["ProbeName"] = df_temp["ProbeName"].to_list()

    return df_agilent.rename(columns=cols)


# read agilent outputs for each gsm
def readagilent(datadir, gsms, scalefactor=1.1, quantile=0.95, ):
    df_raw = agilent_raw(datadir, gsms)
    df = df_raw.drop(columns=["Row", "Col"])
    df_negctl = df[df["ControlType"] == -1]
    df_cutoff = scalefactor * df_negctl.drop(columns=["ControlType"]).quantile(
        quantile, axis=0
    )
    df_bool = df.loc[df["ControlType"] == 0, df_cutoff.index.tolist()].gt(
        df_cutoff, axis=1
    )
    idx_ones = df_bool[df_bool.all(axis=1)].index
    idx_zeros = df_bool[~df_bool.all(axis=1)].index
    df.loc[idx_ones, "Express"] = 1
    df.loc[idx_zeros, "Express"] = 0
    df_results = df.loc[df["ControlType"] == 0, :].copy()
    for gsm in gsms:
        col = "{}.cel.gz.1".format(gsm.lower())
        df_results.loc[:, col] = "A"
        df_results.loc[df_bool.loc[:, gsm], col] = "P"
        col = "{}.cel.gz.2".format(gsm.lower())
        df_results.loc[:, col] = 1.0 - quantile
        col = "{}.cel.gz".format(gsm.lower())
        df_results.rename(columns={gsm: col}, inplace=True)
    # df_results.rename(columns={'ProbeName': 'ID'}, inplace=True)
    # df_results.set_index('ID', inplace=True)
    # df_results.set_index('ProbeName', inplace=True)
    return df_results.drop(["ControlType", "SystematicName"], axis=1)


# convert gene ids to entrez
def fetch_entrez_gene_id(
    input_values,
    input_db="Agilent ID",
    output_db: list[str] = None,
    delay=30,
):
    # Set default values for list of strings
    if output_db is None:
        output_db: list[str] = ["Gene ID", "Ensembl Gene ID"]

    s = BioDBNet()
    # input_db = 'Agilent ID'
    # input_values = df_results.ProbeName.tolist()

    df_maps = pd.DataFrame([], columns=output_db)
    df_maps.index.name = input_db
    i = 0
    # for i in range(0,len(input_values),500):
    while i < len(input_values):
        print("retrieve {}:{}".format(i, min(i + 500, len(input_values))))
        df_test = s.db2db(
            input_db, output_db, input_values[i : min(i + 500, len(input_values))], 9606
        )
        if isinstance(df_test, pd.DataFrame):
            df_maps = pd.concat([df_maps, df_test], sort=False)
        elif df_test == "414":
            print("bioDBnet busy, try again in {} seconds".format(delay))
            time.sleep(delay)
            continue
        i += 500
    return df_maps


if __name__ == '__main__':
    home_dir = os.path.expanduser("~")
    # enable R to py conversions
    pandas2ri.activate()

    # import R libraries
    affy = importr("affy")
    limma = importr("limma")

    # Process the affyio functions
    # This is done using a class because the GSEpipeline also utilizes the R-affyio object
    # This is the best method of keeping this information segregated in a "main" statement, while allowing access to other functions
    affyio = AffyIO().affyio

    # read and translate R functions for handling agilent
    fit_aligent_R = open(f"{home_dir}/work/py/rscripts/fitAgilent.R", "r").read()
    agilentio = SignatureTranslatedAnonymousPackage(fit_aligent_R, "agilentio")
