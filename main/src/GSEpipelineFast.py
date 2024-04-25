import os
import tarfile
import urllib.request

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from fast_bioservices import BioDBNet, Input, Output
from rpy2.robjects import pandas2ri

import instruments
from GSEpipeline import load_gse_soft
from instruments import AffyIO

pandas2ri.activate()


# Input: Extract Gene Info from GEO DataSets

# gse = load_gse_soft(gsename)


def download_gsm_id_maps(datadir, gse, gpls: list[str] = None, vendor="affy"):
    """
    download ID to ENTREZ_GENE_ID maps, create a csv file for each platform, and return dictionary
    :param gpls:
    :param vendor:
    :param gse: gse object
    :param datadir: path to save the csv files
    :return: dictionary of maps indexed by platform
    """
    # Do not allow gpls to be mutable
    # From: https://florimond.dev/en/posts/2018/08/python-mutable-defaults-are-the-source-of-all-evil/
    if gpls is None:
        gpls: list[str] = ["GPL96", "GPL97", "GPL8300"]
    
    biodbnet = BioDBNet()
    
    for gpl in gpls:
        table = gse.gpls[gpl].table.copy()
        if vendor.lower() == "affy":
            temp = table[["ID", "ENTREZ_GENE_ID"]]
        
        elif vendor.lower() == "agilent":
            input_values = table.loc[
                table["CONTROL_TYPE"] == "FALSE", "SPOT_ID"
            ].tolist()
            
            temp = biodbnet.db2db(
                input_values=input_values,
                input_db=Input.AGILENT_ID,
                output_db=[
                    Output.GENE_ID,
                    Output.ENSEMBL_GENE_ID
                ],
            )
            
            temp.drop(columns=["Ensembl Gene ID"], inplace=True)
            temp.reset_index(inplace=True)
            temp.rename(
                columns={
                    Input.AGILENT_ID.value: "ID",
                    Output.GENE_ID.value: "ENTREZ_GENE_ID"
                }, inplace=True
            )
            temp.replace(to_replace="-", value=np.nan, inplace=True)
        
        else:
            print("Unsupported Platform: {}".format(gpl))
            continue
        
        # Save to file
        filefullpath = os.path.join(datadir, "{}entrez.csv".format(gpl.lower()))
        temp.to_csv(filefullpath, index=False)
        # Single Table


class GSEproject:
    def __init__(self, gsename, querytable, rootdir="../"):
        self.gsename = gsename
        # Setup paths
        self.querytable = querytable
        self.rootdir = rootdir
        self.datadir = os.path.join(self.rootdir, "data")
        self.outputdir = os.path.join(self.rootdir, "output")
        self.gene_dir = os.path.join(self.datadir, self.gsename + "_RAW")
        print(
            "Initialize project ({}):\nRoot: {}\nRaw data: {}".format(
                self.gsename, self.rootdir, self.gene_dir
            )
        )
        self.gsm_platform = self.get_gsm_platform()
        pairs = querytable.loc[:, ["GPL ID", "Instrument"]].drop_duplicates()
        gpls = pairs["GPL ID"].tolist()
        vendors = pairs["Instrument"].tolist()
        self.platforms = dict(zip(gpls, vendors))
        self.download_samples()
    
    def organize_gse_raw_data(self):
        """
        Organize raw data at local folder
        :return:
        """
        # create a folder for each platform
        for key in self.platforms.keys():
            platformdir = os.path.join(self.gene_dir, key)
            if not os.path.exists(platformdir):
                os.makedirs(platformdir)
                print("Path created: {}".format(platformdir))
            else:
                print("Path already exist: {}".format(platformdir))
        
        # Move Corresponding Cel files to Folders
        onlyfiles = [f for f in os.listdir(self.gene_dir) if f.endswith(".gz")]
        cnt = 0
        
        for file in onlyfiles:
            filelist = file.split(".")
            prefix = filelist[0]
            if prefix in self.gsm_platform:
                platform = self.gsm_platform[prefix]
                platformdir = os.path.join(self.gene_dir, platform)
                src_path = os.path.join(self.gene_dir, file)
                dst_path = os.path.join(platformdir, file)
                os.rename(src_path, dst_path)
                print("Move {} to {}".format(src_path, dst_path))
                cnt += 1
        print("{} raw data files moved.".format(cnt))
    
    def get_gsm_tables(self):
        """
        get gsm maps in table
        :return:
        """
        gsm_tables = {}
        for gpl, vendor in self.platforms.items():
            filename = "{}entrez.csv".format(gpl.lower())
            filepath = os.path.join(self.datadir, filename)
            if not os.path.isfile(filepath):
                # Could improve to automatic download new tables based on platform
                gse = load_gse_soft(self.gsename)
                download_gsm_id_maps(self.datadir, gse, gpls=[gpl], vendor=vendor)
                print("Skip Unsupported Platform: {}, {}".format(gpl, vendor))
                # continue
            temp = pd.read_csv(filepath)
            df = temp[["ID", "ENTREZ_GENE_ID"]]
            df.set_index("ID", inplace=True)
            gsm_tables[gpl] = df
        
        return gsm_tables
    
    def get_gsm_platform(self):
        """
        Classify the Samples by Platform
        get the platform of each gsm
        :return: dictionary key: gsm, value: platform such as 'GEL96', 'GEL97', 'GEL8300'
        """
        keys = self.querytable["Samples"].str.upper().tolist()
        values = self.querytable["GPL ID"].str.upper().tolist()
        gsm_platform = dict(zip(keys, values))
        
        return gsm_platform
    
    def gsms_included_by(self, df):
        for gsm in self.gsm_platform.keys():
            included = False
            gsml = gsm.lower()
            for key in list(df):
                if gsml == key.lower().split(".")[0]:
                    included = True
                    break
            if not included:
                return False
        
        return True
    
    def get_entrez_table_pipeline(self, fromcsv=True):
        """
        create ENTREZ ID based table from gse
        :return: pandas dataframe for table of GSE
        """
        
        filefullpath = os.path.join(
            self.gene_dir, "{}_sc500_full_table.csv".format(self.gsename)
        )
        if fromcsv and os.path.isfile(filefullpath):
            try:
                df_clean_sc500 = pd.read_csv(filefullpath)
                df_clean_sc500.dropna(axis="columns", how="all", inplace=True)
                df_clean_sc500.dropna(how="all", inplace=True)
                df_clean_sc500 = df_clean_sc500[
                    ~df_clean_sc500.index.duplicated(keep="last")
                ]
                if self.gsms_included_by(df_clean_sc500) and not df_clean_sc500.empty:
                    return df_clean_sc500
                else:
                    print("Need Append GSMs")
            except:
                print("Unable to read {}")
        
        print("Create new table: {}".format(filefullpath))
        gsm_maps = self.get_gsm_tables()
        
        if not any(gsm_maps):
            print("Not available, return empty dataframe")
            return pd.DataFrame([])
        
        # Ready Affy files from folders
        biodbnet = BioDBNet()
        gsm_tables_sc500 = {}
        for key, vendor in self.platforms.items():
            platformdir = os.path.join(self.gene_dir, key)
            print("{} Read Path: {}".format(vendor, platformdir))
            if os.path.exists(platformdir):
                if vendor.lower() == "affy":
                    # Get the AffyIO R-object from instruments.py
                    outputdf = AffyIO().affyio
                    outputdf = outputdf.readaffydir(platformdir)
                    outputdf = ro.conversion.rpy2py(outputdf)
                elif vendor.lower() == "agilent":
                    
                    outputdf = instruments.readagilent(platformdir, list(self.gsm_platform.keys()))
                    
                    gsm_maps[key] = biodbnet.db2db(
                        input_values=list(map(str, list(outputdf["ProbeName"]))),
                        input_db=Input.AGILENT_ID,
                        output_db=[Output.GENE_ID],
                    )
                    
                    gsm_maps[key].rename(columns={"Gene ID": "ENTREZ_GENE_ID"}, inplace=True)
                else:
                    print("Unsupported Platform {} and Vendor {}".format(key, vendor))
                    continue
            else:
                print("Path not exist: {}".format(platformdir))
                continue
            
            drop_idx = np.where(gsm_maps[key]["ENTREZ_GENE_ID"] == "-")[0].tolist()
            outputdf.drop(outputdf.index[drop_idx], axis=0, inplace=True)
            gsm_maps[key].drop(gsm_maps[key].index[drop_idx], axis=0, inplace=True)
            outputdf["ENTREZ_GENE_ID"] = gsm_maps[key]["ENTREZ_GENE_ID"].to_list()
            gsm_tables_sc500[key] = outputdf
        
        # Drop rows without ENTREZ GENE ID, set index to ENTREZ
        for key in self.platforms.keys():
            gsm_tables_sc500[key].dropna(subset=["ENTREZ_GENE_ID"], inplace=True)
            gsm_tables_sc500[key].set_index("ENTREZ_GENE_ID", inplace=True)
            print("gsm table drop: ", gsm_tables_sc500[key])
        
        # Merge tables of platforms
        df_outer_sc500 = None
        for key in self.platforms.keys():
            print("{}: {}".format(key, gsm_tables_sc500[key].shape))
            if df_outer_sc500 is None:
                df_outer_sc500 = gsm_tables_sc500[key]
            else:
                df_outer_sc500 = pd.merge(
                    df_outer_sc500,
                    gsm_tables_sc500[key],
                    on="ENTREZ_GENE_ID",
                    how="outer",
                )
        
        df_outer_sc500.dropna(how="all", inplace=True)
        print("Full: {}".format(df_outer_sc500.shape))
        df_outer_sc500.rename(str.lower, axis="columns", inplace=True)
        keys = []
        vals = []
        gsms_loaded = []
        
        for col in list(df_outer_sc500):
            if ".cel.gz" in col:
                strs = col.split(".cel.gz")
                gsm = strs[0].split("_")[0]
                newcol = "{}.cel.gz{}".format(gsm, strs[-1])
                vals.append(newcol)
                keys.append(col)
                gsms_loaded.append(gsm)
        
        df_outer_sc500.rename(columns=dict(zip(keys, vals)), inplace=True)
        gsms_loaded = list(set(gsms_loaded).union(set(self.gsm_platform.keys())))
        
        # Remove duplicated items, keep largest VALUE for each GSM
        if "df_clean_sc500" not in locals():
            df_clean_sc500 = pd.DataFrame([], index=df_outer_sc500.index)
            df_clean_sc500 = df_clean_sc500[
                ~df_clean_sc500.index.duplicated(keep="first")
            ]
        elif df_clean_sc500.empty:
            df_clean_sc500 = pd.DataFrame([], index=df_outer_sc500.index)
            df_clean_sc500 = df_clean_sc500[
                ~df_clean_sc500.index.duplicated(keep="first")
            ]
        else:
            df_clean_sc500.set_index("ENTREZ_GENE_ID", inplace=True)
            placeholder = pd.DataFrame(
                [], columns=["placeholder"], index=df_outer_sc500.index
            )
            placeholder["placeholder"] = 0
            placeholder.index.name = "ENTREZ_GENE_ID"
            df_clean_sc500 = pd.merge(
                df_clean_sc500, placeholder, on="ENTREZ_GENE_ID", how="outer"
            )
            df_clean_sc500 = df_clean_sc500[
                ~df_clean_sc500.index.duplicated(keep="last")
            ]
        
        for key in gsms_loaded:
            key_low = key.lower()
            col1, col2, col3 = (
                "{}.cel.gz".format(key_low),
                "{}.cel.gz.1".format(key_low),
                "{}.cel.gz.2".format(key_low),
            )
            
            try:
                temp = df_outer_sc500.loc[:, [col1, col2, col3]]
            
            except:
                if key in list(self.gsm_platform.keys()):
                    print("{} not in df_outer_sc500".format(key))
                
                continue
            
            temp.sort_values(by=["ENTREZ_GENE_ID", col1], inplace=True)
            temp = temp[~temp.index.duplicated(keep="last")]
            df_clean_sc500[col1] = temp[col1]
            df_clean_sc500[col2] = temp[col2]
            df_clean_sc500[col3] = temp[col3]
        
        # save to csv file
        try:
            df_clean_sc500.set_index("ENTREZ_GENE_ID", inplace=True)
        
        except:
            pass
        
        df_clean_sc500.dropna(axis="columns", how="all", inplace=True)
        df_clean_sc500.dropna(how="all", inplace=True)
        
        try:
            df_clean_sc500.drop(columns=["placeholder"], inplace=True)
        except:
            pass
        
        df_clean_sc500.sort_index(inplace=True)
        df_clean_sc500.to_csv(filefullpath)
        print("Full table saved to:\n{}".format(filefullpath))
        
        return df_clean_sc500
    
    def download_raw(self, overwrite=False):
        # check if path created
        if (not os.path.isdir(self.gene_dir)) or overwrite:
            os.makedirs(self.gene_dir, exist_ok=True)
            url = (
                "https://www.ncbi.nlm.nih.gov/geo/download/?acc={}&format=file".format(
                    self.gsename
                )
            )
            filefullpath = os.path.join(self.gene_dir, "{}_RAW.tar".format(self.gsename))
            if not os.path.isfile(filefullpath):
                print("Download Raw File: {}".format(filefullpath))
                urllib.request.urlretrieve(url, filefullpath)
            else:
                print("File already exist: {}".format(filefullpath))
            tfile = tarfile.open(filefullpath)
            tfile.extractall(path=self.gene_dir)
            os.remove(filefullpath)
            print("Remove Raw File: {}".format(filefullpath))
            self.organize_gse_raw_data()
        
        else:
            pass
    
    def download_samples(self, overwrite=False):
        os.makedirs(self.gene_dir, exist_ok=True)
        for gsm, gpl in self.gsm_platform.items():
            platformdir = os.path.join(self.gene_dir, gpl)
            os.makedirs(platformdir, exist_ok=True)
            sample_url = (
                "https://www.ncbi.nlm.nih.gov/geo/download/?acc={}&format=file".format(
                    gsm
                )
            )
            filefullpath = os.path.join(self.gene_dir, "{}.tar".format(gsm))
            
            if (not os.path.isfile(filefullpath)) or overwrite:
                urllib.request.urlretrieve(sample_url, filefullpath)
                print("Retrieve Sample: {}".format(filefullpath))
            
            else:
                print("Sample exist: {}".format(filefullpath))
                continue
            
            tfile = tarfile.open(filefullpath)
            tfile.extractall(path=platformdir)
            # os.remove(filefullpath) # keep to avoid re-download
            print("Extract to: {}".format(platformdir))
        print("Retrieve Samples Completed.")
    
    def calculate_z_score(self, df, to_csv=False):
        cols = list(df)
        result = pd.DataFrame([], index=df.index)
        for col in cols:
            if ".gz." not in col:
                score = np.log2(df[col])
                result[col] = (score - score.mean()) / score.std(ddof=0)
        if to_csv:
            filefullpath = os.path.join(
                self.gene_dir, "{}_data_z.csv".format(self.gsename)
            )
            result.to_csv(filefullpath)
        
        return result
