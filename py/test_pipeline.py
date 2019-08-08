#!/usr/bin/env python
import os
import importlib
#importlib.import_module('step1_read_raw_data')
from GSEpipeline import *

currentdir = os.getcwd()
dirlist = currentdir.split('/')
projectdir = '/'.join(dirlist[0:-1])

gse2770 = GSEproject('GSE2770',projectdir)

platforms = ['GPL96','GPL97','GPL8300']
maps_dict = download_gsm_id_maps(gse2770.datadir,gse2770.gse,platforms)


# organize_gse_raw_data(gsename)
gse2770.organize_gse_raw_data()
df_outer = gse2770.create_entrez_table_default()
df_clean = gse2770.create_entrez_table_pipeline()