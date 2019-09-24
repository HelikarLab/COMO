#!/usr/bin/env python3
import os
import importlib
#importlib.import_module('step1_read_raw_data')
from GSEpipeline import *

currentdir = os.getcwd()
dirlist = currentdir.split('/')
projectdir = '/'.join(dirlist[0:-1])

gse2770 = GSEproject('GSE2770',projectdir)

platforms = ['GPL96','GPL97','GPL8300','GPL570','GPL4685']
# maps_dict = download_gsm_id_maps(gse2770.datadir,gse2770.gse,platforms)


# organize_gse_raw_data(gsename)
gse2770.organize_gse_raw_data()
df_outer = gse2770.get_entrez_table_default(fromcsv=False)
df_clean = gse2770.get_entrez_table_pipeline(fromcsv=False)

print('End of Code')