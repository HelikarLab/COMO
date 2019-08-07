#!/usr/bin/env python
import importlib
importlib.import_module('step1_read_raw_data')
#from step1_read_raw_data import *

platforms = ['GPL96','GPL97','GPL8300']
maps_dict = download_gsm_id_maps(datadir,gse,platforms)
