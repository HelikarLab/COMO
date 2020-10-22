# MADRID: a pipeline for MetAbolic Drug Repurposing IDentification

This is the home page for MADRID pipeline.

## How To Run application
- [Install Docker](https://docs.docker.com/install/)
- `docker login`
- `docker pull docker.pkg.github.com/helikarlab/madrid/madrid:r1.22`
- `docker run --cpus=2 -p 4000:8888 docker.pkg.github.com/helikarlab/madrid/madrid:r1.22` Run docker image and assign 2 CPUs to it.
- Open [http://127.0.0.1:4000](http://127.0.0.1:4000) from your Browser, input token shown in command line terminal from previous step
- In your jupyter notebook, open `/pipelines/py/pipeline.ipynb`
- Upload your configuration and data files `/pipelines/data/` according to the instructions in the notebook, update the file names in the jupyter notebook accordingly.
- Run the notebook step by step, or run the step(s) by your needs








## Flow Charts
![1](./doc/IMG_1.jpg)

![2](./doc/IMG_2.jpg)

![3](./doc/IMG_3.jpg)

![4](./doc/IMG_4.jpg)



## Running pipeline: Example

Identify drug targets against Rheumatoid arthritis using GSMN of naïve CD4+ T cell subtype.  

Follow **How To Run application** instructions to run application and open Jupyter notebook.

**Step 1.** In _"/root/pipelines/data/"_ folder, we included GEO accession numbers of transcriptomics data in input file _“transcriptomics_data_inputs.xls”_ and sample names of proteomics data in input file _“proteomics_data_inputs.xls”_. The protein abundance data is provided in input file  _“ProteomicsDataMatrix.xls”_. 

_"transcriptomics_data_inputs.xls"_  includes Microarray samples of naive CD4+ T cells from GSE22886, GSE43005, GSE22045, and GSE24634. The file _"proteomics_data_inputs.xls"_ file contains sample names of Naive CD4+ T cells in _"ProteomicsDataMatrix.xls"_. 

Running Step 1 in _jupyter_ notebook will generate “gene activity” files based on transcriptomics and proteomics data, as described by Puniya et al., 2020. This will save final output in   _"GeneExpression_Naive_Merged.csv"_  and its path in _"step1_results_files.json"

**Step 2.** Our pipeline includes a modified version of the Recon3D model to use as a reference for model contextualization. 
The modified version of Recon3D in _"/root/pipelines/data/"_ folder is available in _"GeneralModel.mat"_ file. 
This step will use _"GeneExpression_Naive_Merged.csv"_ (created in Step 1) with _"GeneralModel.mat"_ and construct a cell-type-specific model of Naive CD4+ cells. This step will save the output file as _"Naive_SpecificModel.json"_ in _"/root/pipeline/data/"_.

Running step 2 with gene activity and Recon3D will generate a model for naive CD4+ T cells. We can use this model in the next steps. However, we advise users to properly investigate, manually curate, and reupload the refined version in  _"/root/pipeline/data"_ to use in Step 4. We provided already curated versions of the Naive CD4+ T cell model uploaded as _"NaiveModel.mat"_ in  _"/root/pipeline/data/"_. 

**Step 3.** We used a dataset (GSE56649) of rheumatoid arthritis to identify differentially expressed genes (disease genes). We defined accession ids of this dataset in input file _"disease_transcriptomics_data_inputs.xlsx"_

This step will generate files _"Disease_UP_GSE56649.txt"_ and _"Disease_DOWN_GSE56649.txt"_ and save their paths in _"step2_results files.json"_ in _"/root/pipeline/data/"_.  Finally, this step will create _"disease_files"_ variable that will include paths of files for up and downregulated genes. 

**Step 4.**  This step will use the model (constructed in Step 2/ uploaded curated version) and perform knock-out simulations of genes overlapping with the drug-target data file obtained from the ConnectivityMap database. We refined the drug target-data file and included in _"/root/pipeline/data/"_ as _"RepurposingHub.txt"_.

This step will use model (either _"Naive_SpecificModel.json",  or an already curated version uploaded as _"NaiveModel.mat"_), _"Disease_UP_GSE56649.txt"_ and _"Disease_DOWN_GSE56649.txt"_, and  _"RepurposingHub.txt"_ . 

The final output files will include drug targets ranked based on Perturbation Effect Score (PES), as described by Puniya et al., 2020.
Final output folder: _"/root/pipelines/output/"_
The output file  _"d_score.csv"_ will contain Entrez ids of ranked genes and corresponding PES.  The file _"drug_score.csv"_ will contain PES ranked drug targets (Entrez ids and gene symbols) with mapped repurposed drugs. 

## Running pipeline: User-defined data
Create config and data files and upload them to _"/root/pipelines/data/"_  of the _Jupyter_. 

 
1. See _“transcriptomics_data_inputs.xls”_ and create a config file for transcriptomics data. Change the names in Step 1 of the Jupyter notebook. For multiple cell-types, use a separate tab for each cell type data in this file. Example: For Th1 and Th2 CD4+ T cell subtypes, the config file should contain two sheets. Each sheet should include samples related to one cell type. 
2. If proteomics data is available, upload a data file with protein abundances (see _“ProteomicsDataMatrix.xls”_) and define sample names in the config file for proteomics data (see _“proteomics_data_inputs.xls”_).  For multiple cell-types, use the same names for the same subtype in transcriptomics (sheet names) and proteomics (headers) config files. 
        
        #If proteomics data is not available, please use the following file names in Step 1 of the Jupyter notebook.
        proteomics_data_file = 'dummy_proteomics_data.xlsx'
        proteomics_config_file = 'dummy_proteomics_config.xlsx' 
3. We provide a generic human metabolic model in _"GeneralModel.mat"_ file in the pipeline. If user provide a new version of the human metabolic model, it should contain Entrez gene IDs in the gene rules.
4. See _"disease_transcriptomics_data_inputs.xlsx"_ and create disease_gene_file. If using a different file name, please update the name in Step 3 of the Jupyter notebook. 
5. For using the automatically created model(s), no change is required in  Step 4. If the user wants to use refined models, upload them in .mat format and uncomment the relevant lines and change file names in substep 3 of Step 4. 

See final outputs in: _"/root/pipelines/output/"_  



## Resources
* https://opencobra.github.io/
* cobrapy [doc](https://cobrapy.readthedocs.io/en/stable/), [installation](https://github.com/opencobra/cobrapy/blob/master/INSTALL.rst)
* Limma R package [User Guide](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf), [installation](https://bioconductor.org/packages/release/bioc/html/limma.html)
* [affy—analysis of Affymetrix GeneChip data at the probe level](./papers/btg405.pdf)
* [R project](https://www.r-project.org/), [Installation of R, R Studio, R Commander](https://www.andrewheiss.com/blog/2012/04/17/install-r-rstudio-r-commander-windows-osx/)
* [rpy2](https://rpy2.readthedocs.io)
* [biopython](https://biopython.org/wiki/Packages)
* Python: [GEOparse](https://geoparse.readthedocs.io/) [PDF](./doc/geoparse.pdf), R: [GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)
* [R Tutorial](https://www.cyclismo.org/tutorial/R/index.html)
* Python RESTful [Tutorial 1](https://kite.com/blog/python/flask-sqlalchemy-tutorial/), [Tutorial 2: API](https://kite.com/blog/python/flask-restful-api-tutorial/), [Tutorial Index](https://realpython.com/tutorials/api/), [Book: RESTful Web Services](http://restfulwebapis.org/RESTful_Web_Services/), [CURL in Python](https://curl.trillworks.com/#python)
* Docker: [Common components](https://www.digitalocean.com/community/tutorials/the-docker-ecosystem-an-introduction-to-common-components), [Containerization](https://www.digitalocean.com/community/tutorials/the-docker-ecosystem-an-overview-of-containerization), [Official Tutorials](https://docs.docker.com/get-started/), [Gitlab Docker Image](https://docs.gitlab.com/ee/ci/docker/using_docker_build.html), [Gitlab Container Registration](https://docs.gitlab.com/ee/user/packages/container_registry/)
* Google Cloud Platform, Bigquery [Public Datasets](https://www.reddit.com/r/bigquery/wiki/datasets)
* Biomart: ENTREZ ID search, [mygene](https://mygene.info/), [BioDBNet](https://biodbnet-abcc.ncifcrf.gov/db/db2db.php)
* limma [R package](http://bioconductor.org/packages/release/bioc/html/limma.html), [AgiMicroRNA](https://bioconductor.org/packages/release/bioc/html/AgiMicroRna.html), [Analyzing Agilent MicroArray Data in R](https://support.bioconductor.org/p/96655/)
* Bioconductor Tutorials: [Courses and Materials](http://master.bioconductor.org/help/course-materials/),
* Python Package: [Bioservices](https://bioservices.readthedocs.io/en/master/)
* SQLite: [database browser](https://sqlitebrowser.org/dl/)
* Matlab Runtime, [Create Python Packages from Matlab code](https://www.mathworks.com/help/compiler_sdk/gs/create-a-python-application-with-matlab-code.html)
* [Drug Repurposing Hub](https://clue.io/repurposing-app)






## Repository Structure

```bash
├── README.md # README
├── code # R & Matlab code
│   ├── Affy_Script.r
│   ├── Agilent\ Analysis\ and\ cutoffs.r
│   ├── Duplicated_cmapFiles_Sript.R
│   ├── GEO_ID_maps.r
│   ├── Knock_out_simulation_1.m
│   ├── SCORE.R
│   ├── SCORE_DOWN.R
│   ├── SCRIPT_ExtractingFRatio_DE_genes_2.R
│   ├── Script_CMap_FCMat_data.R
│   ├── Script_DT.R
│   ├── Script_EntrezWise_DrugsConMap.R
│   ├── Script_Replicating_multiEntrezIDs.r
│   ├── Script_merging_data.txt
│   ├── Validation_withCMap_STAT.R
│   └── merge_R.R
├── data # Gene Data
│   ├── GSE2770_RAW # Folder of RAW data of GSE2770
│   ├── gpl570entrez.csv
│   ├── gpl8300entrez.csv
│   ├── gpl96entrez.csv
│   └── gpl97entrez.csv
├── doc # documents
│   ├── IMG_1.jpg
│   ├── IMG_2.jpg
│   ├── IMG_3.jpg
│   └── IMG_4.jpg
├── output # Outputs & Intermediate Files
│   └── GSE2770
├── papers # papers
│   ├── Methods.pdf
│   ├── btg405.pdf
│   └── usersguide.pdf
└── py # Python code
    ├── GSEpipeline.py # Pipeline Object based on GSEparse package
    ├── GSEpipelineFast.py # Simplified Pipeline Object
    ├── RESTful_GEO_Fetch.ipynb
    ├── RESTful_Scratch_Pad.ipynb
    ├── Scratch_Pad.ipynb
    ├── Sort_CEL_Files.ipynb
    ├── Test_GSEpipeline.ipynb
    ├── templates # html templates for RESTful
    ├── app.py # app for RESTful
    ├── database_setup.py # database setup for RESTful
    ├── instruments.py # affy and agilent functions
    ├── populate.py # database operation for RESTful
    ├── test_pipeline.py # test script of pipeline
    ├── transcriptomic_gen.ipynb
    ├── transcriptomic_gen.py # Step 1 Test Script Entry
    └── transcriptomics.db # Temporary SQLite database
```
