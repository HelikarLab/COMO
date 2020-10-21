# MADRID: a pipeline for MetAbolic Drug Repurposing IDentification

This is the home page for MADRID pipeline.

## How To Run application
- [Install Docker](https://docs.docker.com/install/)
- `docker login`
- `docker pull docker.pkg.github.com/helikarlab/madrid/madrid:r1.21`
- `docker run --cpus=2 -p 4000:8888 docker.pkg.github.com/helikarlab/madrid/madrid:r1.01` Run docker image and assign 2 CPUs to it.
- Open [http://127.0.0.1:4000](http://127.0.0.1:4000) from your Browser, input token shown in command line terminal from previous step
- In your jupyter notebook, open `/pipelines/py/pipeline.ipynb`
- Upload your configuration and data files `/pipelines/data/` according to the instructions in the notebook, update the file names in the jupyter notebook accordingly.
- Run the notebook step by step, or run the step(s) by your needs








## Flow Charts
![1](./doc/IMG_1.jpg)

![2](./doc/IMG_2.jpg)

![3](./doc/IMG_3.jpg)

![4](./doc/IMG_4.jpg)



## Example

Using pipeline to identify drug targets against Rheumatoid arthritis using GSMNs of four CD4+ T cell subtypes: naïve, Th1, Th2, and Th17 developed by Puniya et al., 2020. 

1. In "pipelines/data" folder, we included GEO accession numbers of transcriptomics data in input file 1 (i.e., _“transcriptomics_data_inputs.xls”_) and sample names of proteomics data for each subtype in input file 2 (i.e., _“proteomics_data_inputs.xls”_). The protein abundance file is provided in input file 3 (i.e., _“ProteomicsDataMatrix.xls”_). Running Step 1 in _jupyter_ notebook will generate “gene activity” files based on transcriptomics and proteomics data, as described by Puniya et al., 2020. 

2. Our pipeline includes a modified version of the Recon3D model to use as a reference for model contextualization. Running step 2 with gene activity and Recon3D will generate models for naive, Th1, Th2, and Th17 CD4+ T cell subtypes. These models can be directly used for drug target identification. These models can be further curated and reuploaded. We used already curated versions of developed models provided by Puniya et al., 2020. 

3. We used a dataset of rheumatoid arthritis that is input for Step 3 to identify differentially expressed genes. The up and downregulated genes were used as inputs with curated models in Step 4. 

4. Using constructed models, we perform knock-out simulations based on genes overlapping with input file “RepurposingHub.txt” obtained from the ConnectivityMap database and included in the pipeline. 


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
