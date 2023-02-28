---
title: Repository Structure
sidebar: como_sidebar
permalink: como_repository_structure.html
folder: como
summary: The list of files in COMO as of September 19, 2022.
last_updated: September 19, 2022
---

```
work/
├── data
│   ├── boundary_rxns
│   │   └── bcell_boundary_rxns.csv
│   ├── config_sheets
│   │   ├── disease
│   │   │   ├── disease_data_inputs_naiveB.xlsx
│   │   │   └── disease_data_inputs_smB.xlsx
│   │   ├── microarray_data_inputs_template.xlsx
│   │   ├── mrnaseq_data_inputs_paper_rm.xlsx
│   │   ├── proteomeXchange_urls.csv
│   │   ├── proteomics_data_inputs_paper.xlsx
│   │   ├── scrnaseq_data_inputs_template.xlsx
│   │   └── trnaseq_data_inputs_paper_rm.xlsx
│   ├── data_matrices
│   │   ├── dummy
│   │   │   ├── dummy_microarray_data.csv
│   │   │   ├── dummy_proteomics_data.csv
│   │   │   └── dummy_rnaseq_data.csv
│   │   ├── naiveB
│   │   │   ├── disease
│   │   │   │   ├── gene_counts_matrix_arthritis_naiveB.csv
│   │   │   │   ├── gene_counts_matrix_lupus_a_naiveB.csv
│   │   │   │   └── gene_counts_matrix_lupus_b_naiveB.csv
│   │   │   ├── gene_counts_matrix_full_naiveB.csv
│   │   │   ├── gene_counts_matrix_mrna_naiveB.csv
│   │   │   ├── gene_counts_matrix_total_naiveB.csv
│   │   │   └── protein_abundance_naiveB.csv
│   │   └── smB
│   │       ├── disease
│   │       │   └── gene_counts_matrix_lupus_b_smB.csv
│   │       ├── gene_counts_matrix_full_smB.csv
│   │       ├── gene_counts_matrix_total_smB.csv
│   │       └── protein_abundance_smB.csv
│   ├── force_rxns
│   │   └── bcell_force_rxns.csv
│   ├── GeneralModelUpdatedV2.mat
│   ├── human_proteome_UP000005640_database.fasta
│   ├── iMM_mouse.mat
│   ├── local_files
│   ├── COMO_input
│   │   ├── naiveB
│   │   │   ├── fragmentSizes
│   │   │   │   ├── S1
│   │   │   │   │   ├── naiveB_S1R1_fragment_size.txt
│   │   │   │   │   ├── naiveB_S1R2_fragment_size.txt
│   │   │   │   │   ├── naiveB_S1R3_fragment_size.txt
│   │   │   │   │   └── naiveB_S1R4_fragment_size.txt
│   │   │   │   ├── S2
│   │   │   │   │   ├── naiveB_S2R1_fragment_size.txt
│   │   │   │   │   ├── naiveB_S2R2_fragment_size.txt
│   │   │   │   │   ├── naiveB_S2R3_fragment_size.txt
│   │   │   │   │   └── naiveB_S2R4_fragment_size.txt
│   │   │   │   └── S3
│   │   │   │       ├── naiveB_S3R1_fragment_size.txt
│   │   │   │       ├── naiveB_S3R2_fragment_size.txt
│   │   │   │       └── naiveB_S3R3_fragment_size.txt
│   │   │   ├── geneCounts
│   │   │   │   ├── S1
│   │   │   │   │   ├── naiveB_S1R1.tab
│   │   │   │   │   ├── naiveB_S1R2.tab
│   │   │   │   │   ├── naiveB_S1R3.tab
│   │   │   │   │   └── naiveB_S1R4.tab
│   │   │   │   ├── S2
│   │   │   │   │   ├── naiveB_S2R1.tab
│   │   │   │   │   ├── naiveB_S2R2.tab
│   │   │   │   │   ├── naiveB_S2R3.tab
│   │   │   │   │   └── naiveB_S2R4.tab
│   │   │   │   └── S3
│   │   │   │       ├── naiveB_S3R1.tab
│   │   │   │       ├── naiveB_S3R2.tab
│   │   │   │       └── naiveB_S3R3.tab
│   │   │   ├── insertSizeMetrics
│   │   │   │   ├── S1
│   │   │   │   │   ├── naiveB_S1R1_insert_size.txt
│   │   │   │   │   ├── naiveB_S1R2_insert_size.txt
│   │   │   │   │   ├── naiveB_S1R3_insert_size.txt
│   │   │   │   │   └── naiveB_S1R4_insert_size.txt
│   │   │   │   ├── S2
│   │   │   │   │   ├── naiveB_S2R1_insert_size.txt
│   │   │   │   │   ├── naiveB_S2R2_insert_size.txt
│   │   │   │   │   ├── naiveB_S2R3_insert_size.txt
│   │   │   │   │   └── naiveB_S2R4_insert_size.txt
│   │   │   │   └── S3
│   │   │   │       ├── naiveB_S3R1_insert_size.txt
│   │   │   │       ├── naiveB_S3R2_insert_size.txt
│   │   │   │       └── naiveB_S3R3_insert_size.txt
│   │   │   ├── layouts
│   │   │   │   ├── S1
│   │   │   │   │   ├── naiveB_S1R1_layout.txt
│   │   │   │   │   ├── naiveB_S1R2_layout.txt
│   │   │   │   │   ├── naiveB_S1R3_layout.txt
│   │   │   │   │   └── naiveB_S1R4_layout.txt
│   │   │   │   ├── S2
│   │   │   │   │   ├── naiveB_S2R1_layout.txt
│   │   │   │   │   ├── naiveB_S2R2_layout.txt
│   │   │   │   │   ├── naiveB_S2R3_layout.txt
│   │   │   │   │   └── naiveB_S2R4_layout.txt
│   │   │   │   └── S3
│   │   │   │       ├── naiveB_S3R1_layout.txt
│   │   │   │       ├── naiveB_S3R2_layout.txt
│   │   │   │       └── naiveB_S3R3_layout.txt
│   │   │   ├── prepMethods
│   │   │   │   ├── naiveB_S1R1_prep_method.txt
│   │   │   │   ├── naiveB_S1R2_prep_method.txt
│   │   │   │   ├── naiveB_S1R3_prep_method.txt
│   │   │   │   ├── naiveB_S1R4_prep_method.txt
│   │   │   │   ├── naiveB_S2R1_prep_method.txt
│   │   │   │   ├── naiveB_S2R2_prep_method.txt
│   │   │   │   ├── naiveB_S2R3_prep_method.txt
│   │   │   │   ├── naiveB_S2R4_prep_method.txt
│   │   │   │   ├── S1
│   │   │   │   │   ├── naiveB_S1R1_prep_method.txt
│   │   │   │   │   ├── naiveB_S1R2_prep_method.txt
│   │   │   │   │   ├── naiveB_S1R3_prep_method.txt
│   │   │   │   │   └── naiveB_S1R4_prep_method.txt
│   │   │   │   ├── S2
│   │   │   │   │   ├── naiveB_S2R1_prep_method.txt
│   │   │   │   │   ├── naiveB_S2R2_prep_method.txt
│   │   │   │   │   ├── naiveB_S2R3_prep_method.txt
│   │   │   │   │   └── naiveB_S2R4_prep_method.txt
│   │   │   │   └── S3
│   │   │   │       ├── naiveB_S3R1_prep_method.txt
│   │   │   │       ├── naiveB_S3R2_prep_method.txt
│   │   │   │       └── naiveB_S3R3_prep_method.txt
│   │   │   └── strandedness
│   │   │       ├── S1
│   │   │       │   ├── naiveB_S1R1_strandedness.txt
│   │   │       │   ├── naiveB_S1R2_strandedness.txt
│   │   │       │   ├── naiveB_S1R3_strandedness.txt
│   │   │       │   └── naiveB_S1R4_strandedness.txt
│   │   │       ├── S2
│   │   │       │   ├── naiveB_S2R1_strandedness.txt
│   │   │       │   ├── naiveB_S2R2_strandedness.txt
│   │   │       │   ├── naiveB_S2R3_strandedness.txt
│   │   │       │   └── naiveB_S2R4_strandedness.txt
│   │   │       └── S3
│   │   │           ├── naiveB_S3R1_strandedness.txt
│   │   │           ├── naiveB_S3R2_strandedness.txt
│   │   │           └── naiveB_S3R3_strandedness.txt
│   │   └── smB
│   │       ├── fragmentSizes
│   │       │   ├── S1
│   │       │   │   ├── smB_S1R1_fragment_size.txt
│   │       │   │   ├── smB_S1R2_fragment_size.txt
│   │       │   │   ├── smB_S1R3_fragment_size.txt
│   │       │   │   └── smB_S1R4_fragment_size.txt
│   │       │   └── S2
│   │       │       ├── smB_S2R1_fragment_size.txt
│   │       │       ├── smB_S2R2_fragment_size.txt
│   │       │       └── smB_S2R3_fragment_size.txt
│   │       ├── geneCounts
│   │       │   ├── S1
│   │       │   │   ├── smB_S1R1.tab
│   │       │   │   ├── smB_S1R2.tab
│   │       │   │   ├── smB_S1R3.tab
│   │       │   │   └── smB_S1R4.tab
│   │       │   └── S2
│   │       │       ├── smB_S2R1.tab
│   │       │       ├── smB_S2R2.tab
│   │       │       └── smB_S2R3.tab
│   │       ├── insertSizeMetrics
│   │       │   ├── S1
│   │       │   │   ├── smB_S1R1_insert_size.txt
│   │       │   │   ├── smB_S1R2_insert_size.txt
│   │       │   │   ├── smB_S1R3_insert_size.txt
│   │       │   │   └── smB_S1R4_insert_size.txt
│   │       │   └── S2
│   │       │       ├── smB_S2R1_insert_size.txt
│   │       │       ├── smB_S2R2_insert_size.txt
│   │       │       └── smB_S2R3_insert_size.txt
│   │       ├── layouts
│   │       │   ├── S1
│   │       │   │   ├── smB_S1R1_layout.txt
│   │       │   │   ├── smB_S1R2_layout.txt
│   │       │   │   ├── smB_S1R3_layout.txt
│   │       │   │   └── smB_S1R4_layout.txt
│   │       │   └── S2
│   │       │       ├── smB_S2R1_layout.txt
│   │       │       ├── smB_S2R2_layout.txt
│   │       │       └── smB_S2R3_layout.txt
│   │       ├── prepMethods
│   │       │   ├── S1
│   │       │   │   ├── smB_S1R1_prep_method.txt
│   │       │   │   ├── smB_S1R2_prep_method.txt
│   │       │   │   ├── smB_S1R3_prep_method.txt
│   │       │   │   └── smB_S1R4_prep_method.txt
│   │       │   └── S2
│   │       │       ├── smB_S2R1_prep_method.txt
│   │       │       ├── smB_S2R2_prep_method.txt
│   │       │       └── smB_S2R3_prep_method.txt
│   │       └── strandedness
│   │           ├── S1
│   │           │   ├── smB_S1R1_strandedness.txt
│   │           │   ├── smB_S1R2_strandedness.txt
│   │           │   ├── smB_S1R3_strandedness.txt
│   │           │   └── smB_S1R4_strandedness.txt
│   │           └── S2
│   │               ├── smB_S2R1_strandedness.txt
│   │               ├── smB_S2R2_strandedness.txt
│   │               └── smB_S2R3_strandedness.txt
│   ├── proteomics_entrez_map.csv
│   └── Repurposing_Hub_export.txt
└── py
    ├── async_bioservices
    │   ├── database_convert.py
    │   ├── __init__.py
    │   ├── input_database.py
    │   ├── output_database.py
    │   └── taxon_ids.py
    ├── cluster_rnaseq.py
    ├── cluster_sources.py
    ├── configlib.py
    ├── create_context_specific_model.py
    ├── disease_analysis .py
    ├── GSEpipelineFast.py
    ├── GSEpipeline.py
    ├── __init__.py
    ├── instruments.py
    ├── knock_out_simulation.py
    ├── merge_xomics.py
    ├── microarray.db
    ├── microarray_gen.py
    ├── pipeline.ipynb
    ├── pipeline_paper_demo.ipynb
    ├── project.py
    ├── proteomics
    │   ├── Crux.py
    │   ├── FileInformation.py
    │   ├── FTPManager.py
    │   ├── __init__.py
    │   └── proteomics_preprocess.py
    ├── proteomics_gen.py
    ├── rnaseq_gen.py
    ├── rnaseq_preprocess.py
    └── rscripts
        ├── cluster_samples.R
        ├── cluster_sources.R
        ├── combine_distributions.R
        ├── DGE.R
        ├── fitAffy.R
        ├── fitAgilent.R
        ├── generate_counts_matrix.R
        ├── protein_transform.R
        ├── rnaseq.R
        └── transform.Rmd
```
