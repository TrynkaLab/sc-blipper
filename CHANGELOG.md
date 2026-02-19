# v0.0.4-alpha (tbd)

### Major updates
- Added process cnmf_qc for program level QC metrics as well as additional clustering metrics
- Migrated decoupler to python, easing python 3.8 dependency
- Bumped main env to python 3.10

### Minor changes
- Added environment.yml and enviroment_scvi.yml to ease installation
- Removed requirements.txt
- Added installation.md
- Fell back to fgsea v1.32 from v1.35 to enable conda install
- Added script that produces e-distance plot based on clustering
- Added summary file & plot to k_selection with addtional QC metrics

### Bugfixes
- Resolved https://github.com/TrynkaLab/sc-blipper/issues/2
- Resolved https://github.com/TrynkaLab/sc-blipper/issues/3

### Documentation
- Fixed issue with incorrect install instructions

### Configuration
- Added parameter qc_skip_calc_varexp with default false to skip estimation of variance explained
- Retired params.enrich.omnipath_cache_dir
  
  
---


# v0.0.3-alpha (10-02-2026)

### Major changes
- Added support for scVI
- Added support for filtering biotypes using regex pattern (e.g. "protein_coding|lincRNA")
  - Filtering is applied per batch prior to merging
- Added option to completely override cNMF input files (skips all preprocessing)
- Updated way gene id conversion is configued so its much more intuitve. It now works with input and output namespaces instead of booleans.
  - manifest convert_ids > namespace (either 'ensembl' or 'gene_name')
  - config now has 'output_namespace' (either 'ensembl' or 'gene_name')
  - convert.convert_gene_ids is retired as conversion can be skipped if the namespace is matched
  - enrich.input_namespace controls ID conversion for enrich workflow (has no manifest)
- Updated maximal pathway size for GSEA & ORA from 500 to 2000

### Bugfixes
- Fixed bug when specifying filter list for link_h5ad
- Fixed bug in create_cnmf_summary for cases with zero usage causing NAs in scaling leading to a crash

### Documentation
- Updated parameter documentation
- Moved documentation from readme to wiki

### Configuration
- Moved harmony_vars > preprocess.harmony.harmony_vars
- Moved cnmf.n_variable > preprocess.n_variable
- Added preprocess.batch_correction
- Added rn_biotype_filter
- Added enrich.input_namespace
- Added configuration block 'preprocess'
- Added configuration block 'preprocess.scvi'
- Added configuration block 'cnmf.input'
- Added configuration 'enrich.max_pathway_size'
- Removed convert.convert_gene_ids
- Removed is_ensembl_id
- Changed default cnmf.summarize.topn from 5 to 10

### Code strucutre
- Created new subworkflow cnmf_stage to handle batch correction

---

# v0.0.2-alpha (20-11-2025)

## Major Changes
- Added gnomad v4 geneset
- Added function to create a summary table for cnmf
- Added option to invert enrichment in enrich workflow, to enable input of signed pvalues
- Added automated estimation of n_workers
- Added cytokine & TF lists
- Added script to create GEP profile overview plots
- Updated submit script with option for queue and time limit for submitting the nextflow job
- Updated to only report significant postive enrichments in the top summary table
- Updated the way gmt files are provided. By default all gmt's in the assets folder are used
- Updated CD4 marker list


## Bugfixes
- Fixed bug caused by argparse in nmf summarize
- Fixed import issue in enrich workflow
- Fixed issue with staging of output when no enrichment is run
- Fixed an issue with ORA results with zero overlap creating a new test group as the overlap gene column was empty


---


# v0.0.1-alpha

This was too long ago and I don't remember :) 
