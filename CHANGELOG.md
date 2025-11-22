


# 22-11-2025

### Major changes
- Added support for scVI
- Added support for filtering biotypes using regex pattern (e.g. "protein_coding|lincRNA")
  - Filtering is applied per batch prior to merging
- Added option to completely override cNMF input files
- Updated way gene id conversion is configued so its much more intuitve. It now works with input and output namespaces instead of booleans.
  - manifest convert_ids > namespace (either 'ensembl' or 'gene_name')
  - config now has 'output_namespace' (either 'ensembl' or 'gene_name')
  - convert.convert_gene_ids is retired as conversion can be skipped if the namespace is matched
  - enrich.input_namespace controls ID conversion for enrich workflow (has no manifest)

### Bugfixes
- Fixed bug when specifying filter list for link_h5ad

### Documentation
- Updated parameter documentation

### Configuration
- Moved harmony_vars > preprocess.harmony.harmony_vars
- Moved cnmf.n_variable > preprocess.n_variable
- Added preprocess.batch_correction
- Added rn_biotype_filter
- Added enrich.input_namespace
- Added configuration block 'preprocess'
- Added configuration block 'preprocess.scvi'
- Added configuration block 'cnmf.input'
- Removed convert.convert_gene_ids
- Removed is_ensembl_id

### Code strucutre
- Created new subworkflow cnmf_stage to handle batch correction