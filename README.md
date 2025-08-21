# sc-blipper

This is a nextflow pipeline for post-proccessing single cell RNAseq datasets and performing gene set enrichments.

There are currently 3 workflows available in the pipeline, below is an overview of the minimal I/O, it can change depending on which options are provided

1. cnmf > runs consensus non-negative matrix factorization and annotation of the outputs 
   - Input:
     - manifest.tsv pointing to seurat.rds and or .h5ad with raw counts
     - params.config file setting parameters
     - Optional: .gmt files for geneset and summary statistics for gwas enrichment
   - Output:
     - CNMF factors for specified K
     - Formatted & merged .h5ad with counts 
     - Ensembl reference file and ID linkers
     - Optional Enrichments for specified K
       - GSEA, ORA, Magma, DecoupleR
2. enrich > runs various enrichment approaches on a matrix
    - Input:
      - gene x condition matrix (can be tranposed as well using `enrich.transpose=true`)
      - params.config file setting parameters
      - .gmt files for geneset and summary statistics for gwas enrichment
    - Output:
      - Merged enrichment file with FDR correction
      - Individual enrichment files
      - Ensembl reference file and ID linkers
3. convert > merge and convert a mixture of seurat and anndata objects (counts)
    - Input:
     - manifest.tsv pointing to seurat.rds and or .h5ad with raw counts
     - params.config file setting parameters
   - Output:
    

Once installed and setup I reccomend running through the wrapper script `sc-blipper`.
General pipeline usage for sanger farm-22 is as such 

```
Usage: sc-blipper <cnmf|enrich|convert> [-c <file.nf>] [-l] -- [nextflow pipeline args]
<cnmf|enrich|convert>     The workflow to run
-c                        <path/to/config.nf> Nextflow config file for the run
-l                        Run nextflow locally instead of submitting to oversubscribed
--                        The rest is passed to nextlfow and overrides -c

Examples:
sc-blipper enrich -c conf.nf -l
sc-blipper enrich -c conf.nf -- --rn_runname hello_world --enrich.input_matrix matrix.tsv
```

For non sanger setups, the runner is fully configurable, mostly pointing to nextflow installs and setting env variables.
Alternatively, it is also fully possible to run the pipeline directly through nextflow with your own runner scripts.
But I will point to the nextflow documentation for setting this up


# Installing (Sanger farm22)
Add the following to your .bashrc
```
export $PATH="$PATH:/software/teamtrynka/installs/sc-blipper/"
```

Then 
```
source ~/.bashrc
```
And thats it!

# Installing (other configurations)

TBD


# Note on gene ids
The pipeline runs on either gene symbols if `convert.is_ensembl_id=false` or on ensembl id if `convert.is_ensembl_id=true`.
You can convert between these two by manipulating these parameters, as well as the manifest.

For example tp merge one h5ad with ensembl ids with a seurat file with gene symbols and end up with gene symbols:

1. for h5ad set convert_ids in manifest as true
2. for seurat set convert_ids in manifest false
3. set `convert.is_ensembl_id=false` (you want gene symbols)
4. set `convert.convert_gene_names=true`

If you indead want to run with ensembl_ids:

1. for h5ad set convert_ids in manifest as false
2. for seurat set convert_ids in manifest true
3. set `convert.is_ensembl_id=true` (you want gene symbols)
4. set `convert.convert_gene_names=true`

By default, gene-ensembl links are downloaded from biomart, the version is controllable through `rn_ensembl_version` (currently 114 is the higest)
A custom id linking file with two columns, old id, new id can be specified with `convert.id_linker`, but this is mostly untested, it should work if the
target is ensembl id or gene symbol, but any other id will not work. I.e. entrez > ensembl id is ok, but ensembl id > entrez may not currently work for some steps.


## CNMF
At the moment the starting point will be QC'ed seurat or scanpy objects.
One or more can be supplied in the manifest. Currently in phase 1 which aims to:

1. Convert to anndata (only if Seurat)
2. Merge objects (only if manifest has more then 1 row)
3. Run cNMF
4. cNMF annotation
    - Progeny
    - Decoupler
    - Gene set enrichment
    - Heatmap of predifned genes
5. Merge cNMF annotations, create anndata
6. Convert cNMF to seurat (optional)
