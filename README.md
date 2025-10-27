# sc-blipper

This is a nextflow pipeline for post-proccessing single cell RNAseq datasets and performing gene set enrichments.

On the whishlist to implement:
- starCAT > inferring nmf usages based on a reference run
- scCellFie > metabolite pathway activitiy potential inference 
- singularity container for the software (currently only supports conda)

# General usage

If you are not familiar with nextflow, While not needed, I would reccomend reading up a little first as it will make more sense, especially when it comes to resolving errors. I tried to write this so you don't need to get too deep into the weeds tough, nextflow can be a bit of a rabbithole. The pipeline loosely follows nf-core principles.

Some resources
- https://training.nextflow.io/2.1/#nextflow-for-science
- https://nf-co.re/docs/usage/getting_started/terminology
- https://nf-co.re/docs/usage/troubleshooting/basics

There are currently 3 workflows available in the pipeline, below is an overview of the minimal I/O, it can change depending on which options are provided

1. cnmf > runs consensus non-negative matrix factorization and annotation of the outputs 
2. enrich > runs various enrichment approaches on a matrix
3. convert > merge and convert a mixture of seurat and anndata objects (counts)    

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
# Set default LSF group, important for NF to work, change for your group here
export LSB_DEFAULTGROUP=teamtrynka

export PATH="$PATH:/software/teamtrynka/installs/sc-blipper/"
```

Then 
```
source ~/.bashrc
```
And thats it, you are good to go!

# Installing (other configurations)

1. Make sure you have nextflow available (>=25.04.6)
2. Clone the repo
3. Create a conda env following the instructions in requirements.txt (needs some manual pacakges, future will add singularity containers)
4. Update `nextflow.config` or override with your config the path to the conda env (`params.rn_conda=/path/to/env"`)
5. Add a profile to work with your cluster configuration (can be put in './conf' folder). Also check if any enviroment variables need to be set for your scheduler.
6. Add the new profile to the `nextflow.config` `profiles{}` block
7. (optional) Update the runner script `sc-blipper` as the primary entry point (By default works with LSF, easy to update to SLURM)
8. (optional) Add the runner script `sc-blipper` to your path


# Configuring a pipeline run

A run can be configured using a nextflow config file supplied to `sc-blipper -c <config.config>`. Please see `conf/example_<x>.config` for some examples. Also see `nextflow.config` for a full list of available options.

# Note on gmt files

Reference pathways are provided as .gmt files in the assets folder. The ones you want to use can bet set with `enrich.gmt_files`. By default none are set. They are bundled in two flavours, as ensembl id or as gene symbols. 

# Note on gene ids
The pipeline runs on either gene symbols if `convert.is_ensembl_id=false` or on ensembl id if `convert.is_ensembl_id=true`.
You can convert between these two by manipulating these parameters, as well as the manifest.
Ensembl > gene name mapping is not unique. To make sure gene names are unique and all genes are preserved, the following strategy is applied:

1. Sort ensembl on chromosome 1-22,X,Y, MT
2. Sort so biotype "protein_coding" apears first
3. Set missing gene names to NO_NAME
4. Make genes unique by appending a number at the end if duplicated `<gene_name>_<number>`

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
target is ensembl id or gene symbol, but others are untested and might fail for some of the steps.


# Note on configuring resource limits

The default resource labels have been tested and optimized to run with an object of ~50k cells. If a job crashes, Nextflow will attempt it again, doubling resource requirements where possible up to 2 times. However this can be quite wastefull, so with larger objects of million+ cells you may need to change resource labels, particularly for the `convert` and `cnmf` proccesses. The file `conf/processes.config` has a list of available resource labels. The resource labels for a config can be adjused by the `label` parameter. For instance to give the cnmf proccesses more memory set `cnmf.label='medium'`. Should a suitable resource label not be available, you can define your own in your config file, see the `conf/processes.config` for examples.


# Note on finalizing output and cleaning up after a successfull run

The pipeline results uses Nextflow publishDir directive, this means the output in results is linked to the process output in the workdir. This is nice and efficient for organizing output, especially when you are messing with the pipeline settings, as it avoids duplicating things. However, its not so handy for archiving or finalizing results. For this reason NEVER REMOVE THE WORKDIR BEFORE FINALIZING. To finalize the results and make them ready for backup etc:

```
rsync -rP --copy-links results results_final

```

This will make a deep copy of the results, copying all the symlinked files as files, not links. Then its safe to remove the results and workdir folders. This makes it impossible to pick a run up halfway, and you will need to start fresh after doing this.

```
rm -r results
rm -r workdir
```

-------------------------------------------------------

# Wokrflow: CNMF

## What it does
This workflow takes raw counts from one or multiple files, optionally harmony corrects them, runs cNMF and then performs various enrichment analysis on the gene spectra scores from the cNMF to aid in annotation. 

1. Convert to anndata (only if Seurat)
2. Merge objects (only if manifest has more then 1 row)
3. Run cNMF
4. cNMF annotation
    - Progeny
    - Decoupler
    - Gene set enrichment (GSEA / ORA)
    - Magma (genomic enrichment of GWAS results)
5. Merge cNMF annotations, create h5ad

## General IO
The input is one or multiple seurat and or h5ad files with counts / .X set. The pipeline works with .h5ad so seurat objects are first converted

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


## Output files

The output for the pipeline is organized into a couple of folders
```
├── cnmf
│   └── consensus > output of the cNMF consensus process, the intermediate cnmf outputs are not stages and just live in workdir for now (they are not very usefull)
├── h5ad
│   ├── cnmf > h5ad with HVGs with (optionally) harmony corrected counts
│   ├── merged > Merged h5ad optionally gene name converted and subsetted to genes of interest
│   └── per_batch > Input converted to h5ad (if input is h5ad, this is just linked)
├── magma 
│   └── per_trait > per trais association results from magma
└── reference
    ├── ensembl > ensembl reference files
    └── magma > the magma reference files (pre-calculated gene scores for each trait, can be re-used between pipeline runs)
```

#### cnmf
```
└── consensus
    └── <rn_runname>
        ├── enrichments > Merged enrichment results for all K per database
        ├── final > Merged enrichment resuls for all K and all databases
        ├── k_<X>
        │   ├── enrichments > Enrichment results per database
        │   └── final > Enrichment results for all databases merged
        └── k_selection > files to aid in selecting optimal K value
```

So for k=18 the output might look like this. The output should be fairly self explanatory. For a detailled breakdown on the nuances of the cNMF output please refer to:

- https://github.com/dylkot/cNMF
- https://github.com/dylkot/cNMF/issues/59#issuecomment-2099453595

```
.
├── enrichments
│   ├── k_18_collectri_activities.tsv.gz
│   ├── k_18_collectri_activity_matrix.tsv.gz
│   ├── k_18_collectri.enrich.gz
│   ├── k_18_fgsea_results.enrich.gz
│   ├── k_18_fgsea_results_top5.pdf
│   ├── k_18_fgsea_results_top5.tsv
│   ├── k_18_fgsea_results.tsv.gz
│   ├── k_18_merged_magma_results.enrich.gz
│   ├── k_18_merged_magma_results.tsv.gz
│   ├── k_18_ora_results.enrich.gz
│   ├── k_18_ora_results_top5.pdf
│   ├── k_18_ora_results_top5.tsv
│   ├── k_18_ora_results.tsv.gz
│   ├── k_18_pathways.pdf
│   ├── k_18_pathways_scaled.pdf
│   ├── k_18_progeny_activities.tsv.gz
│   ├── k_18_progeny_activity_matrix.tsv.gz
│   ├── k_18_progeny.enrich.gz
│   ├── k_18_signif_tfs.pdf
│   └── k_18_signif_tfs_scaled.pdf
├── final
│   ├── k_18_merged_nominal.tsv
│   └── k_18_merged.tsv.gz
├── immerse_pilot2_all_modalities_lrpprc_only.clustering.k_18.dt_0_1.png
├── immerse_pilot2_all_modalities_lrpprc_only.gene_spectra_score.k_18.dt_0_1.txt.gz
├── immerse_pilot2_all_modalities_lrpprc_only.gene_spectra_tpm.k_18.dt_0_1.txt.gz
├── immerse_pilot2_all_modalities_lrpprc_only.k_18.dt_0_1.h5ad
├── immerse_pilot2_all_modalities_lrpprc_only.spectra.k_18.dt_0_1.consensus.txt.gz
├── immerse_pilot2_all_modalities_lrpprc_only.starcat_spectra.k_18.dt_0_1.txt.gz
└── immerse_pilot2_all_modalities_lrpprc_only.usages.k_18.dt_0_1.consensus.txt.gz
```

#### H5ad

Optionally an h5ad file is created which has the usages as `X` and the spectra score as `.varm['spectra']`. The gene names for the spectra score can be found in `.uns['spectra_names']`. The `.obsm` contains the merged annotations from the original input for ease of use.



# Workflow: Enrich


## What it does

This workflow takes a numeric matrix (lfc/beta's/pvalues/-1,0,1) and runs enrichment analysis on it. It has strong overlaps with cNMF, but is slightly different in that the starting point is generalized to work with any numeric matrix. It is very configurable, so lots of options.


## General IO
- Input:
  - gene x condition matrix (can be tranposed if it is condition x gene, using `enrich.transpose=true`)
  - params.config file setting parameters
  - .gmt files for geneset and summary statistics for gwas enrichment
- Output:
  - Merged enrichment file with FDR correction
  - Individual enrichment files
  - Ensembl reference file and ID linkers

## Settings

### Required input
Input matrix is provided with `enrich.input_matrix=<path/to/file>`

Reference databases are provided with `enrich.gmt_files="/path/to/file.gmt,/path/to/file2.gmt"`
The gmt files must be in the target namespace. Reference databases in ensembl ids and HCNG symbols are provided in the assets folder.

If not running magma (`enrich.run_magma=false`) this is all you need to set
If you are running magma, you also need to either provide:

- `magma.manifest_sumstats` A manifest tsv file with trait name, snpcol, pvalcol, path to sumarry stats
-  `magma.ld_reference` A suitable LD reference panel for inferring the gene-gene correlations

If using pre-computed magma scores:
- `magma.manifest_magma` A path to a previous pipeline run magma manifest


### Gene universe (optional)
By default, the universe is all genes included in the matrix, this works well for DE tests for instance, as the background forms all tested genes. However there might be cases where the tested geneset is pre-enriched for something, or you tested all genes, so LFCs might have an expression bias. 

In these cases you can set `enrich.universe=<path/to/file>` to a file that has one column with the gene names/gene ids in the output namespace. So if `convert.convert_gene_names=true` make sure you provide this in the target namespace, not the input namespace.

NOTE: GSEA runs on numeric values, so if your universe is much larger then your tested genes, its best not to run this test by setting `enrich.run_gsea=false`


### Enrichment tests to run (optional)
There are 4 enrichment tests implemented, these can be toggeled on and off with the following flags
```
enrich {
    run_gsea = true
    run_ora = true
    run_decoupler = true
    run_magma = true
}
```

### Annotate the results table (optional)
You may wish to annotate the final results table with additional metadata. In this case these can be provided with `enrich.annotate` which should be the path to a tsv file containing the column annotations to add. Are assumed to be in the same order as colnames (I think). 


# Workflow: Convert


## What it does
Takes one or more seurat or h5ad files, optionally coverts gene names or ids and optionally run harmony batch correction.


## General IO
- Input:
 - manifest.tsv pointing to seurat.rds and or .h5ad with raw counts
 - params.config file setting parameters
- Output:
  - Merged, id converted h5ad
  - Optionally harmony corrected counts using the same process as described for cNMF



# References / acknowledgements

This is a wrapper pipline that depends on previously developped tools cNMF, GSEA, decoupleR and MAGMA. Please cite the original publications if you use this:
- cNMF: https://elifesciences.org/articles/43803
- cNMF: https://github.com/dylkot/cNMF
- GSEA: https://www.pnas.org/doi/10.1073/pnas.0506580102
- fgsea: https://github.com/alserglab/fgsea
- fgsea: https://www.biorxiv.org/content/10.1101/060012v3
- decoupleR: https://saezlab.github.io/decoupleR/
- decoupleR: https://doi.org/10.1093/bioadv/vbac016
- MAGMA: https://doi.org/10.1371/journal.pcbi.1004219
- MAGMA: https://cncr.nl/research/magma/

