# sc-blipper

sc-blipper is a Nextflow pipeline for post-QC analysis of (single cell) RNAseq data optimized for HPC enviroments. The core features include running consensus non-negative matrix factorization (cnmf), preprocessing (gene id/gene name/format conversion, merging), batch correction (harmony (cnmf adaptation), scvi) and gene set enrichment (fgsea, ora) and enrichment of GWAS signals through MAGMA. For enrichment, it comes bundled with common genesets for humans. The starting point is QCed raw counts data in a .h5ad or .rds seurat object.

The pipeline is highly configurable, but defaults are designed to work out of the box for datasets of 50-100k cells. Larger objects completely possible, but you will likely need to tweak the resource 'label' arguments for each process to allow for more memory. While the pipeline is easy to run with no real Nextflow experiance, instalattion may require some tweaking depending on your HPC configuration. There is a guide for instalation on this wiki. Currently we do not offer singularity containers, but this is on the todo list which should help greatly with portability.

For a detailled usage, install instructions and docs please see the [wiki](https://github.com/TrynkaLab/sc-blipper/wiki)

For the latest version see the dev branch (might be less stable). Versioning is done through tags. For a list of changes see CHANGELOG.md

# Quick usage

Once fully installed and configured, the pipeline can be run by creating a manifest a config file and submitting it through the sc-blipper runner. All configuration can be done through the config file. For a full list of options, please see nextflow.config (parameter docs in the works).

The manifest should have the columns id, file,	namespace and each row represents a .h5ad or .rds dataset to be merged and analyzed. Seurat and h5ad's can be mixed. They must have raw counts data available, any normalizations are ignored, but cell meta on .obs and @meta are carried forward. 
- id: is the name of the batch to run
- file: is the full path to the input data, must end in either .h5ad (counts in .X) or .rds (counts on RNA@counts)
- namespace: is the type of gene ids the file has. Either 'ensembl' or 'gene_name' 

The manifest is provided as the parameter 'rn_manifest'.

Several example config files are provided in conf/examples

```
Usage: /software/teamtrynka/installs/sc-blipper/sc-blipper <cnmf|enrich|convert> [-c <file.nf>] [-lqtw] [-w workdir] -- [nextflow pipeline args]
<cnmf|enrich|convert>           The workflow to run
-c                              <path/to/config.nf> Nextflow config file for the run
-l                              Run nextflow locally instead of submitting to oversubscribed
-w                              Set the nextflow work directory (default: ../workdir)
-q                              Set the queue for the nextflow job (default: oversubscribed)
-w                              Set the time limit for the nextflow job (hours) (default: 120)
-- The rest is passed to nextlfow and overrides what you provide with -c

Examples:
sc-blipper enrich -c conf.nf -l
sc-blipper enrich -c conf.nf -w /path/to/workdir -- --rn_runname hello_world --enrich.input_matrix matrix.tsv
```

# On the whishlist to implement
- starCAT > inferring nmf usages based on a reference run
- scCellFie > metabolite pathway activitiy potential inference 
- singularity container for the software (currently only supports conda)

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
- scVI: https://github.com/scverse/scvi-tools
- scVI: https://www.nature.com/articles/s41592-018-0229-2

