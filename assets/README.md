
# Markers


## CD4_markers
These have been manually curated 

## lambert_2018_tfs
List of transcription factors curated by lambert et al. https://pubmed.ncbi.nlm.nih.gov/29425488/

## cui_2023_cytokines
List of 86 cytokines from supplementary table 1: https://www.nature.com/articles/s41586-023-06816-9#Sec34
plus some extra chemokine and GZM additions
GZMB
GZMA
GZMK
CCL1
CCL3
CCL4
CCL4L2
CCL3L1
CXCL8
CCL5


# Genesets


## List of genesets

| File | Description |
|------|-------------|
| c1.all.v2025.1.Hs.symbols.gmt | Chromosomal location gene sets |
| c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt | Updated version of KEGG pathways |
| c2.cp.pid.v2025.1.Hs.symbols.gmt | Primary immunodeficiency genes from MSigDB |
| c2.cp.reactome.v2025.1.Hs.symbols.gmt | Reactome pathway gene sets |
| c5.go.bp.v2025.1.Hs.symbols.gmt | Gene Ontology (GO) Biological Process terms |
| c5.go.cc.v2025.1.Hs.symbols.gmt | Gene Ontology (GO) Cellular Component terms |
| c5.go.mf.v2025.1.Hs.symbols.gmt | Gene Ontology (GO) Molecular Function terms |
| c5.hpo.v2025.1.Hs.symbols.gmt | Human Phenotype Ontology terms (rare genetic diseases) |
| essentiality.depmap.hart.symbols.gmt | Cell-essential genes identified from DepMap CRISPR screens |
| gel.pid.ibd.v8.32.symbols.gmt | Primary immunodeficiency genes curated by Genomics England |
| gnomad.v4.lof.symbols.gmt | pLI (Loss of Function Intolerant) genes from GnomAD v4 |
| h.all.v2025.1.Hs.symbols.gmt | Hallmark pathway gene sets from MSigDB |
| otar.drug.gene.20250205.symbols.gmt | Drug target genes curated by Open Targets (release 20250205) |
| sccellfie.metabolic.pathways.symbols.gmt | Human metabolic pathway genes from scCellFie |


## msigdb genesets
Enrichment will be tested against .gmt files in this folder. gmt files were downloaded from msigdb.org on 20205-08-11
https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2

Gensesets exist in gene symbol and ensembl id versions.
Ensembl id conversion has been done using ensembl v114 by bin/fetch_ensembl_genes.r and bin/convert_gmt_ids.r


## gel.pid.ibd.v8.32
This has v8.32 from the genomics england panel for primariy immunodeficiencies and monogenic IBD panel and was manually generated
https://panelapp.genomicsengland.co.uk/panels/398/


## otar.drug.gene.20250205
Drug-gene and phase information downloaded from OpenTargets https://platform.opentargets.org/downloads
Has been filtered into different drug phases as well as into cancer and immune mediated disease drugs.


## sccellfie.metabolic.pathways
Metabolic system annotations obtained from https://github.com/earmingol/scCellFie/tree/main/task_data/homo_sapiens


## essentiality.depmap.hart
Essential genes from depmap and Hart et al 2015

Union of:
- CRISPRInferredCommonEssentials.csv https://depmap.org/portal/data_page/?tab=currentRelease
- AchillesCommonEssentialControls.csv https://depmap.org/portal/data_page/?tab=currentRelease
- Hart et al 2015 https://pubmed.ncbi.nlm.nih.gov/26627737/


## gnomad.v4.lof
This was created using the script 'create_pli_genes.r' which loads the gnomad v4.1 file.
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
 
It contains gene sets of PlI score filtered at 0.8, 0.9, 0.95 and 0.99 as well as the top 50, 100, 250 and 500 genes based 
on the LOF intolerance z-score. The gensets are restricted to cannonical transcripts, and those with an ensembl ID. 