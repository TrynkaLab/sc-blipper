

# Genesets


## msigdb genesets
Enrichment will be tested against .gmt files in this folder. gmt files were downloaded from msigdb.org on 20205-08-11
https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2

Gensesets exist in gene symbol and ensembl id versions.
Ensembl id conversion has been done using ensembl v114 by bin/fetch_ensembl_genes.r and bin/convert_gmt_ids.r


## gel.pid.ibd.v8.32.symbols
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