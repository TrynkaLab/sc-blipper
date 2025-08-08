# sc-blipper

This is a nextflow pipeline for analyzing single cell RNAseq datasets. At the moment the starting point will be QC'ed seurat or scanpy objects.
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
