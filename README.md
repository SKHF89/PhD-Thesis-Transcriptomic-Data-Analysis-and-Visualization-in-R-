This repository contains comprehensive transcriptome analysis comparing 
PC_Insert (Primary cell culture inserts), Orga_Insert (organoids-derived inserts), and Tissue conditions, including:
DESeq2 differential expression analysis identifying key regulated genes
GSEA (Gene set enrichment analysis) performed on three gene sets:
1) Correlation-filtered genes 
2) Differentially expressed genes (DEGs)
3) Hub genes identified through WGCNA; WGCNA results identifying gene modules and hub genes
Visualization outputs:
PCA plots showing sample relationships
Volcano plots highlighting DEGs
Heatmaps of key epithelial markers
The R script performs extensive bioinformatic analyses with results in organized folders:
GSEA_Results, WGCNA_Results, Visualization_Results. 
Each GSEA analysis examines enrichment across Gene Ontology (biological processes, cellular components, molecular functions) 
and KEGG pathways, with network and dotplot visualizations.
