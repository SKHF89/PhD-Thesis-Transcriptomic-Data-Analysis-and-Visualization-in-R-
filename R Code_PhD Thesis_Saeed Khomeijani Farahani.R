############    Start of the Code  #############

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(zFPKM)
library(DESeq2)
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(BiocParallel)
library(enrichplot)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(circlize)

# Load Sheep OrgDb from AnnotationHub
ah <- AnnotationHub()
sheep_orgdb <- ah[["AH114632"]]  # Load sheep OrgDb
print(sheep_orgdb)
cat("Number of genes in sheep_orgdb:", length(keys(sheep_orgdb)), "\n")

# ------------------------#
# Step 1: Load and Filter FPKM Data with zFPKM
# ------------------------#

# Load the FPKM data
fpkm_data <- read.csv("fpkm_rumen.csv", row.names = 1)
cat("Number of genes in original FPKM data:", nrow(fpkm_data), "\n")

# Remove invalid gene names
fpkm_data <- fpkm_data[!(fpkm_data$gene_name %in% c("-", "", NA)), ]
cat("Number of genes in FPKM data after removing invalid names:", nrow(fpkm_data), "\n")

# Extract gene names separately and process FPKM values only
gene_names_fpkm <- fpkm_data$gene_name  
fpkm_values <- fpkm_data[, -ncol(fpkm_data)]  

# Calculate zFPKM values
zFPKM_matrix <- zFPKM(as.data.frame(fpkm_values))

# Plot zFPKM distribution for visualization
zFPKM_long <- tidyr::pivot_longer(as.data.frame(zFPKM_matrix), cols = everything(), names_to = "Sample", values_to = "zFPKM")
zFPKM_long <- zFPKM_long %>% filter(is.finite(zFPKM))
p <- ggplot(zFPKM_long, aes(x = zFPKM)) +
  geom_histogram(bins = 100, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = -3, color = "red", linetype = "dashed") +  
  labs(title = "Distribution of zFPKM Values", x = "zFPKM", y = "Frequency") +
  theme_minimal()
# Save the plot as JPEG
ggsave("zFPKM_distribution.jpg", 
       plot = p,
       width = 10,  # width in inches
       height = 7,  # height in inches
       dpi = 300)   # resolution

# Apply zFPKM cutoff to select active genes
zFPKM_cutoff <- -3
active_genes <- rowSums(zFPKM_matrix > zFPKM_cutoff) > (ncol(zFPKM_matrix) / 4)
fpkm_filtered <- fpkm_values[active_genes, ]
gene_names_filtered <- gene_names_fpkm[active_genes]

# Create a named vector for gene names to ensure alignment with gene IDs
gene_names_fpkm_named <- setNames(gene_names_filtered, rownames(fpkm_filtered))

cat("Number of genes after zFPKM filtering:", nrow(fpkm_filtered), "\n")

# ------------------------#
# Step 2: Load and Filter Raw Count Data Based on zFPKM Filtering
# ------------------------#

# Load raw count data
raw_count_data <- read.csv("raw count_rumen.csv", row.names = 1)
cat("Number of genes in original raw count data:", nrow(raw_count_data), "\n")

# Filter out invalid gene names
raw_count_data <- raw_count_data[!(raw_count_data$gene_name %in% c("-", "", NA)), ]
cat("Number of genes in raw count data after removing invalid names:", nrow(raw_count_data), "\n")

# Filter raw counts based on zFPKM active genes
filtered_gene_ids <- rownames(fpkm_filtered)  # Genes that passed zFPKM filtering
raw_counts_filtered <- raw_count_data[rownames(raw_count_data) %in% filtered_gene_ids, ]
cat("Number of genes in raw count data after zFPKM filtering:", nrow(raw_counts_filtered), "\n")

# Apply additional filter based on raw count threshold
raw_counts_filtered <- raw_counts_filtered[rowSums(raw_counts_filtered[, -ncol(raw_counts_filtered)]) >= 120, ]
cat("Number of genes after applying raw count threshold:", nrow(raw_counts_filtered), "\n")

# Ensure FPKM data matches the final filtered raw counts
fpkm_final_filtered <- fpkm_filtered[rownames(fpkm_filtered) %in% rownames(raw_counts_filtered), ]

# Retrieve the final gene names for these filtered genes
gene_names_final_filtered <- gene_names_fpkm_named[rownames(fpkm_final_filtered)]

# Check if the lengths of final filtered data match
cat("Number of genes in final filtered FPKM data:", nrow(fpkm_final_filtered), "\n")
cat("Number of genes in final filtered raw count data:", nrow(raw_counts_filtered), "\n")

# Save final filtered FPKM data with gene names
fpkm_final_filtered$gene_name <- gene_names_final_filtered
write.csv(fpkm_final_filtered, "final_filtered_fpkm_values_with_gene_names.csv", row.names = TRUE)
write.csv(raw_counts_filtered, "filtered_raw_counts.csv", row.names = TRUE)

# ------------------------#
# Step 4: DESeq2 Analysis on Filtered Raw Counts
# ------------------------#

# Remove gene names for DESeq2 analysis, as they are not required
raw_counts_filtered <- raw_counts_filtered[, -ncol(raw_counts_filtered)]

# Define experimental conditions for DESeq2
conditions <- factor(c(rep("PC_Insert", 4), rep("Orga_Insert", 4), rep("Tissue", 4)))
col_data <- data.frame(condition = conditions)
rownames(col_data) <- colnames(raw_counts_filtered)

# Create DESeq2 dataset and perform analysis
dds <- DESeqDataSetFromMatrix(countData = raw_counts_filtered, colData = col_data, design = ~ condition)
dds <- DESeq(dds)

# Retrieve DESeq2 results for each comparison
res_pc_vs_tissue <- results(dds, contrast = c("condition", "PC_Insert", "Tissue"), pAdjustMethod = "BH")
res_orga_vs_tissue <- results(dds, contrast = c("condition", "Orga_Insert", "Tissue"), pAdjustMethod = "BH")
res_pc_vs_orga <- results(dds, contrast = c("condition", "PC_Insert", "Orga_Insert"), pAdjustMethod = "BH")

# Convert DESeq2 results to data frames and save
res_pc_vs_tissue_df <- as.data.frame(res_pc_vs_tissue)
res_orga_vs_tissue_df <- as.data.frame(res_orga_vs_tissue)
res_pc_vs_orga_df <- as.data.frame(res_pc_vs_orga)

cat("Number of genes in DESeq2 results (PC vs Tissue):", nrow(res_pc_vs_tissue_df), "\n")
cat("Number of genes in DESeq2 results (Orga vs Tissue):", nrow(res_orga_vs_tissue_df), "\n")
cat("Number of genes in DESeq2 results (PC vs Orga):", nrow(res_pc_vs_orga_df), "\n")

# ------------------------#
# Step 5: Prepare VST-Transformed Data for Downstream Analysis
# ------------------------#

# Perform variance-stabilizing transformation
vst_data <- vst(dds, blind = FALSE)

# Use vst_data for further analysis, correlation studies, and plotting
cat("VST data ready for downstream analysis.\n")

# Save VST-transformed data for use in other scripts or analyses
vst_counts <- assay(vst_data)
write.csv(vst_counts, "vst_transformed_counts.csv", row.names = TRUE)

# ------------------------#
# CORRELATION ANALYSIS AND FILTERING (ALL FILTERED GENES)
# ------------------------#
# Use VST-transformed data for correlation analysis
vst_data_matrix <- assay(vst_data)

# Double-check the dimensions of the VST data matrix
cat("Number of genes in VST data matrix:", nrow(vst_data_matrix), "\n")

# Define group indices for each condition
pc_insert_idx <- which(conditions == "PC_Insert")
tissue_idx <- which(conditions == "Tissue")
orga_insert_idx <- which(conditions == "Orga_Insert")

# Function to calculate correlations between groups
calculate_correlations <- function(data, group1_idx, group2_idx) {
  apply(data, 1, function(x) {
    if (var(x[group1_idx]) == 0 || var(x[group2_idx]) == 0) {
      return(NA)  # Handle cases where variance is zero
    } else {
      return(cor(x[group1_idx], x[group2_idx], method = "pearson"))
    }
  })
}

# Calculate correlations for each comparison using VST data
correlations_pc_vs_tissue <- calculate_correlations(vst_data_matrix, pc_insert_idx, tissue_idx)
correlations_orga_vs_tissue <- calculate_correlations(vst_data_matrix, orga_insert_idx, tissue_idx)
correlations_pc_vs_orga <- calculate_correlations(vst_data_matrix, pc_insert_idx, orga_insert_idx)

# Create correlation data frames for each comparison
correlations_pc_vs_tissue_df <- data.frame(
  gene_id = rownames(vst_data_matrix), 
  gene_name = gene_names_final_filtered, 
  correlation = correlations_pc_vs_tissue
) %>% filter(!is.na(correlation))

correlations_orga_vs_tissue_df <- data.frame(
  gene_id = rownames(vst_data_matrix), 
  gene_name = gene_names_final_filtered, 
  correlation = correlations_orga_vs_tissue
) %>% filter(!is.na(correlation))

correlations_pc_vs_orga_df <- data.frame(
  gene_id = rownames(vst_data_matrix), 
  gene_name = gene_names_final_filtered, 
  correlation = correlations_pc_vs_orga
) %>% filter(!is.na(correlation))

# Print the number of genes with valid correlations for each comparison
cat("Number of genes with valid correlations (PC vs Tissue):", nrow(correlations_pc_vs_tissue_df), "\n")
cat("Number of genes with valid correlations (Orga vs Tissue):", nrow(correlations_orga_vs_tissue_df), "\n")
cat("Number of genes with valid correlations (PC vs Orga):", nrow(correlations_pc_vs_orga_df), "\n")

# Save the correlation results for reference
write.csv(correlations_pc_vs_tissue_df, "correlations_PC_vs_Tissue.csv", row.names = FALSE)
write.csv(correlations_orga_vs_tissue_df, "correlations_Orga_vs_Tissue.csv", row.names = FALSE)
write.csv(correlations_pc_vs_orga_df, "correlations_PC_vs_Orga.csv", row.names = FALSE)

# ------------------------#
# MERGE CORRELATION RESULTS WITH DESEQ2 RESULTS
# ------------------------#

# Merge the correlation data frames with DESeq2 results to include log2FoldChange
merged_correlations_pc_vs_tissue <- merge(correlations_pc_vs_tissue_df, res_pc_vs_tissue_df, by.x = "gene_id", by.y = "row.names", all.x = TRUE)
merged_correlations_orga_vs_tissue <- merge(correlations_orga_vs_tissue_df, res_orga_vs_tissue_df, by.x = "gene_id", by.y = "row.names", all.x = TRUE)
merged_correlations_pc_vs_orga <- merge(correlations_pc_vs_orga_df, res_pc_vs_orga_df, by.x = "gene_id", by.y = "row.names", all.x = TRUE)

# Remove rows with NA values in the padj column to retain only significant genes
merged_correlations_pc_vs_tissue <- merged_correlations_pc_vs_tissue[!is.na(merged_correlations_pc_vs_tissue$padj), ]
merged_correlations_orga_vs_tissue <- merged_correlations_orga_vs_tissue[!is.na(merged_correlations_orga_vs_tissue$padj), ]
merged_correlations_pc_vs_orga <- merged_correlations_pc_vs_orga[!is.na(merged_correlations_pc_vs_orga$padj), ]

# Save merged correlation data frames
write.csv(merged_correlations_pc_vs_tissue, "merged_correlations_PC_vs_Tissue.csv", row.names = FALSE)
write.csv(merged_correlations_orga_vs_tissue, "merged_correlations_Orga_vs_Tissue.csv", row.names = FALSE)
write.csv(merged_correlations_pc_vs_orga, "merged_correlations_PC_vs_Orga.csv", row.names = FALSE)

# Print the number of genes after merging with DESeq2 results for each comparison
cat("Number of genes after merging with DESeq2 (PC vs Tissue):", nrow(merged_correlations_pc_vs_tissue), "\n")
cat("Number of genes after merging with DESeq2 (Orga vs Tissue):", nrow(merged_correlations_orga_vs_tissue), "\n")
cat("Number of genes after merging with DESeq2 (PC vs Orga):", nrow(merged_correlations_pc_vs_orga), "\n")

# ------------------------#
# FILTER FOR HIGHLY CORRELATED GENES
# ------------------------#

# Set correlation threshold for filtering
correlation_threshold <- 0.7

# Filter for highly correlated genes (|correlation| > correlation_threshold) for each group by group comparison
highly_correlated_pc_vs_tissue <- merged_correlations_pc_vs_tissue %>%
  filter(abs(correlation) > correlation_threshold)

highly_correlated_orga_vs_tissue <- merged_correlations_orga_vs_tissue %>%
  filter(abs(correlation) > correlation_threshold)

highly_correlated_pc_vs_orga <- merged_correlations_pc_vs_orga %>%
  filter(abs(correlation) > correlation_threshold)

# Print the number of genes after filtering for high correlation
cat("Number of highly correlated genes (PC vs Tissue):", nrow(highly_correlated_pc_vs_tissue), "\n")
cat("Number of highly correlated genes (Orga vs Tissue):", nrow(highly_correlated_orga_vs_tissue), "\n")
cat("Number of highly correlated genes (PC vs Orga):", nrow(highly_correlated_pc_vs_orga), "\n")

# Save filtered highly correlated genes data frames
write.csv(highly_correlated_pc_vs_tissue, "filtered_highly_correlated_genes_PC_vs_Tissue.csv", row.names = FALSE)
write.csv(highly_correlated_orga_vs_tissue, "filtered_highly_correlated_genes_Orga_vs_Tissue.csv", row.names = FALSE)
write.csv(highly_correlated_pc_vs_orga, "filtered_highly_correlated_genes_PC_vs_Orga.csv", row.names = FALSE)

# ------------------------#
# PREPARE GENE LISTS FOR GSEA AFTER APPLYING CORRELATION FILTER
# ------------------------#

# Rank genes by log2 fold change after filtering for high correlation
geneList_pc_vs_tissue <- highly_correlated_pc_vs_tissue$log2FoldChange
names(geneList_pc_vs_tissue) <- highly_correlated_pc_vs_tissue$gene_id
geneList_pc_vs_tissue <- sort(geneList_pc_vs_tissue, decreasing = TRUE)

geneList_orga_vs_tissue <- highly_correlated_orga_vs_tissue$log2FoldChange
names(geneList_orga_vs_tissue) <- highly_correlated_orga_vs_tissue$gene_id
geneList_orga_vs_tissue <- sort(geneList_orga_vs_tissue, decreasing = TRUE)

geneList_pc_vs_orga <- highly_correlated_pc_vs_orga$log2FoldChange
names(geneList_pc_vs_orga) <- highly_correlated_pc_vs_orga$gene_id
geneList_pc_vs_orga <- sort(geneList_pc_vs_orga, decreasing = TRUE)

# Print the number of genes in each final gene list for GSEA analysis
cat("Number of genes in geneList for GSEA (PC vs Tissue):", length(geneList_pc_vs_tissue), "\n")
cat("Number of genes in geneList for GSEA (Orga vs Tissue):", length(geneList_orga_vs_tissue), "\n")
cat("Number of genes in geneList for GSEA (PC vs Orga):", length(geneList_pc_vs_orga), "\n")


# ================================================== #
# COMPLETE GSEA ANALYSIS AND VISUALIZATION WORKFLOW
# ================================================== #

# Set global options
options(ggrepel.max.overlaps = Inf)

# ========================= #
# 1. CORE FUNCTION DEFINITIONS
# ========================= #

# Function to run individual GSEA with proper simplification sequence
run_individual_gsea <- function(gene_list, analysis_type, ont = NULL) {
  if(is.null(gene_list) || length(gene_list) == 0) {
    message(sprintf("Invalid gene list provided for %s analysis", analysis_type))
    return(NULL)
  }
  
  tryCatch({
    message(sprintf("\nStarting %s analysis...", analysis_type))
    if(analysis_type == "KEGG") {
      result <- gseKEGG(
        geneList = jitter(gene_list, factor = 1e-10),
        organism = "oas",
        minGSSize = 6,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        eps = 0,
        BPPARAM = SnowParam(exportglobals = FALSE)
      )
    } else {
      result <- gseGO(
        geneList = jitter(gene_list, factor = 1e-10),
        OrgDb = sheep_orgdb,
        keyType = "ENTREZID",
        ont = ont,
        minGSSize = 6,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        eps = 0,
        BPPARAM = SnowParam(exportglobals = FALSE)
      )
      
      # Simplify GO results immediately after getting them
      if(!is.null(result) && nrow(result@result) > 0) {
        message(sprintf("Simplifying %s pathways...", analysis_type))
        result <- clusterProfiler::simplify(result,
                                            cutoff = 0.7,
                                            by = "p.adjust",
                                            select_fun = min,
                                            measure = "Wang",
                                            semData = NULL)
      }
    }
    
    # Check if any enrichment was found
    if(is.null(result) || nrow(result@result) == 0) {
      message(sprintf("No enrichment found for %s analysis", analysis_type))
      return(NULL)
    }
    
    message(sprintf("Found %d enriched terms for %s analysis", nrow(result@result), analysis_type))
    return(result)
    
  }, error = function(e) {
    message(sprintf("Error in %s analysis: %s", analysis_type, e$message))
    return(NULL)
  })
}

# Function to select balanced pathways
select_balanced_pathways <- function(gsea_result) {
  if(is.null(gsea_result)) {
    message("No GSEA result provided")
    return(NULL)
  }
  
  if(nrow(gsea_result@result) == 0) {
    message("No significant pathways found")
    return(NULL)
  }
  
  tryCatch({
    result_df <- as.data.frame(gsea_result@result)
    
    # Base exclusion pattern
    exclude_pattern <- paste(
      "Folate|COVID|group of donors|Staphylo|Diabetic|Parkinson|papilloma|",
      "Malaria|fatty liver|Huntington|muscle|cytomegalo|diabetic|Drug|",
      "carcinoma|arthritis|cancer|carcinogenesis|cardio|virus|Axon|Viral|viral|",
      "xenobiotics|Relaxin|Steroid|disease|Disease|Alzheimer|Tuberculosis|", 
      "lupus|polysaccharide|MHC|Platelet|antigen|Intestinal|IgA|Th1|Th17|",
      "Allograft|Toxoplasmosis|toxoplasmosis|Insulin|insulin|secretion|plasmosis|Yersinia|infection|",
      sep=""
    )
    
    # Filter pathways
    initial_count <- nrow(result_df)
    result_df <- result_df[!grepl(exclude_pattern, result_df$Description), ]
    message(sprintf("Filtered out %d pathways based on exclusion criteria", initial_count - nrow(result_df)))
    
    if(nrow(result_df) == 0) {
      message("No pathways remaining after filtering")
      return(NULL)
    }
    
    # Split into activated and suppressed pathways
    activated <- result_df[result_df$NES > 0, ]
    suppressed <- result_df[result_df$NES < 0, ]
    
    message("\nInitial pathway counts:")
    message(sprintf("Activated pathways: %d", nrow(activated)))
    message(sprintf("Suppressed pathways: %d", nrow(suppressed)))
    
    # Calculate balanced selection
    n_activated <- min(4, nrow(activated))
    n_suppressed <- min(4, nrow(suppressed))
    
    # Adjust counts for imbalanced cases
    if(n_activated < 4 && nrow(suppressed) > 4) {
      n_suppressed <- min(8 - n_activated, nrow(suppressed))
    }
    if(n_suppressed < 4 && nrow(activated) > 4) {
      n_activated <- min(8 - n_suppressed, nrow(activated))
    }
    
    message("\nFinal selection:")
    message(sprintf("Selected activated pathways: %d", n_activated))
    message(sprintf("Selected suppressed pathways: %d", n_suppressed))
    
    selected_activated <- head(activated$ID, n_activated)
    selected_suppressed <- head(suppressed$ID, n_suppressed)
    
    return(c(selected_activated, selected_suppressed))
    
  }, error = function(e) {
    message(sprintf("Error in pathway selection: %s", e$message))
    return(NULL)
  })
}

# Function to check and validate gene list
validate_gene_list <- function(gene_list, name = "") {
  if(is.null(gene_list)) {
    message(sprintf("Gene list %sis NULL", ifelse(name == "", "", paste(name, " "))))
    return(FALSE)
  }
  
  if(length(gene_list) == 0) {
    message(sprintf("Gene list %sis empty", ifelse(name == "", "", paste(name, " "))))
    return(FALSE)
  }
  
  if(is.null(names(gene_list))) {
    message(sprintf("Gene list %slacks names", ifelse(name == "", "", paste(name, " "))))
    return(FALSE)
  }
  
  return(TRUE)
}

# Function to check directory structure
check_and_create_dirs <- function() {
  tryCatch({
    dir.create("GSEA_Results", showWarnings = FALSE)
    dir.create(file.path("GSEA_Results", "Plots"), showWarnings = FALSE)
    dir.create(file.path("GSEA_Results", "Data"), showWarnings = FALSE)
    message("Directory structure created successfully")
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error creating directory structure: %s", e$message))
    return(FALSE)
  })
}

# ========================= #
# 2. PATHWAY ANALYSIS AND RESULTS SAVING FUNCTIONS
# ========================= #

# Function to save detailed pathway information 
save_pathway_details <- function(gsea_result, type, comparison, output_dir = ".") {
  if(is.null(gsea_result)) {
    message(sprintf("No GSEA result provided for %s %s", type, comparison))
    return(NULL)
  }
  
  if(nrow(gsea_result@result) == 0) {
    message(sprintf("No enriched pathways found for %s %s", type, comparison))
    return(NULL)
  }
  
  tryCatch({
    # Create main results dataframe
    main_results <- as.data.frame(gsea_result@result)
    
    # Create detailed gene info with error checking for each pathway
    gene_details_list <- lapply(1:nrow(main_results), function(i) {
      pathway_id <- main_results$ID[i]
      
      if(is.null(main_results$core_enrichment[i]) || 
         is.na(main_results$core_enrichment[i]) || 
         main_results$core_enrichment[i] == "") {
        message(sprintf("Warning: No core enrichment genes for pathway %s", pathway_id))
        return(NULL)
      }
      
      tryCatch({
        core_genes <- strsplit(main_results$core_enrichment[i], "/")[[1]]
        
        # Check if gene names are available
        gene_names <- gene_names_fpkm_named[core_genes]
        if(all(is.na(gene_names))) {
          message(sprintf("Warning: No gene names found for pathway %s", pathway_id))
        }
        
        data.frame(
          Pathway_ID = pathway_id,
          Pathway_Name = main_results$Description[i],
          Gene_ID = core_genes,
          Gene_Name = gene_names,
          NES = main_results$NES[i],
          pvalue = main_results$pvalue[i],
          p.adjust = main_results$p.adjust[i],
          qvalue = main_results$qvalue[i],
          Leading_Edge = main_results$leading_edge[i],
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        message(sprintf("Error processing pathway %s: %s", pathway_id, e$message))
        return(NULL)
      })
    })
    
    # Remove NULL entries and combine results
    gene_details <- do.call(rbind, Filter(Negate(is.null), gene_details_list))
    
    if(is.null(gene_details) || nrow(gene_details) == 0) {
      message(sprintf("No valid gene details to save for %s %s", type, comparison))
      return(NULL)
    }
    
    # Create directory if it doesn't exist
    dir.create(file.path(output_dir, "GSEA_Results", "Data"), recursive = TRUE, showWarnings = FALSE)
    
    # Save files
    filename_base <- file.path(output_dir, "GSEA_Results", "Data", 
                               paste0("GSEA_", type, "_", comparison))
    
    write.csv(main_results, paste0(filename_base, "_pathway_summary.csv"), row.names = FALSE)
    write.csv(gene_details, paste0(filename_base, "_gene_details.csv"), row.names = FALSE)
    
    message(sprintf("Saved detailed results for %s %s", type, comparison))
    
    # Return results invisibly
    return(invisible(list(
      summary = main_results,
      details = gene_details,
      stats = list(
        total_pathways = nrow(main_results),
        total_genes = length(unique(gene_details$Gene_ID)),
        named_genes = sum(!is.na(gene_details$Gene_Name))
      )
    )))
    
  }, error = function(e) {
    message(sprintf("Error saving pathway details for %s %s: %s", type, comparison, e$message))
    return(NULL)
  })
}

# Function to analyze genes in pathways 
analyze_gsea_genes <- function(gsea_result, pathway_id, geneList) {
  if(!validate_gene_list(geneList)) return(NULL)
  
  tryCatch({
    idx <- which(gsea_result@result$ID == pathway_id)
    if(length(idx) == 0) {
      message(sprintf("Pathway %s not found in results", pathway_id))
      return(NULL)
    }
    
    message("\n=== ANALYZING GSEA PATHWAY GENES ===")
    message(sprintf("Pathway: %s", pathway_id))
    message(sprintf("Description: %s", gsea_result@result$Description[idx]))
    message(sprintf("SetSize: %d", gsea_result@result$setSize[idx]))
    
    # Get and validate core enrichment genes
    if(is.null(gsea_result@result$core_enrichment[idx]) || 
       is.na(gsea_result@result$core_enrichment[idx])) {
      message("No core enrichment genes found")
      return(NULL)
    }
    
    core_genes <- strsplit(gsea_result@result$core_enrichment[idx], "/")[[1]]
    
    # Get genes in ranked list
    ranked_genes <- names(geneList)
    genes_in_analysis <- intersect(core_genes, ranked_genes)
    
    if(length(genes_in_analysis) == 0) {
      message("No genes found in analysis")
      return(NULL)
    }
    
    # Calculate and display gene statistics
    message("\nGene Statistics:")
    message(sprintf("Total genes in pathway: %d", length(core_genes)))
    message(sprintf("Genes found in ranked list: %d", length(genes_in_analysis)))
    
    # Get fold changes for found genes
    fc_values <- geneList[genes_in_analysis]
    
    # Sort genes by absolute fold change
    sorted_idx <- order(abs(fc_values), decreasing = TRUE)
    sorted_genes <- genes_in_analysis[sorted_idx]
    sorted_fc <- fc_values[sorted_idx]
    
    message("\nGenes ordered by absolute fold change:")
    for(i in seq_along(sorted_genes)) {
      message(sprintf("%s (%s): %.3f", 
                      sorted_genes[i],
                      gene_names_fpkm_named[sorted_genes[i]],
                      sorted_fc[i]))
    }
    
    # Return detailed results
    return(invisible(list(
      genes = sorted_genes,
      fold_changes = sorted_fc,
      gene_names = gene_names_fpkm_named[sorted_genes],
      stats = list(
        total_genes = length(core_genes),
        found_genes = length(genes_in_analysis),
        max_fc = max(abs(sorted_fc)),
        min_fc = min(abs(sorted_fc))
      )
    )))
    
  }, error = function(e) {
    message(sprintf("Error analyzing genes for pathway %s: %s", pathway_id, e$message))
    return(NULL)
  })
}

# ========================= #
# 3. VISUALIZATION FUNCTIONS
# ========================= #
# Helper function to merge GSEA results while preserving activation status
merge_gsea_results <- function(gsea_list) {
  # Combine all results into a single data frame
  result_df <- do.call(rbind, lapply(names(gsea_list), function(name) {
    result <- gsea_list[[name]]@result
    result$Comparison <- name
    result$.sign <- ifelse(result$NES > 0, "Activated", "Suppressed")
    result
  }))
  
  # Create a new gseaResult object with the combined results
  merged <- gsea_list[[1]]
  merged@result <- result_df
  return(merged)
}

# Function to process and plot dotplots 
process_and_plot <- function(gsea_list, title) {
  # Remove NULL or empty results
  gsea_list <- Filter(function(x) !is.null(x) && nrow(x@result) > 0, gsea_list)
  
  if(length(gsea_list) == 0) {
    message(sprintf("No valid results for %s after filtering empty results", title))
    return(NULL)
  }
  
  tryCatch({
    processed_list <- list()
    for(name in names(gsea_list)) {
      message(sprintf("\nProcessing %s for %s", name, title))
      selected_paths <- select_balanced_pathways(gsea_list[[name]])  
      if(!is.null(selected_paths) && length(selected_paths) > 0) {
        current_result <- gsea_list[[name]]
        current_result@result <- current_result@result[current_result@result$ID %in% selected_paths, ]
        processed_list[[name]] <- current_result
      }
    }
    
    if(length(processed_list) == 0) {
      message(sprintf("No pathways remain after filtering for %s", title))
      return(NULL)
    }
    
    merged_result <- merge_gsea_results(processed_list)
    
    if(nrow(merged_result@result) == 0) {
      message(sprintf("No results to plot for %s after merging", title))
      return(NULL)
    }
    
    # Create dotplot with modified layout
    p <- dotplot(merged_result,
                 showCategory = nrow(merged_result@result),
                 split = ".sign") +
      facet_grid(.~.sign) +  # Split by activation status
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "gray95"),
        legend.position = "right",
        legend.box = "vertical",
        axis.text.y = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "cm")
      ) +
      ggtitle(title) +
      # Move x-axis to Comparison and GeneRatio to legend
      aes(x = Comparison) +
      guides(
        size = guide_legend(title = "GeneRatio", order = 1),
        color = guide_colorbar(title = "p.adjust", order = 2)
      )
    
    # Save plot
    plot_path <- file.path("GSEA_Results", "Plots", paste0("GSEA_", gsub(" ", "_", title), "_dotplot.jpeg"))
    jpeg(plot_path, width = 15*300, height = 8*300, res = 300)
    print(p)
    dev.off()
    message(sprintf("Saved dotplot: %s", plot_path))
    
    # Display in R
    print(p)
    
    return(merged_result)
    
  }, error = function(e) {
    message(sprintf("Error in dotplot creation for %s: %s", title, e$message))
    message("\nDebug information:")
    message("Number of valid results: ", length(gsea_list))
    for(name in names(gsea_list)) {
      message(sprintf("%s: %d pathways", name, 
                      ifelse(is.null(gsea_list[[name]]), 0, 
                             nrow(gsea_list[[name]]@result))))
    }
    return(NULL)
  })
}

# Function to create cnetplots
create_cnetplot <- function(gsea_result, title, gene_list) {
  if(is.null(gsea_result)) {
    message(sprintf("No GSEA results available for %s", title))
    return(NULL)
  }
  
  if(nrow(gsea_result@result) == 0) {
    message(sprintf("No enriched pathways found for %s", title))
    return(NULL)
  }
  
  tryCatch({
    selected_paths <- select_balanced_pathways(gsea_result)
    if(is.null(selected_paths)) {
      message(sprintf("No pathways selected for cnetplot: %s", title))
      return(NULL)
    }
    
    # Print pathway info
    message(sprintf("\nSelected pathways for %s:", title))
    result_df <- as.data.frame(gsea_result)
    selected_df <- result_df[result_df$ID %in% selected_paths, ]
    print(selected_df[, c("ID", "Description", "setSize", "NES")])
    
    # Analyze genes for each pathway
    message("\nAnalyzing genes in pathways:")
    for(path_id in selected_paths) {
      analyze_gsea_genes(gsea_result, path_id, gene_list)
    }
    
    # Make gene IDs readable
    result_readable <- setReadable(gsea_result, OrgDb = sheep_orgdb, 'ENTREZID')
    result_subset <- result_readable
    result_subset@result <- result_readable@result[result_readable@result$ID %in% selected_paths, ]
    
    # Create plot
    p <- cnetplot(result_subset,
                  node_label = "all",
                  showCategory = length(selected_paths),
                  cex.params = list(
                    category_node = 1.5,
                    gene_node = 0.4,
                    category_label = 1.2,
                    gene_label = 0.5
                  ),
                  max.overlaps = Inf,
                  layout = "dh",
                  node.shape = "sphere",
                  edge.width = 0.5,
                  node.color.alpha = 0.7,
                  color.params = list(foldChange = gene_list)) +
      scale_color_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        name = "Log2 Fold Change"
      ) +
      theme(legend.position = "right") +
      ggtitle(title)
    
    return(p)
    
  }, error = function(e) {
    message(sprintf("Error creating cnetplot for %s: %s", title, e$message))
    return(NULL)
  })
}

# Function to save cnetplots 
save_cnetplot <- function(plot, filename) {
  if(is.null(plot)) {
    message(sprintf("No plot to save for %s", filename))
    return(FALSE)
  }
  
  tryCatch({
    plot_path <- file.path("GSEA_Results", "Plots", filename)
    jpeg(plot_path, width = 15*300, height = 12*300, res = 300)
    print(plot)
    dev.off()
    message(sprintf("Saved cnetplot: %s", plot_path))
    print(plot)  # Display in R
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error saving cnetplot %s: %s", filename, e$message))
    return(FALSE)
  })
}

# ========================= #
# 4. ANALYSIS EXECUTION
# ========================= #

# Create directory structure
check_and_create_dirs()

# Run GSEA for each comparison and analysis type
message("\n=== Running GSEA Analyses ===")

# BP analyses
message("\nRunning BP analyses...")
gsea_BP_pc_tissue <- run_individual_gsea(geneList_pc_vs_tissue, "GO", "BP")
gsea_BP_orga_tissue <- run_individual_gsea(geneList_orga_vs_tissue, "GO", "BP")
gsea_BP_pc_orga <- run_individual_gsea(geneList_pc_vs_orga, "GO", "BP")

# Save BP results
save_pathway_details(gsea_BP_pc_tissue, "BP", "PC_vs_Tissue")
save_pathway_details(gsea_BP_orga_tissue, "BP", "Orga_vs_Tissue")
save_pathway_details(gsea_BP_pc_orga, "BP", "PC_vs_Orga")

# CC analyses
message("\nRunning CC analyses...")
gsea_CC_pc_tissue <- run_individual_gsea(geneList_pc_vs_tissue, "GO", "CC")
gsea_CC_orga_tissue <- run_individual_gsea(geneList_orga_vs_tissue, "GO", "CC")
gsea_CC_pc_orga <- run_individual_gsea(geneList_pc_vs_orga, "GO", "CC")

# Save CC results
save_pathway_details(gsea_CC_pc_tissue, "CC", "PC_vs_Tissue")
save_pathway_details(gsea_CC_orga_tissue, "CC", "Orga_vs_Tissue")
save_pathway_details(gsea_CC_pc_orga, "CC", "PC_vs_Orga")

# MF analyses
message("\nRunning MF analyses...")
gsea_MF_pc_tissue <- run_individual_gsea(geneList_pc_vs_tissue, "GO", "MF")
gsea_MF_orga_tissue <- run_individual_gsea(geneList_orga_vs_tissue, "GO", "MF")
gsea_MF_pc_orga <- run_individual_gsea(geneList_pc_vs_orga, "GO", "MF")

# Save MF results
save_pathway_details(gsea_MF_pc_tissue, "MF", "PC_vs_Tissue")
save_pathway_details(gsea_MF_orga_tissue, "MF", "Orga_vs_Tissue")
save_pathway_details(gsea_MF_pc_orga, "MF", "PC_vs_Orga")

# KEGG analyses
message("\nRunning KEGG analyses...")
gsea_KEGG_pc_tissue <- run_individual_gsea(geneList_pc_vs_tissue, "KEGG")
gsea_KEGG_orga_tissue <- run_individual_gsea(geneList_orga_vs_tissue, "KEGG")
gsea_KEGG_pc_orga <- run_individual_gsea(geneList_pc_vs_orga, "KEGG")

# Save KEGG results
save_pathway_details(gsea_KEGG_pc_tissue, "KEGG", "PC_vs_Tissue")
save_pathway_details(gsea_KEGG_orga_tissue, "KEGG", "Orga_vs_Tissue")
save_pathway_details(gsea_KEGG_pc_orga, "KEGG", "PC_vs_Orga")

# Create lists for merging results
message("\nPreparing result lists for visualization...")

# Create original lists
gsea_list_BP <- list(
  "Orga_vs_Tissue" = gsea_BP_orga_tissue,
  "PC_vs_Tissue" = gsea_BP_pc_tissue,
  "PC_vs_Orga" = gsea_BP_pc_orga
)

gsea_list_CC <- list(
  "Orga_vs_Tissue" = gsea_CC_orga_tissue,
  "PC_vs_Tissue" = gsea_CC_pc_tissue,
  "PC_vs_Orga" = gsea_CC_pc_orga
)

gsea_list_MF <- list(
  "Orga_vs_Tissue" = gsea_MF_orga_tissue,
  "PC_vs_Tissue" = gsea_MF_pc_tissue,
  "PC_vs_Orga" = gsea_MF_pc_orga
)

gsea_list_KEGG <- list(
  "Orga_vs_Tissue" = gsea_KEGG_orga_tissue,
  "PC_vs_Tissue" = gsea_KEGG_pc_tissue,
  "PC_vs_Orga" = gsea_KEGG_pc_orga
)

# Create visualizations
message("\n=== Creating Visualizations ===")

# Create dotplots
message("\nCreating dotplots...")
bp_dot <- process_and_plot(gsea_list_BP, "GO Biological Process")
cc_dot <- process_and_plot(gsea_list_CC, "GO Cellular Component")
mf_dot <- process_and_plot(gsea_list_MF, "GO Molecular Function")
kegg_dot <- process_and_plot(gsea_list_KEGG, "KEGG Pathways")

# Create and save cnetplots
message("\nCreating cnetplots...")

# BP cnetplots
message("\nProcessing BP cnetplots...")
p_cnet_BP_orga_tissue <- create_cnetplot(gsea_BP_orga_tissue, "BP Orga vs Tissue", geneList_orga_vs_tissue)
p_cnet_BP_pc_tissue <- create_cnetplot(gsea_BP_pc_tissue, "BP PC vs Tissue", geneList_pc_vs_tissue)
p_cnet_BP_pc_orga <- create_cnetplot(gsea_BP_pc_orga, "BP PC vs Orga", geneList_pc_vs_orga)

save_cnetplot(p_cnet_BP_orga_tissue, "GSEA_cnetplot_BP_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_BP_pc_tissue, "GSEA_cnetplot_BP_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_BP_pc_orga, "GSEA_cnetplot_BP_PC_vs_Orga.jpeg")

# CC cnetplots
message("\nProcessing CC cnetplots...")
p_cnet_CC_orga_tissue <- create_cnetplot(gsea_CC_orga_tissue, "CC Orga vs Tissue", geneList_orga_vs_tissue)
p_cnet_CC_pc_tissue <- create_cnetplot(gsea_CC_pc_tissue, "CC PC vs Tissue", geneList_pc_vs_tissue)
p_cnet_CC_pc_orga <- create_cnetplot(gsea_CC_pc_orga, "CC PC vs Orga", geneList_pc_vs_orga)

save_cnetplot(p_cnet_CC_orga_tissue, "GSEA_cnetplot_CC_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_CC_pc_tissue, "GSEA_cnetplot_CC_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_CC_pc_orga, "GSEA_cnetplot_CC_PC_vs_Orga.jpeg")

# MF cnetplots
message("\nProcessing MF cnetplots...")
p_cnet_MF_orga_tissue <- create_cnetplot(gsea_MF_orga_tissue, "MF Orga vs Tissue", geneList_orga_vs_tissue)
p_cnet_MF_pc_tissue <- create_cnetplot(gsea_MF_pc_tissue, "MF PC vs Tissue", geneList_pc_vs_tissue)
p_cnet_MF_pc_orga <- create_cnetplot(gsea_MF_pc_orga, "MF PC vs Orga", geneList_pc_vs_orga)

save_cnetplot(p_cnet_MF_orga_tissue, "GSEA_cnetplot_MF_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_MF_pc_tissue, "GSEA_cnetplot_MF_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_MF_pc_orga, "GSEA_cnetplot_MF_PC_vs_Orga.jpeg")

# KEGG cnetplots
message("\nProcessing KEGG cnetplots...")
p_cnet_KEGG_orga_tissue <- create_cnetplot(gsea_KEGG_orga_tissue, "KEGG Orga vs Tissue", geneList_orga_vs_tissue)
p_cnet_KEGG_pc_tissue <- create_cnetplot(gsea_KEGG_pc_tissue, "KEGG PC vs Tissue", geneList_pc_vs_tissue)
p_cnet_KEGG_pc_orga <- create_cnetplot(gsea_KEGG_pc_orga, "KEGG PC vs Orga", geneList_pc_vs_orga)

save_cnetplot(p_cnet_KEGG_orga_tissue, "GSEA_cnetplot_KEGG_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_KEGG_pc_tissue, "GSEA_cnetplot_KEGG_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_KEGG_pc_orga, "GSEA_cnetplot_KEGG_PC_vs_Orga.jpeg")

message("\nAll analyses and visualizations complete!")
message("Results have been saved in the GSEA_Results directory")

# Print final summary
message("\n=== Final Analysis Summary ===")
message("\nResults directories created:")
message("- GSEA_Results/Plots: Contains all visualization files")
message("- GSEA_Results/Data: Contains all data files")


# ------------------------#
# PREPARE DEG-BASED GENE LISTS FOR GSEA
# ------------------------#

message("\n=== Preparing DEG-Based Gene Lists for GSEA ===")

# Function to filter DEGs
filter_degs <- function(res, comparison_name, padj_cutoff = 0.05, log2fc_cutoff = 1) {
  filtered_res <- as.data.frame(res) %>%
    filter(!is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) > log2fc_cutoff)
  
  # Create gene list directly matching your correlation approach
  gene_list <- filtered_res$log2FoldChange
  names(gene_list) <- rownames(filtered_res)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  message(sprintf("\nDEG Analysis for %s:", comparison_name))
  message(sprintf("Number of significant DEGs (padj < %s, |log2FC| > %s): %d", 
                  padj_cutoff, log2fc_cutoff, length(gene_list)))
  
  return(gene_list)
}

# Prepare gene lists for each comparison
geneList_pc_vs_tissue_deg <- filter_degs(res_pc_vs_tissue, "PC vs Tissue")
geneList_orga_vs_tissue_deg <- filter_degs(res_orga_vs_tissue, "Orga vs Tissue")
geneList_pc_vs_orga_deg <- filter_degs(res_pc_vs_orga, "PC vs Orga")

# ================================================== #
# DEG-BASED GSEA ANALYSIS AND VISUALIZATION WORKFLOW
# ================================================== #

message("\n=== Running DEG-Based GSEA Analyses ===")

# BP analyses for DEGs
message("\nRunning BP analyses for DEGs...")
gsea_BP_pc_tissue_deg <- run_individual_gsea(geneList_pc_vs_tissue_deg, "GO", "BP")
gsea_BP_orga_tissue_deg <- run_individual_gsea(geneList_orga_vs_tissue_deg, "GO", "BP")
gsea_BP_pc_orga_deg <- run_individual_gsea(geneList_pc_vs_orga_deg, "GO", "BP")

# Save BP results for DEGs
save_pathway_details(gsea_BP_pc_tissue_deg, "BP_DEG", "PC_vs_Tissue")
save_pathway_details(gsea_BP_orga_tissue_deg, "BP_DEG", "Orga_vs_Tissue")
save_pathway_details(gsea_BP_pc_orga_deg, "BP_DEG", "PC_vs_Orga")

# CC analyses for DEGs
message("\nRunning CC analyses for DEGs...")
gsea_CC_pc_tissue_deg <- run_individual_gsea(geneList_pc_vs_tissue_deg, "GO", "CC")
gsea_CC_orga_tissue_deg <- run_individual_gsea(geneList_orga_vs_tissue_deg, "GO", "CC")
gsea_CC_pc_orga_deg <- run_individual_gsea(geneList_pc_vs_orga_deg, "GO", "CC")

# Save CC results for DEGs
save_pathway_details(gsea_CC_pc_tissue_deg, "CC_DEG", "PC_vs_Tissue")
save_pathway_details(gsea_CC_orga_tissue_deg, "CC_DEG", "Orga_vs_Tissue")
save_pathway_details(gsea_CC_pc_orga_deg, "CC_DEG", "PC_vs_Orga")

# MF analyses for DEGs
message("\nRunning MF analyses for DEGs...")
gsea_MF_pc_tissue_deg <- run_individual_gsea(geneList_pc_vs_tissue_deg, "GO", "MF")
gsea_MF_orga_tissue_deg <- run_individual_gsea(geneList_orga_vs_tissue_deg, "GO", "MF")
gsea_MF_pc_orga_deg <- run_individual_gsea(geneList_pc_vs_orga_deg, "GO", "MF")

# Save MF results for DEGs
save_pathway_details(gsea_MF_pc_tissue_deg, "MF_DEG", "PC_vs_Tissue")
save_pathway_details(gsea_MF_orga_tissue_deg, "MF_DEG", "Orga_vs_Tissue")
save_pathway_details(gsea_MF_pc_orga_deg, "MF_DEG", "PC_vs_Orga")

# KEGG analyses for DEGs
message("\nRunning KEGG analyses for DEGs...")
gsea_KEGG_pc_tissue_deg <- run_individual_gsea(geneList_pc_vs_tissue_deg, "KEGG")
gsea_KEGG_orga_tissue_deg <- run_individual_gsea(geneList_orga_vs_tissue_deg, "KEGG")
gsea_KEGG_pc_orga_deg <- run_individual_gsea(geneList_pc_vs_orga_deg, "KEGG")

# Save KEGG results for DEGs
save_pathway_details(gsea_KEGG_pc_tissue_deg, "KEGG_DEG", "PC_vs_Tissue")
save_pathway_details(gsea_KEGG_orga_tissue_deg, "KEGG_DEG", "Orga_vs_Tissue")
save_pathway_details(gsea_KEGG_pc_orga_deg, "KEGG_DEG", "PC_vs_Orga")

# Create lists for merging results
message("\nPreparing DEG-based result lists for visualization...")

# Create original lists for DEG-based results
gsea_list_BP_deg <- list(
  "Orga_vs_Tissue" = gsea_BP_orga_tissue_deg,
  "PC_vs_Tissue" = gsea_BP_pc_tissue_deg,
  "PC_vs_Orga" = gsea_BP_pc_orga_deg
)

gsea_list_CC_deg <- list(
  "Orga_vs_Tissue" = gsea_CC_orga_tissue_deg,
  "PC_vs_Tissue" = gsea_CC_pc_tissue_deg,
  "PC_vs_Orga" = gsea_CC_pc_orga_deg
)

gsea_list_MF_deg <- list(
  "Orga_vs_Tissue" = gsea_MF_orga_tissue_deg,
  "PC_vs_Tissue" = gsea_MF_pc_tissue_deg,
  "PC_vs_Orga" = gsea_MF_pc_orga_deg
)

gsea_list_KEGG_deg <- list(
  "Orga_vs_Tissue" = gsea_KEGG_orga_tissue_deg,
  "PC_vs_Tissue" = gsea_KEGG_pc_tissue_deg,
  "PC_vs_Orga" = gsea_KEGG_pc_orga_deg
)

# Create visualizations
message("\n=== Creating DEG-Based Visualizations ===")

# Create dotplots
message("\nCreating DEG-based dotplots...")
bp_dot_deg <- process_and_plot(gsea_list_BP_deg, "DEG-Based GO Biological Process")
cc_dot_deg <- process_and_plot(gsea_list_CC_deg, "DEG-Based GO Cellular Component")
mf_dot_deg <- process_and_plot(gsea_list_MF_deg, "DEG-Based GO Molecular Function")
kegg_dot_deg <- process_and_plot(gsea_list_KEGG_deg, "DEG-Based KEGG Pathways")

# Create and save cnetplots
message("\nCreating DEG-based cnetplots...")

# BP cnetplots
message("\nProcessing DEG-based BP cnetplots...")
p_cnet_BP_orga_tissue_deg <- create_cnetplot(gsea_BP_orga_tissue_deg, "DEG-Based BP Orga vs Tissue", geneList_orga_vs_tissue_deg)
p_cnet_BP_pc_tissue_deg <- create_cnetplot(gsea_BP_pc_tissue_deg, "DEG-Based BP PC vs Tissue", geneList_pc_vs_tissue_deg)
p_cnet_BP_pc_orga_deg <- create_cnetplot(gsea_BP_pc_orga_deg, "DEG-Based BP PC vs Orga", geneList_pc_vs_orga_deg)

save_cnetplot(p_cnet_BP_orga_tissue_deg, "GSEA_DEG_cnetplot_BP_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_BP_pc_tissue_deg, "GSEA_DEG_cnetplot_BP_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_BP_pc_orga_deg, "GSEA_DEG_cnetplot_BP_PC_vs_Orga.jpeg")

# CC cnetplots
message("\nProcessing DEG-based CC cnetplots...")
p_cnet_CC_orga_tissue_deg <- create_cnetplot(gsea_CC_orga_tissue_deg, "DEG-Based CC Orga vs Tissue", geneList_orga_vs_tissue_deg)
p_cnet_CC_pc_tissue_deg <- create_cnetplot(gsea_CC_pc_tissue_deg, "DEG-Based CC PC vs Tissue", geneList_pc_vs_tissue_deg)
p_cnet_CC_pc_orga_deg <- create_cnetplot(gsea_CC_pc_orga_deg, "DEG-Based CC PC vs Orga", geneList_pc_vs_orga_deg)

save_cnetplot(p_cnet_CC_orga_tissue_deg, "GSEA_DEG_cnetplot_CC_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_CC_pc_tissue_deg, "GSEA_DEG_cnetplot_CC_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_CC_pc_orga_deg, "GSEA_DEG_cnetplot_CC_PC_vs_Orga.jpeg")

# MF cnetplots
message("\nProcessing DEG-based MF cnetplots...")
p_cnet_MF_orga_tissue_deg <- create_cnetplot(gsea_MF_orga_tissue_deg, "DEG-Based MF Orga vs Tissue", geneList_orga_vs_tissue_deg)
p_cnet_MF_pc_tissue_deg <- create_cnetplot(gsea_MF_pc_tissue_deg, "DEG-Based MF PC vs Tissue", geneList_pc_vs_tissue_deg)
p_cnet_MF_pc_orga_deg <- create_cnetplot(gsea_MF_pc_orga_deg, "DEG-Based MF PC vs Orga", geneList_pc_vs_orga_deg)

save_cnetplot(p_cnet_MF_orga_tissue_deg, "GSEA_DEG_cnetplot_MF_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_MF_pc_tissue_deg, "GSEA_DEG_cnetplot_MF_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_MF_pc_orga_deg, "GSEA_DEG_cnetplot_MF_PC_vs_Orga.jpeg")

# KEGG cnetplots
message("\nProcessing DEG-based KEGG cnetplots...")
p_cnet_KEGG_orga_tissue_deg <- create_cnetplot(gsea_KEGG_orga_tissue_deg, "DEG-Based KEGG Orga vs Tissue", geneList_orga_vs_tissue_deg)
p_cnet_KEGG_pc_tissue_deg <- create_cnetplot(gsea_KEGG_pc_tissue_deg, "DEG-Based KEGG PC vs Tissue", geneList_pc_vs_tissue_deg)
p_cnet_KEGG_pc_orga_deg <- create_cnetplot(gsea_KEGG_pc_orga_deg, "DEG-Based KEGG PC vs Orga", geneList_pc_vs_orga_deg)

save_cnetplot(p_cnet_KEGG_orga_tissue_deg, "GSEA_DEG_cnetplot_KEGG_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_KEGG_pc_tissue_deg, "GSEA_DEG_cnetplot_KEGG_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_KEGG_pc_orga_deg, "GSEA_DEG_cnetplot_KEGG_PC_vs_Orga.jpeg")

message("\nAll DEG-based analyses and visualizations complete!")


# ========================= #
# weighted correlation network analysis (WGCNA)
# ========================= #

# Load required libraries
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(igraph)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Create directory for WGCNA plots
dir.create("WGCNA_Results", showWarnings = FALSE)
dir.create("WGCNA_Results/Plots", showWarnings = FALSE)

# Step 1: Prepare expression data from DESeq2
# Get VST transformed data
vst_data <- assay(vst(dds, blind=FALSE))
print("Initial VST data dimensions:")
print(dim(vst_data))

# Transpose data (samples should be rows for WGCNA)
vst_data <- t(vst_data)
print("Transposed VST data dimensions:")
print(dim(vst_data))

# Load and prepare trait data
traits_data <- data.frame(
  Trait = factor(c(rep("Group1", 4), rep("Group2", 4), rep("Group3", 4)))
)
rownames(traits_data) <- c(paste0("PC_Insert", 1:4), 
                           paste0("Orga_Insert", 1:4),
                           paste0("Tissue", 1:4))

# Check sample matching and order
print("\nChecking sample name matching:")
print(all(rownames(traits_data) %in% rownames(vst_data)))
print(all(rownames(vst_data) %in% rownames(traits_data)))
vst_data <- vst_data[rownames(traits_data), ]

# Step 2: Pick soft-threshold power
powers <- c(1:30)  # Testing powers 1-30

# Pick soft threshold
sft <- pickSoftThreshold(
  vst_data,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

# Plot scale independence and mean connectivity
jpeg(file.path("WGCNA_Results", "Plots", "Scale_independence.jpg"),
     width = 15 * 300, height = 12 * 300, res = 300, quality = 100)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     main="Scale independence",
     type="n")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
abline(h=0.85,col="red")
dev.off()

par(mfrow = c(1,2))  # Add this again
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     main="Scale independence",
     type="n")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
abline(h=0.85,col="red")

# Plot mean connectivity
jpeg(file.path("WGCNA_Results", "Plots", "Mean_connectivity.jpg"),
     width = 10 * 300, height = 8 * 300, res = 300, quality = 100)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers, col="red")
dev.off()

# Find the best power
best_power <- min(sft$fitIndices$Power[sft$fitIndices$SFT.R.sq > 0.85])
print(paste("Best power:", best_power))

# Step 3: Create the network
# Calculate adjacency matrix
adjacency <- adjacency(vst_data,
                       power = best_power,
                       type = "signed")

# Convert adjacency matrix to TOM
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# Clustering using TOM
# Hierarchical clustering of genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module identification
minModuleSize <- 60
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

# Convert numeric labels to colors
dynamicColors <- labels2colors(dynamicMods)

# Plot dendrogram with colors
jpeg(file.path("WGCNA_Results", "Plots", "Dendrogram_with_colors.jpg"), 
     width = 10 * 300, height = 8 * 300, res = 300, quality = 100)
par(mfrow=c(1,1))
plotDendroAndColors(geneTree, 
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)
dev.off()

# Step 4: Module eigengenes analysis and merging
# Calculate initial eigengenes
MEList <- moduleEigengenes(vst_data, colors = dynamicColors)
MEs <- MEList$eigengenes

# Look at initial module correlations
print("\nInitial Module Eigengene correlations:")
print(cor(MEs))

# Calculate initial dissimilarity and clustering
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot initial module eigengene clustering
jpeg(file.path("WGCNA_Results", "Plots", "Initial_modules_clustering.jpg"),
     width = 10 * 300, height = 8 * 300, res = 300, quality = 100)
par(mfrow=c(1,1))
plot(METree, main = "Clustering of module eigengenes before merging",
     xlab = "", sub = "")
abline(h=0.25, col = "red")
dev.off()

# Plot initial module correlations as heatmap
jpeg(file.path("WGCNA_Results", "Plots", "Initial_modules_correlation.jpg"),
     width = 15 * 300, height = 12 * 300, res = 300, quality = 100)
heatmap(cor(MEs), main="Module correlations before merging")
dev.off()

# Merge close modules
merge <- mergeCloseModules(vst_data, 
                           dynamicColors, 
                           cutHeight = 0.25, 
                           verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# Set the nModules variable
nModules <- length(unique(mergedColors))

# Calculate and visualize dissimilarity after merging
MEDissThresh <- 1-cor(mergedMEs)
METreeThresh <- hclust(as.dist(MEDissThresh), method = "average")

# Plot dendrogram of merged modules
jpeg(file.path("WGCNA_Results", "Plots", "Merged_modules_clustering.jpg"),
     width = 10 * 300, height = 8 * 300, res = 300, quality = 100)
par(mfrow=c(1,1))
plot(METreeThresh, main = "Clustering of module eigengenes after merging",
     xlab = "", sub = "")
abline(h=0.25, col = "red")
dev.off()

# Plot correlation heatmap of merged modules
jpeg(file.path("WGCNA_Results", "Plots", "Merged_modules_correlation.jpg"),
     width = 15 * 300, height = 12 * 300, res = 300, quality = 100)
heatmap(cor(mergedMEs), main="Module correlations after merging")
dev.off()

# Plot dendrogram with merged colors
jpeg(file.path("WGCNA_Results", "Plots", "Merged_modules_dendrogram.jpg"), 
     width = 10 * 300, height = 8 * 300, res = 300, quality = 100)
par(mfrow=c(1,1))
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

# Step 5: Relate modules to traits and create optimized heatmap
# Convert categorical trait to dummy variables
traitGroups <- model.matrix(~ 0 + Trait, data = traits_data)
colnames(traitGroups) <- c("Group1", "Group2", "Group3")


# Start JPEG device 
jpeg(file.path("WGCNA_Results", "Plots", "Module_trait_relationships_with_legend.jpg"), 
     width = 18 * 300,  
     height = 12 * 300,   
     res = 300,           
     quality = 100)        

# Calculate correlations and p-values
moduleTraitCor <- cor(mergedMEs, traitGroups, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(vst_data))

# Create text matrix for heatmap values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", 
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Get module colors and labels
modColors <- substring(names(mergedMEs), 3)

# Set up the plot layout - one row, two columns
layout(matrix(c(1, 2), 1, 2), widths = c(4, 1.2))

# Set margins for heatmap
par(mar = c(5, 12, 4, 1))

# Create the heatmap
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = c("PC_Insert", "Orga_Insert", "Tissue"),
  yLabels = names(mergedMEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.4,       
  cex.lab = 1.5,              
  xLabelsAngle = 0,         
  zlim = c(-1, 1),            
  main = "Module-trait relationships"
)

# Set margins for module color legend panel
par(mar = c(5, 0, 4, 2))

# Create module color legend panel
plot(0, type = "n",
     xlab = "", ylab = "",
     xlim = c(0, 1),
     ylim = c(0, nModules),
     xaxt = "n", yaxt = "n",
     bty = "n")

# Add colored rectangles and module names
for(i in 1:nModules) {
  rect(0, i-1, 0.3, i,
       col = modColors[nModules-i+1],
       border = "black",
       lwd = 1)
  
  text(0.35, i-0.5,
       paste0(modColors[nModules-i+1], " (", table(mergedColors)[modColors[nModules-i+1]], ")"),
       adj = 0,
       cex = 1.3)
}

# Reset layout
layout(1)

# Close JPEG device
dev.off()

# Step 6: Export gene module assignments
moduleGenes <- list()
for(module in unique(mergedColors)) {
  moduleGenes[[module]] <- colnames(vst_data)[mergedColors == module]
}

# Print module sizes
print("\nModule sizes:")
for(module in names(moduleGenes)) {
  print(paste("Module", module, ":", length(moduleGenes[[module]]), "genes"))
}

# Step 7: Calculate Gene Significance and Module Membership
# Calculate gene significance for each trait
GS <- as.data.frame(cor(vst_data, traitGroups, use = "p"))
# Calculate module membership
MM <- as.data.frame(cor(vst_data, mergedMEs, use = "p"))

# Function to find hub genes for a single group
find_group_hub_genes <- function(moduleTraitCor, trait_column, GS, MM, 
                                 mergedColors, vst_data,
                                 cor_threshold = 0.6, 
                                 MM_threshold = 0.8, 
                                 GS_threshold = 0.5) {
  
  # Find modules with correlation > 0.7
  sig_modules <- which(moduleTraitCor[, trait_column] > cor_threshold)
  
  if(length(sig_modules) == 0) {
    message(paste("No modules found with correlation >", cor_threshold, "for", trait_column))
    return(NULL)
  }
  
  # Store results
  hub_genes_list <- list()
  
  # For each significant module
  for(mod_idx in sig_modules) {
    module_color <- substring(rownames(moduleTraitCor)[mod_idx], 3)
    message(paste("\nAnalyzing module:", module_color))
    
    # Get genes in this module
    module_genes <- mergedColors == module_color
    
    # Get MM and GS values for these genes
    module_MM <- abs(MM[module_genes, paste0("ME", module_color)])
    module_GS <- abs(GS[module_genes, trait_column])
    
    # Plot MM vs GS for this module with threshold lines
    jpeg(file.path("WGCNA_Results", "Plots", paste0("MM_GS_", trait_column, "_", module_color, ".jpg")),
         width = 8*300, height = 6*300, res = 300)
    par(mar=c(5,5,4,2))
    verboseScatterplot(module_MM,
                       module_GS,
                       xlab = paste("Module Membership in", module_color, "module"),
                       ylab = paste("Gene significance for", trait_column),
                       main = paste("Module membership vs. gene significance\n",
                                    module_color, "module -", trait_column),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
                       col = module_color)
    # Add threshold lines
    abline(h = GS_threshold, col = "red", lty = 2)
    abline(v = MM_threshold, col = "red", lty = 2)
    dev.off()
    
    # Find hub genes based on thresholds
    hub_genes <- which(module_MM >= MM_threshold & module_GS >= GS_threshold)
    
    if(length(hub_genes) > 0) {
      # Get gene names
      gene_names <- colnames(vst_data)[module_genes][hub_genes]
      
      # Store hub gene info
      hub_genes_list[[module_color]] <- data.frame(
        gene = gene_names,
        MM = module_MM[hub_genes],
        GS = module_GS[hub_genes],
        module = module_color,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all hub genes for this group
  if(length(hub_genes_list) > 0) {
    all_hub_genes <- do.call(rbind, hub_genes_list)
    # Save all hub genes for this group
    write.csv(all_hub_genes,
              file = paste0("hub_genes_", trait_column, ".csv"),
              row.names = FALSE)
    return(all_hub_genes)
  } else {
    return(NULL)
  }
}

# Function to compare hub genes between groups and find DEGs
compare_group_hub_genes <- function(group1_hubs, group2_hubs, dds, 
                                    group1_name, group2_name) {
  if(is.null(group1_hubs) || is.null(group2_hubs)) {
    message("One or both groups have no hub genes")
    return(NULL)
  }
  
  # Combine hub genes from both groups
  combined_genes <- unique(c(group1_hubs$gene, group2_hubs$gene))
  
  # Get DESeq2 results
  res_deg <- results(dds, contrast = c("condition", group1_name, group2_name))
  
  # Create data frame with DEG results for combined hub genes
  hub_degs <- data.frame(
    gene = combined_genes,
    log2FoldChange = res_deg[combined_genes, "log2FoldChange"],
    padj = res_deg[combined_genes, "padj"],
    stringsAsFactors = FALSE
  )
  
  # Add group information
  hub_degs$in_group1 <- hub_degs$gene %in% group1_hubs$gene
  hub_degs$in_group2 <- hub_degs$gene %in% group2_hubs$gene
  
  # Filter for significant DEGs
  sig_hub_degs <- hub_degs[!is.na(hub_degs$padj) & 
                             hub_degs$padj < 0.05 & 
                             abs(hub_degs$log2FoldChange) > 1, ]
  
  # Save results
  write.csv(hub_degs,
            file = paste0("hub_genes_", group1_name, "_vs_", group2_name, "_all.csv"),
            row.names = FALSE)
  
  if(nrow(sig_hub_degs) > 0) {
    write.csv(sig_hub_degs,
              file = paste0("hub_genes_", group1_name, "_vs_", group2_name, "_DEGs.csv"),
              row.names = FALSE)
  }
  
  return(list(all = hub_degs, deg = sig_hub_degs))
}

# Main analysis
# First find hub genes for each group
message("\nFinding hub genes for each group...")
group1_hubs <- find_group_hub_genes(
  moduleTraitCor = moduleTraitCor,
  trait_column = "Group1",
  GS = GS,
  MM = MM,
  mergedColors = mergedColors,
  vst_data = vst_data
)

group2_hubs <- find_group_hub_genes(
  moduleTraitCor = moduleTraitCor,
  trait_column = "Group2",
  GS = GS,
  MM = MM,
  mergedColors = mergedColors,
  vst_data = vst_data
)

group3_hubs <- find_group_hub_genes(
  moduleTraitCor = moduleTraitCor,
  trait_column = "Group3",
  GS = GS,
  MM = MM,
  mergedColors = mergedColors,
  vst_data = vst_data
)

# Then compare between groups
message("\nComparing hub genes between groups...")

# PC_Insert vs Orga_Insert
pc_vs_orga <- compare_group_hub_genes(
  group1_hubs = group1_hubs,
  group2_hubs = group2_hubs,
  dds = dds,
  group1_name = "PC_Insert",
  group2_name = "Orga_Insert"
)

# Orga_Insert vs Tissue
orga_vs_tissue <- compare_group_hub_genes(
  group1_hubs = group2_hubs,
  group2_hubs = group3_hubs,
  dds = dds,
  group1_name = "Orga_Insert",
  group2_name = "Tissue"
)

# PC_Insert vs Tissue
pc_vs_tissue <- compare_group_hub_genes(
  group1_hubs = group1_hubs,
  group2_hubs = group3_hubs,
  dds = dds,
  group1_name = "PC_Insert",
  group2_name = "Tissue"
)

# Print summary
message("\nSummary of results:")
message("\nHub genes in each group:")
message("Group1 (PC_Insert): ", ifelse(!is.null(group1_hubs), nrow(group1_hubs), 0))
message("Group2 (Orga_Insert): ", ifelse(!is.null(group2_hubs), nrow(group2_hubs), 0))
message("Group3 (Tissue): ", ifelse(!is.null(group3_hubs), nrow(group3_hubs), 0))

message("\nDifferentially expressed hub genes between groups:")
message("PC_Insert vs Orga_Insert: ", ifelse(!is.null(pc_vs_orga$deg), nrow(pc_vs_orga$deg), 0))
message("Orga_Insert vs Tissue: ", ifelse(!is.null(orga_vs_tissue$deg), nrow(orga_vs_tissue$deg), 0))
message("PC_Insert vs Tissue: ", ifelse(!is.null(pc_vs_tissue$deg), nrow(pc_vs_tissue$deg), 0))

print("Hub gene analysis complete!")

# ------------------------#
# Prepare Hub DEG-Based Gene Lists for GSEA Analysis
# ------------------------#

message("\n=== Preparing Hub DEG-Based Gene Lists for GSEA ===")

# Define a function to filter and prepare hub DEGs for GSEA
filter_hub_degs <- function(hub_deg_result, comparison_name) {
  # Ensure there are DEGs in the result
  if (is.null(hub_deg_result$deg) || nrow(hub_deg_result$deg) == 0) {
    message(sprintf("No significant hub DEGs found for %s", comparison_name))
    return(NULL)
  }
  
  # Sort the DEGs by log2 fold change in decreasing order for GSEA input
  geneList <- hub_deg_result$deg$log2FoldChange
  names(geneList) <- hub_deg_result$deg$gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Output detailed statistics
  message(sprintf("\nHub DEG Analysis for %s:", comparison_name))
  message(sprintf("Total number of significant hub DEGs: %d", length(geneList)))
  message(sprintf("Number of upregulated genes: %d", sum(geneList > 0)))
  message(sprintf("Number of downregulated genes: %d", sum(geneList < 0)))
  
  return(geneList)
}

# Apply filtering to each hub DEG comparison result
geneList_hub_degs_pc_vs_orga <- filter_hub_degs(pc_vs_orga, "PC vs Orga")
geneList_hub_degs_orga_vs_tissue <- filter_hub_degs(orga_vs_tissue, "Orga vs Tissue")
geneList_hub_degs_pc_vs_tissue <- filter_hub_degs(pc_vs_tissue, "PC vs Tissue")

# Check the structure and length of each hub DEG-based gene list
message("\nSummary of hub DEG lists:")
print(head(geneList_hub_degs_pc_vs_orga))
print(head(geneList_hub_degs_orga_vs_tissue))
print(head(geneList_hub_degs_pc_vs_tissue))

# Run GSEA analysis on hub DEGs
message("\n=== Running Hub DEG-Based GSEA Analyses ===")

# BP analyses for hub DEGs
message("\nRunning BP analyses for hub DEGs...")
gsea_BP_pc_tissue_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_tissue, "GO", "BP")
gsea_BP_orga_tissue_hub <- run_individual_gsea(geneList_hub_degs_orga_vs_tissue, "GO", "BP")
gsea_BP_pc_orga_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_orga, "GO", "BP")

# Save BP results for hub DEGs
save_pathway_details(gsea_BP_pc_tissue_hub, "BP_HubDEG", "PC_vs_Tissue")
save_pathway_details(gsea_BP_orga_tissue_hub, "BP_HubDEG", "Orga_vs_Tissue")
save_pathway_details(gsea_BP_pc_orga_hub, "BP_HubDEG", "PC_vs_Orga")

# CC analyses for hub DEGs
message("\nRunning CC analyses for hub DEGs...")
gsea_CC_pc_tissue_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_tissue, "GO", "CC")
gsea_CC_orga_tissue_hub <- run_individual_gsea(geneList_hub_degs_orga_vs_tissue, "GO", "CC")
gsea_CC_pc_orga_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_orga, "GO", "CC")

# Save CC results for hub DEGs
save_pathway_details(gsea_CC_pc_tissue_hub, "CC_HubDEG", "PC_vs_Tissue")
save_pathway_details(gsea_CC_orga_tissue_hub, "CC_HubDEG", "Orga_vs_Tissue")
save_pathway_details(gsea_CC_pc_orga_hub, "CC_HubDEG", "PC_vs_Orga")

# MF analyses for hub DEGs
message("\nRunning MF analyses for hub DEGs...")
gsea_MF_pc_tissue_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_tissue, "GO", "MF")
gsea_MF_orga_tissue_hub <- run_individual_gsea(geneList_hub_degs_orga_vs_tissue, "GO", "MF")
gsea_MF_pc_orga_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_orga, "GO", "MF")

# Save MF results for hub DEGs
save_pathway_details(gsea_MF_pc_tissue_hub, "MF_HubDEG", "PC_vs_Tissue")
save_pathway_details(gsea_MF_orga_tissue_hub, "MF_HubDEG", "Orga_vs_Tissue")
save_pathway_details(gsea_MF_pc_orga_hub, "MF_HubDEG", "PC_vs_Orga")

# KEGG analyses for hub DEGs
message("\nRunning KEGG analyses for hub DEGs...")
gsea_KEGG_pc_tissue_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_tissue, "KEGG")
gsea_KEGG_orga_tissue_hub <- run_individual_gsea(geneList_hub_degs_orga_vs_tissue, "KEGG")
gsea_KEGG_pc_orga_hub <- run_individual_gsea(geneList_hub_degs_pc_vs_orga, "KEGG")

# Save KEGG results for hub DEGs
save_pathway_details(gsea_KEGG_pc_tissue_hub, "KEGG_HubDEG", "PC_vs_Tissue")
save_pathway_details(gsea_KEGG_orga_tissue_hub, "KEGG_HubDEG", "Orga_vs_Tissue")
save_pathway_details(gsea_KEGG_pc_orga_hub, "KEGG_HubDEG", "PC_vs_Orga")

# Create lists for merging results for dotplots
gsea_list_BP_hub <- list(
  "Orga_vs_Tissue" = gsea_BP_orga_tissue_hub,
  "PC_vs_Tissue" = gsea_BP_pc_tissue_hub,
  "PC_vs_Orga" = gsea_BP_pc_orga_hub
)

gsea_list_CC_hub <- list(
  "Orga_vs_Tissue" = gsea_CC_orga_tissue_hub,
  "PC_vs_Tissue" = gsea_CC_pc_tissue_hub,
  "PC_vs_Orga" = gsea_CC_pc_orga_hub
)

gsea_list_MF_hub <- list(
  "Orga_vs_Tissue" = gsea_MF_orga_tissue_hub,
  "PC_vs_Tissue" = gsea_MF_pc_tissue_hub,
  "PC_vs_Orga" = gsea_MF_pc_orga_hub
)

gsea_list_KEGG_hub <- list(
  "Orga_vs_Tissue" = gsea_KEGG_orga_tissue_hub,
  "PC_vs_Tissue" = gsea_KEGG_pc_tissue_hub,
  "PC_vs_Orga" = gsea_KEGG_pc_orga_hub
)

# Create visualizations
message("\n=== Creating Hub DEG-Based Visualizations ===")

# Create dotplots
message("\nCreating hub DEG-based dotplots...")
bp_dot_hub <- process_and_plot(gsea_list_BP_hub, "Hub DEG GO Biological Process")
cc_dot_hub <- process_and_plot(gsea_list_CC_hub, "Hub DEG GO Cellular Component")
mf_dot_hub <- process_and_plot(gsea_list_MF_hub, "Hub DEG GO Molecular Function")
kegg_dot_hub <- process_and_plot(gsea_list_KEGG_hub, "Hub DEG KEGG Pathways")

# Create and save cnetplots
message("\nCreating hub DEG-based cnetplots...")

# BP cnetplots
message("\nProcessing hub DEG-based BP cnetplots...")
p_cnet_BP_orga_tissue_hub <- create_cnetplot(gsea_BP_orga_tissue_hub, "Hub DEG-Based BP Orga vs Tissue", geneList_hub_degs_orga_vs_tissue)
p_cnet_BP_pc_tissue_hub <- create_cnetplot(gsea_BP_pc_tissue_hub, "Hub DEG-Based BP PC vs Tissue", geneList_hub_degs_pc_vs_tissue)
p_cnet_BP_pc_orga_hub <- create_cnetplot(gsea_BP_pc_orga_hub, "Hub DEG-Based BP PC vs Orga", geneList_hub_degs_pc_vs_orga)

save_cnetplot(p_cnet_BP_orga_tissue_hub, "GSEA_HubDEG_cnetplot_BP_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_BP_pc_tissue_hub, "GSEA_HubDEG_cnetplot_BP_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_BP_pc_orga_hub, "GSEA_HubDEG_cnetplot_BP_PC_vs_Orga.jpeg")

# CC cnetplots
message("\nProcessing hub DEG-based CC cnetplots...")
p_cnet_CC_orga_tissue_hub <- create_cnetplot(gsea_CC_orga_tissue_hub, "Hub DEG-Based CC Orga vs Tissue", geneList_hub_degs_orga_vs_tissue)
p_cnet_CC_pc_tissue_hub <- create_cnetplot(gsea_CC_pc_tissue_hub, "Hub DEG-Based CC PC vs Tissue", geneList_hub_degs_pc_vs_tissue)
p_cnet_CC_pc_orga_hub <- create_cnetplot(gsea_CC_pc_orga_hub, "Hub DEG-Based CC PC vs Orga", geneList_hub_degs_pc_vs_orga)

save_cnetplot(p_cnet_CC_orga_tissue_hub, "GSEA_HubDEG_cnetplot_CC_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_CC_pc_tissue_hub, "GSEA_HubDEG_cnetplot_CC_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_CC_pc_orga_hub, "GSEA_HubDEG_cnetplot_CC_PC_vs_Orga.jpeg")

# MF cnetplots
message("\nProcessing hub DEG-based MF cnetplots...")
p_cnet_MF_orga_tissue_hub <- create_cnetplot(gsea_MF_orga_tissue_hub, "Hub DEG-Based MF Orga vs Tissue", geneList_hub_degs_orga_vs_tissue)
p_cnet_MF_pc_tissue_hub <- create_cnetplot(gsea_MF_pc_tissue_hub, "Hub DEG-Based MF PC vs Tissue", geneList_hub_degs_pc_vs_tissue)
p_cnet_MF_pc_orga_hub <- create_cnetplot(gsea_MF_pc_orga_hub, "Hub DEG-Based MF PC vs Orga", geneList_hub_degs_pc_vs_orga)

save_cnetplot(p_cnet_MF_orga_tissue_hub, "GSEA_HubDEG_cnetplot_MF_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_MF_pc_tissue_hub, "GSEA_HubDEG_cnetplot_MF_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_MF_pc_orga_hub, "GSEA_HubDEG_cnetplot_MF_PC_vs_Orga.jpeg")

# KEGG cnetplots
message("\nProcessing hub DEG-based KEGG cnetplots...")
p_cnet_KEGG_orga_tissue_hub <- create_cnetplot(gsea_KEGG_orga_tissue_hub, "Hub DEG-Based KEGG Orga vs Tissue", geneList_hub_degs_orga_vs_tissue)
p_cnet_KEGG_pc_tissue_hub <- create_cnetplot(gsea_KEGG_pc_tissue_hub, "Hub DEG-Based KEGG PC vs Tissue", geneList_hub_degs_pc_vs_tissue)
p_cnet_KEGG_pc_orga_hub <- create_cnetplot(gsea_KEGG_pc_orga_hub, "Hub DEG-Based KEGG PC vs Orga", geneList_hub_degs_pc_vs_orga)

save_cnetplot(p_cnet_KEGG_orga_tissue_hub, "GSEA_HubDEG_cnetplot_KEGG_Orga_vs_Tissue.jpeg")
save_cnetplot(p_cnet_KEGG_pc_tissue_hub, "GSEA_HubDEG_cnetplot_KEGG_PC_vs_Tissue.jpeg")
save_cnetplot(p_cnet_KEGG_pc_orga_hub, "GSEA_HubDEG_cnetplot_KEGG_PC_vs_Orga.jpeg")

message("\nAll hub DEG-based analyses and visualizations complete!")

# ========================= #
# PCA and Volcano Plot Analysis 
# ========================= #

# Create directory for plots 
dir.create("Visualization_Results", showWarnings = FALSE)
dir.create("Visualization_Results/PCA", showWarnings = FALSE)
dir.create("Visualization_Results/Volcano", showWarnings = FALSE)


# Volcano Plot function 
create_volcano_plot <- function(dds, condition1, condition2, genes_to_label = NULL) {
  tryCatch({
    res <- results(dds, contrast=c("condition", condition1, condition2))
    res <- na.omit(res)
    
    volcano_data <- as.data.frame(res)
    volcano_data$gene_name <- gene_names_fpkm_named[rownames(volcano_data)]
    
    volcano_data$significance <- case_when(
      volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1 ~ 
        ifelse(volcano_data$log2FoldChange > 1, "Up-regulated", "Down-regulated"),
      TRUE ~ "Not significant"
    )
    
    gene_counts <- table(volcano_data$significance)
    
    # Set conditional positions
    y_pos <- if(condition1 == "PC_Insert" && condition2 == "Orga_Insert") {
      110  # Lower position for PC vs Orga fold change labels
    } else {
      210  # Original position for other comparisons
    }
    
    x_pos_padj <- if(condition1 == "PC_Insert" && condition2 == "Orga_Insert") {
      8.5   # Adjusted x position for PC vs Orga p.adj text
    } else {
      10   # Original position for other comparisons
    }
    
    volcano_plot <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
      geom_point(alpha=0.7, size=2) +
      scale_color_manual(
        values=c("Up-regulated"="red", "Down-regulated"="blue", "Not significant"="grey"),
        labels=paste0(names(gene_counts), " (n=", gene_counts, ")")
      ) +
      geom_vline(xintercept=c(-1, 1), linetype="dashed", color="darkgrey", linewidth=1.5) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey", linewidth=1.5) +
      
      annotate("text", x=-1.5, y=-5, 
               label="-1", size=8, fontface="bold") +
      annotate("text", x=1.3, y=-5, 
               label="1", size=8, fontface="bold") +
      
      # Fold change labels
      annotate("text", x=-1.5, y=y_pos, 
               label=expression(bold("Fold Change = " ~ bold(frac(1, 2)))), 
               size=7, angle=90, parse=TRUE) +
      annotate("text", x=1.5, y=y_pos, 
               label=expression(bold("Fold Change = " ~ bold(2))), 
               size=7, angle=90, parse=TRUE) +
      
      annotate("text", x=-20, y=-log10(0.05)+ 4, 
               label="1.3", size=8, fontface="bold") +
      # Use conditional x_pos_padj for p.adj text
      annotate("text", x=x_pos_padj, y=-log10(0.05) + 4, 
               label="P.adj value = 0.05", size=8, fontface="bold") +
      
      labs(
        title=paste(condition1, "vs", condition2),
        x=expression(log[2]~"Fold Change"),
        y=expression(-log[10]~"Adjusted P-value"),
        color="Differential Expression"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size=28, face="bold", hjust=0.5),
        axis.title = element_text(size=24, face="bold"),
        axis.text = element_text(size=20),
        legend.title = element_text(size=24, face="bold"),
        legend.text = element_text(size=20)
      ) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    
    if (!is.null(genes_to_label)) {
      matched_genes <- volcano_data$gene_name[toupper(volcano_data$gene_name) %in% toupper(genes_to_label)]
      genes_to_show <- volcano_data[volcano_data$gene_name %in% matched_genes, ]
      
      volcano_plot <- volcano_plot +
        geom_text_repel(
          data = genes_to_show,
          aes(label = gene_name),
          size = 5,
          box.padding = 1,
          point.padding = 0.5,
          force = 10,
          max.overlaps = Inf
        )
    }
    
    ggsave(file.path("Visualization_Results", "Volcano",  
                     paste0("Volcano_", condition1, "_vs_", condition2, ".jpg")),
           volcano_plot, width = 12, height = 8, dpi = 300)
    return(volcano_plot)
  }, error = function(e) {
    message("Error in volcano plot creation: ", e$message)
    return(NULL)
  })
}

# PCA Plot function 
create_pca_plot <- function(dds) {
  tryCatch({
    vsd <- vst(dds, blind=FALSE)
    pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    pca_plot <- ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, label=name)) +
      geom_point(size=6, alpha=0.8) +
      geom_text_repel(
        size=6,
        box.padding = 1.5,
        point.padding = 0.7,
        force = 12,
        max.overlaps = Inf
      ) +
      scale_color_brewer(palette="Set1") +
      scale_x_continuous(
        breaks = scales::pretty_breaks(n = 12),
        limits = function(x) c(floor(min(x)), ceiling(max(x)))
      ) +
      scale_y_continuous(
        breaks = scales::pretty_breaks(n = 12),
        limits = function(x) c(floor(min(x)), ceiling(max(x)))
      ) +
      labs(
        title="Principal Component Analysis",
        subtitle=paste0("Total Variance Explained: ", sum(percentVar[1:2]), "%"),
        x=paste0("PC1 (", percentVar[1], "% variance)"),
        y=paste0("PC2 (", percentVar[2], "% variance)"),
        color="Condition"
      ) +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      theme_minimal() +
      theme(
        plot.title = element_text(size=24, face="bold", hjust=0.5),
        plot.subtitle = element_text(size=18, hjust=0.5),
        axis.title = element_text(size=20, face="bold"),
        axis.text = element_text(size=16),
        legend.title = element_text(size=18, face="bold"),
        legend.text = element_text(size=16),
        legend.position = "right",
        panel.grid.major = element_line(color="grey90"),
        panel.grid.minor = element_line(color="grey95")
      )
    
    ggsave(file.path("Visualization_Results", "PCA", "PCA_plot.jpg"),
           pca_plot, width = 12, height = 8, dpi = 300)
    return(pca_plot)
  }, error = function(e) {
    message("Error in PCA plot creation: ", e$message)
    return(NULL)
  })
}

message("\n=== Creating Visualizations PCA and Volcano Plots===")

pca_plot <- create_pca_plot(dds)

comparisons <- list(
  c("PC_Insert", "Tissue"),
  c("Orga_Insert", "Tissue"),
  c("PC_Insert", "Orga_Insert")
)

for(comp in comparisons) {
  volcano_plot <- create_volcano_plot(dds, comp[1], comp[2])
}

# ========================= #
# Heatmap Analysis 
# ========================= #

# Create directories

dir.create("Visualization_Results/Heatmap", showWarnings = FALSE)
dir.create("Visualization_Results/Heatmap/extra", showWarnings = FALSE)
# Log2 FPKM version
create_log2fpkm_heatmap <- function(genes_of_interest, filename_prefix) {
  tryCatch({
    matched_ids <- names(gene_names_fpkm_named)[toupper(gene_names_fpkm_named) %in% toupper(genes_of_interest)]
    
    fpkm_with_pseudo <- fpkm_final_filtered[matched_ids, 1:(ncol(fpkm_final_filtered)-1)] + 1
    expr_data <- as.matrix(log2(fpkm_with_pseudo))
    rownames(expr_data) <- gene_names_fpkm_named[matched_ids]
    
    conditions <- factor(c(rep("PC_Insert", 4), rep("Orga_Insert", 4), rep("Tissue", 4)))
    
    ha <- HeatmapAnnotation(
      condition = conditions,
      col = list(condition = c(
        "PC_Insert" = "#E41A1C",
        "Orga_Insert" = "#4DAF4A",
        "Tissue" = "#377EB8"
      )),
      annotation_name_side = "left",
      # The word "condition" above heatmap
      annotation_name_gp = gpar(fontsize = 22),
      # Settings for condition legend (top-right)
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 22),    # "condition" in legend
        labels_gp = gpar(fontsize = 20)     # PC_Insert, Tissue, etc in legend
      )
    )
    
    ht <- Heatmap(
      expr_data,
      name = "log2(FPKM)",
      col = colorRamp2(c(min(expr_data), mean(expr_data), max(expr_data)), c("blue", "white", "red")),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      # Gene names on right side
      row_names_gp = gpar(fontsize = 22),
      # Sample names at bottom
      column_names_gp = gpar(fontsize = 22),
      top_annotation = ha,
      height = unit(20, "cm"),
      width = unit(25, "cm"),
      # Settings for log2(FPKM) legend
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 20),     # "log2(FPKM)" text
        labels_gp = gpar(fontsize = 20),     # Numbers in legend
        legend_height = unit(8, "cm"),
        grid_width = unit(2, "cm")
      )
    )
    
    jpeg(file.path("Visualization_Results", "Heatmap","extra", paste0(filename_prefix, "_log2fpkm_heatmap.jpg")), 
         width = 5000, height = 3600, res = 300)
    draw(ht, padding = unit(c(3, 12, 3, 12), "cm"))
    dev.off()
  })
}

# VST version
create_vst_heatmap <- function(dds, genes_of_interest, filename_prefix) {
  tryCatch({
    vsd <- vst(dds, blind=FALSE)
    expr_data <- assay(vsd)
    
    matched_ids <- names(gene_names_fpkm_named)[toupper(gene_names_fpkm_named) %in% toupper(genes_of_interest)]
    expr_data <- expr_data[matched_ids, ]
    rownames(expr_data) <- gene_names_fpkm_named[matched_ids]
    
    ha <- HeatmapAnnotation(
      condition = colData(dds)$condition,
      col = list(condition = c(
        "PC_Insert" = "#E41A1C",
        "Orga_Insert" = "#4DAF4A",
        "Tissue" = "#377EB8"
      )),
      annotation_name_side = "left",
      # The word "condition" above heatmap
      annotation_name_gp = gpar(fontsize = 22),
      # Settings for condition legend (top-right)
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 22),    # "condition" in legend
        labels_gp = gpar(fontsize = 20)     # PC_Insert, Tissue, etc in legend
      )
    )
    
    ht <- Heatmap(
      expr_data,
      name = "VST",
      col = colorRamp2(c(min(expr_data), mean(expr_data), max(expr_data)), c("blue", "white", "red")),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      # Gene names on right side
      row_names_gp = gpar(fontsize = 22),
      # Sample names at bottom
      column_names_gp = gpar(fontsize = 22),
      top_annotation = ha,
      height = unit(20, "cm"),
      width = unit(25, "cm"),
      # Settings for VST legend
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 20),     # "VST" text
        labels_gp = gpar(fontsize = 20),     # Numbers in legend
        legend_height = unit(8, "cm"),
        grid_width = unit(2, "cm")
      )
    )
    
    jpeg(file.path("Visualization_Results", "Heatmap","extra", paste0(filename_prefix, "_vst_heatmap.jpg")), 
         width = 5000, height = 3600, res = 300)
    draw(ht, padding = unit(c(3, 12, 3, 12), "cm"))
    dev.off()
  })
}

# Define gene lists
gene_lists <- list(
  Stratum-Basale = c("IGFBP5", "COL17A1", "SYNE2", "FABP5", "UHRF1", "HELLS", "THEM5", "PDLIM4", "PNISR", "KRT5", "KRT14", "SELENBP1"),
  Stratum-Spinosum = c("S100A8", "S100A9", "CA1", "CA12", "ACAA2", "BDH1", "ACAT1", "SLC16A1", "TST", "KRT6A", "KRT15", "KRT17", "GSTA1", "LYPD3", "YBX1", "DHRS7", "RRD-SPRR"),
  Differentiating-Keratinocytes = c("DSG2", "DSG3", "DSC2", "DSC3", "KRT4", "KRT13", "KRT15", "KRT17", "KRT36", "FABP5", "KRTDAP", "SPINK5", "CSTA", "GSTA1", "IGF2", "PI3", "TFF3", "RBP2", "MGST1"),
  Terminally-differentiated-Keratinocyes = c("DSG1", "DSC1", "DSC3", "KRT23", "EVPL", "PPL", "PKP1", "IVL", "ANXA1", "RPS2", "RPL8", "AHNAK", "DSP", "JUP", "S100A9", "CHMP2A", "DSTN"),
  Fibroblasts = c("COL3A1", "COL1A1", "POSTN", "DCN", "MFAP5", "APOE", "KRT8"),
  Immune-related = c("C1QC","C1QA","C1QB","CD4","CD69", "ISG15","MX1", "MX2", "RSAD2", "IFIT3","IFIT5"),
  Macrophages = c("C1QC","C1QA","C1QB"),
  Claudins = c("CLDN1", "CLDN4", "CLDN3", "CLDN5", "CLDN6", "CLDN7", "CLDN12", "CLDN17", "CLDN23"),
  TRP-channels = c("TRPM4", "TRPV4", "TRPV3", "TRPM6", "TRPM7", "TRPV2", "TRPV6", "TRPC1"),
  Gap-junctions = c("GJA1", "GJA9", "GJB2", "GJB3", "GJB5", "GJB6", "GJC2"),
  Endothelials = c("PECAM1","CDH5","CLDN5")
  
)

# Run both log2FPKM and VST versions for the gene lists
for(name in names(gene_lists)) {
  create_log2fpkm_heatmap(gene_lists[[name]], name)
  create_vst_heatmap(dds, gene_lists[[name]], name)
}

############    End of the Code  #############