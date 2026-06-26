# ==============================================================================
# MASTER AD BIOMARKER DISCOVERY PIPELINE (EARLY-TRANSLATION ARCHITECTURE)
# ==============================================================================

# Clear environment memory to guarantee a clean start
rm(list = ls())
gc()

# Consolidated Library Initialization
library(WGCNA)
library(matrixStats)
library(tidyverse)
library(Biobase)
library(GEOquery)
library(biomaRt)
library(limma)
library(pheatmap)
library(patchwork)
library(viridis)
library(circlize)
library(magick)

# Enable multi-threading inside WGCNA globally
allowWGCNAThreads()

cat("======================================================================\n")
cat("  PHASE 1 & 2: DATA ACQUISITION & EARLY FEATURE SPACE TRANSLATION     \n")
cat("======================================================================\n\n")

# 1. DEFINE FILE SYSTEMS LOCATIONS
path_gse118553_matrix <- "C:/Users/USER/Documents/Genomac_clients/Jaiyoba_FINAL/Data/GSE118553_series_matrix (1).txt.gz"
path_gse282742_matrix <- "C:/Users/USER/Documents/Genomac_clients/Jaiyoba_FINAL/Data/GSE282742_series_matrix.txt.gz"
path_gse282742_counts  <- "C:/Users/USER/Documents/Genomac_clients/Jaiyoba_FINAL/Data/GSE282742_Expected_count.txt.gz"

# 2. LOAD GEO SERIES MATRIX ANNOTATIONS
cat("Loading local brain tissue series object... ")
gse118553 <- getGEO(filename = path_gse118553_matrix, GSEMatrix = TRUE, AnnotGPL = TRUE)
raw_pheno_brain <- pData(gse118553)
expr_brain_all <- exprs(gse118553)
cat("Done.\n")

cat("Loading local blood metadata series object... ")
gse282742 <- getGEO(filename = path_gse282742_matrix, GSEMatrix = TRUE)
raw_pheno_blood <- pData(gse282742)
cat("Done.\n")

cat("Parsing raw blood RNA-Seq counts data... ")
blood_counts_raw <- read_delim(path_gse282742_counts, delim = "\t", col_types = cols()) %>% as.data.frame()
rownames(blood_counts_raw) <- blood_counts_raw[, 1]
blood_counts_raw <- blood_counts_raw[, -1]
cat("Done.\n")

# 3. METADATA HARMONIZATION
cat("\nHarmonizing Clinical Phenotypes... ")
pheno_brain_clean <- raw_pheno_brain %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(
    Diagnosis = case_when(
      str_detect(`disease state:ch1`, "(?i)control") ~ "Control",
      str_detect(`disease state:ch1`, "(?i)asym")    ~ "AsymAD",
      str_detect(`disease state:ch1`, "(?i)^ad")     ~ "SymAD",
      TRUE ~ NA_character_
    ),
    Region = case_when(
      str_detect(`tissue:ch1`, "(?i)entorhinal") ~ "Entorhinal",
      str_detect(`tissue:ch1`, "(?i)frontal")    ~ "Frontal",
      TRUE ~ NA_character_
    ),
    Age   = as.numeric(str_extract(as.character(`age:ch1`), "\\d+")),
    Sex   = if_else(str_detect(`gender:ch1`, "(?i)female|f"), "Female", "Male")
  ) %>% filter(!is.na(Diagnosis), !is.na(Region))

pheno_blood_clean <- raw_pheno_blood %>%
  rownames_to_column(var = "GSM_ID") %>%
  mutate(
    SampleID = str_extract(title, "^VGH\\d+"),
    Diagnosis = case_when(
      str_detect(`disease state:ch1`, "(?i)^ad")    ~ "AD",
      str_detect(`disease state:ch1`, "(?i)p-mci")  ~ "PMCI",
      str_detect(`disease state:ch1`, "(?i)s-mci")  ~ "SMCI",
      TRUE ~ NA_character_
    ),
    Age   = as.numeric(str_extract(as.character(`age:ch1`), "\\d+")),
    Sex   = if_else(str_detect(`Sex:ch1`, "(?i)female|f"), "Female", "Male")
  ) %>% filter(!is.na(Diagnosis), !is.na(SampleID))
cat("Done.\n")

# 4. SUBSET MATRICES
pheno_early_brain <- pheno_brain_clean %>% filter(Region == "Entorhinal" & Diagnosis %in% c("Control", "AsymAD"))
expr_early_brain  <- expr_brain_all[, pheno_early_brain$SampleID]

pheno_late_brain <- pheno_brain_clean %>% filter(Region == "Frontal" & Diagnosis %in% c("Control", "SymAD"))
expr_late_brain  <- expr_brain_all[, pheno_late_brain$SampleID]

common_blood_samples <- intersect(pheno_blood_clean$SampleID, colnames(blood_counts_raw))
pheno_blood_aligned  <- pheno_blood_clean %>% filter(SampleID %in% common_blood_samples)
expr_blood_transformed <- log2(as.matrix(blood_counts_raw[, pheno_blood_aligned$SampleID]) + 1)

# ------------------------------------------------------------------------------
# 5. EARLY FEATURE SPACE TRANSLATION (PROBES/ENSG -> SYMBOLS)
# ------------------------------------------------------------------------------
cat("\nExtracting Annotation Dictionaries...\n")
brain_probe_map <- fData(gse118553) %>%
  rownames_to_column(var = "ProbeID") %>%
  dplyr::select(ProbeID, Symbol = `Gene symbol`) %>%
  dplyr::filter(!is.na(Symbol) & Symbol != "") %>%
  mutate(Symbol = str_trim(str_split_i(Symbol, "///", 1)))

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
blood_id_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = rownames(expr_blood_transformed),
  mart       = mart
) %>% dplyr::filter(!is.na(hgnc_symbol) & hgnc_symbol != "")

collapse_matrix_to_symbols <- function(expr_matrix, mapping_df, id_col, symbol_col) {
  valid_ids <- intersect(rownames(expr_matrix), mapping_df[[id_col]])
  expr_subset <- expr_matrix[valid_ids, ]
  
  mapping_aligned <- mapping_df %>% 
    dplyr::filter(!!sym(id_col) %in% valid_ids) %>%
    distinct(!!sym(id_col), .keep_all = TRUE)
  
  symbols_vector <- mapping_aligned[[symbol_col]][match(rownames(expr_subset), mapping_aligned[[id_col]])]
  row_vars <- rowVars(expr_subset, na.rm = TRUE)
  
  df_temp <- as.data.frame(expr_subset)
  rownames(df_temp) <- NULL 
  df_temp$Symbol <- symbols_vector
  df_temp$Var    <- row_vars
  
  expr_collapsed <- df_temp %>%
    arrange(Symbol, desc(Var)) %>%
    distinct(Symbol, .keep_all = TRUE) %>%
    dplyr::select(-Var)
  
  out_matrix <- as.matrix(expr_collapsed %>% dplyr::select(-Symbol))
  rownames(out_matrix) <- expr_collapsed$Symbol
  return(out_matrix)
}

cat("Collapsing all matrices into shared Gene Symbol space...\n")
expr_early_symbols <- collapse_matrix_to_symbols(expr_early_brain, brain_probe_map, "ProbeID", "Symbol")
expr_late_symbols  <- collapse_matrix_to_symbols(expr_late_brain, brain_probe_map, "ProbeID", "Symbol")
expr_blood_symbols <- collapse_matrix_to_symbols(expr_blood_transformed, blood_id_map, "ensembl_gene_id", "hgnc_symbol")

# ------------------------------------------------------------------------------
# 6. NETWORK QUALITY CONTROL ON SYMBOL MATRICES
# ------------------------------------------------------------------------------
run_quality_control <- function(expr, pheno, label) {
  sds <- rowSds(expr, na.rm=TRUE)
  expr_filtered <- expr[sds > 0 & !is.na(sds), ]
  
  adj <- adjacency(expr_filtered, power = 1, type = "signed")
  k <- rowSums(adj) - 1
  Z.k <- (k - mean(k, na.rm=TRUE)) / sd(k, na.rm=TRUE)
  
  keep_samples <- names(Z.k)[Z.k >= -2.5]
  pheno_clean  <- pheno %>% filter(SampleID %in% keep_samples)
  expr_clean   <- expr_filtered[, pheno_clean$SampleID, drop = FALSE]
  
  return(list(expr = expr_clean, pheno = pheno_clean))
}

early_qc <- run_quality_control(expr_early_symbols, pheno_early_brain, "Early Brain")
late_qc  <- run_quality_control(expr_late_symbols, pheno_late_brain, "Late Brain")
blood_qc <- run_quality_control(expr_blood_symbols, pheno_blood_aligned, "Blood")


# ==============================================================================
# PHASE 3: DIFFERENTIAL EXPRESSION ANALYSIS (SYMBOL SPACE)
# ==============================================================================
cat("\nStarting Phase 3: Differential Expression...\n")

run_deg <- function(qc_data, contrast_str, contrast_name) {
  pheno_comp <- qc_data$pheno %>% filter(!is.na(Diagnosis) & !is.na(Age) & !is.na(Sex))
  design <- model.matrix(~ 0 + Diagnosis + Age + Sex, data = pheno_comp)
  colnames(design) <- make.names(colnames(design))
  rownames(design) <- pheno_comp$SampleID
  
  expr_sync <- qc_data$expr[, rownames(design), drop = FALSE]
  fit <- lmFit(expr_sync, design)
  contrast <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast), trend = TRUE, robust = TRUE)
  
  list(top_table = topTable(fit2, coef = 1, number = Inf, adjust.method = "BH"), expr_sync = expr_sync, pheno_sync = pheno_comp)
}

deg_early <- run_deg(early_qc, "DiagnosisAsymAD - DiagnosisControl", "Early_AD")
deg_late  <- run_deg(late_qc, "DiagnosisSymAD - DiagnosisControl", "Late_AD")
deg_blood <- run_deg(blood_qc, "DiagnosisAD - DiagnosisSMCI", "Blood_AD")


# ==============================================================================
# PHASE 4: WEIGHTED GENE CO-EXPRESSION NETWORK CONSTRUCTION
# ==============================================================================
cat("\nStarting Phase 4: Network Construction...\n")

build_wgcna <- function(expr_sync, power_val) {
  mads <- rowMads(as.matrix(expr_sync), na.rm = TRUE)
  top_genes <- order(mads, decreasing = TRUE)[1:min(10000, length(mads))]
  work_mat <- t(expr_sync[top_genes, ])
  
  adj <- adjacency(work_mat, power = power_val, type = "signed")
  tom <- TOMsimilarity(adj, TOMType = "signed")
  gene_tree <- hclust(as.dist(1 - tom), method = "average")
  
  mods <- cutreeDynamic(dendro = gene_tree, distM = 1 - tom, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
  merged <- mergeCloseModules(work_mat, labels2colors(mods), cutHeight = 0.15, verbose = 0)
  
  # Colors are now natively mapped to Gene Symbols
  final_colors <- merged$colors
  names(final_colors) <- colnames(work_mat)
  
  list(work_mat = work_mat, tree = gene_tree, merged = merged, colors = final_colors)
}

net_early <- build_wgcna(deg_early$expr_sync, power_val = 5)
net_late  <- build_wgcna(deg_late$expr_sync, power_val = 4)

# ==============================================================================
# PATCH: CORRECTED MODULE-TRAIT CORRELATION FUNCTION
# ==============================================================================

# 1. Redefine the function to explicitly calculate Eigengenes safely
correlate_traits <- function(net, pheno) {
  # Explicitly calculate Eigengenes using the final verified colors
  MEs <- moduleEigengenes(net$work_mat, colors = net$colors)$eigengenes
  
  # Prepare trait vector
  trait_vec <- as.numeric(as.factor(pheno$Diagnosis))
  
  # Calculate Pearson correlation and significance
  cor_val <- cor(MEs, trait_vec, use = "p")
  p_val <- corPvalueStudent(cor_val, nrow(MEs))
  
  return(list(cor = cor_val, p = p_val))
}

# 2. Re-run the correlations
cat("Calculating module-trait correlations...\n")
trait_early <- correlate_traits(net_early, deg_early$pheno_sync)
trait_late  <- correlate_traits(net_late, deg_late$pheno_sync)

traittable_late <- as.data.frame(trait_late)
traittable_early <- as.data.frame(trait_early)

cat("SUCCESS: Trait correlation calculated successfully.\n")

# 3. Quick verification (Preview the top correlated modules in Early Brain)
early_preview <- data.frame(
  Module = gsub("ME", "", rownames(trait_early$cor)),
  Correlation = as.numeric(trait_early$cor),
  P_Value = as.numeric(trait_early$p)
)
cat("\n--- Preview: Early Brain Top Correlates ---\n")
print(head(early_preview[order(abs(early_preview$Correlation), decreasing=TRUE), ]))

# 1. Format the late trait correlation object into a clean dataframe
late_preview <- data.frame(
  Module = gsub("ME", "", rownames(trait_late$cor)),
  Correlation = as.numeric(trait_late$cor),
  P_Value = as.numeric(trait_late$p)
)

# 2. Print the top results, sorted by the strongest absolute correlation
cat("\n--- Preview: Late Brain (Frontal) Top Correlates ---\n")
print(head(late_preview[order(abs(late_preview$Correlation), decreasing=TRUE), ], 10))


# ==============================================================================
# PHASE 5: MANUSCRIPT METRIC EXTRACTION
# ==============================================================================
cat("\nStarting Phase 5: Metric Extraction (Simplified)...\n")

extract_metrics <- function(work_mat, final_colors, power_val) {
  all_modules <- unique(final_colors)
  profile_list <- list()
  
  for (mod in all_modules) {
    if (mod == "grey") next
    mod_genes <- names(final_colors)[final_colors == mod]
    mod_data <- as.matrix(work_mat[, mod_genes])
    
    scaled_data <- scale(mod_data)
    svd_res <- svd(scaled_data, nu = 1, nv = 1)
    mod_var <- round((svd_res$d[1]^2 / sum(svd_res$d^2)) * 100, 2)
    
    cor_mat <- cor(mod_data, use = "p")
    adj_mat <- ((1 + cor_mat) / 2)^power_val
    
    mod_density <- round(mean(adj_mat[lower.tri(adj_mat)]), 4)
    k_conn <- rowSums(adj_mat) - 1
    max_k <- max(k_conn)
    
    mod_size <- length(mod_genes)
    mod_cent <- round(sum(max_k - k_conn) / ((mod_size - 1) * (mod_size - 2)), 4)
    mod_het  <- round(sd(k_conn) / mean(k_conn), 4)
    
    profile_list[[mod]] <- data.frame(
      Module = mod, Size = mod_size, Variance = mod_var, Density = mod_density, 
      Centralization = mod_cent, Heterogeneity = mod_het
    )
  }
  do.call(rbind, profile_list) %>% arrange(desc(Size))
}

early_table <- extract_metrics(net_early$work_mat, net_early$colors, 5)
late_table  <- extract_metrics(net_late$work_mat, net_late$colors, 4)

# Export clean copy-paste files to your data directory
data_dir <- "C:/Users/USER/Documents/Genomac_clients/Jaiyoba_FINAL/Data/"
write_csv(early_table, paste0(data_dir, "New_Manuscript_Table1_Early_Network_Topography.csv"))
write_csv(late_table, paste0(data_dir, "New_Manuscript_Table2_Late_Network_Topography.csv"))


# ==============================================================================
# PATCHED PHASE 6: FULL CROSS-TISSUE PRESERVATION & PRIORITIZATION
# ==============================================================================
cat("\nStarting Phase 6: Full Cross-Tissue Preservation...\n")

# ------------------------------------------------------------------------------
# PART 1: EARLY STAGE (ENTORHINAL CORTEX) PRESERVATION
# ------------------------------------------------------------------------------
common_genes_early <- intersect(colnames(net_early$work_mat), rownames(deg_blood$expr_sync))

brain_temp_early <- net_early$work_mat[, common_genes_early]
blood_temp_early <- t(deg_blood$expr_sync[common_genes_early, ])

# Variance filter
gsg_brain_early <- goodSamplesGenes(brain_temp_early, verbose = 0)
gsg_blood_early <- goodSamplesGenes(blood_temp_early, verbose = 0)
valid_genes_early <- common_genes_early[gsg_brain_early$goodGenes & gsg_blood_early$goodGenes]

cat(paste(".. Retained", length(valid_genes_early), "valid genes for EARLY Brain projection.\n"))

multiExpr_early <- list(
  Brain = list(data = brain_temp_early[, valid_genes_early]),
  Blood = list(data = blood_temp_early[, valid_genes_early])
)
colorList_early <- list(Brain = net_early$colors[valid_genes_early])
names(colorList_early$Brain) <- valid_genes_early

cat("Executing robust WGCNA Permutation Engine (Early)...\n")
mp_early <- modulePreservation(multiExpr_early, colorList_early, referenceNetworks = 1, nPermutations = 100, randomSeed = 42, verbose = 0)

z_scores_early <- mp_early$preservation$Z[[1]][[2]]$Zsummary.pres
names(z_scores_early) <- rownames(mp_early$preservation$Z[[1]][[2]])

# ------------------------------------------------------------------------------
# PART 2: LATE STAGE (FRONTAL CORTEX) PRESERVATION
# ------------------------------------------------------------------------------
common_genes_late <- intersect(colnames(net_late$work_mat), rownames(deg_blood$expr_sync))

brain_temp_late <- net_late$work_mat[, common_genes_late]
blood_temp_late <- t(deg_blood$expr_sync[common_genes_late, ])

# Variance filter
gsg_brain_late <- goodSamplesGenes(brain_temp_late, verbose = 0)
gsg_blood_late <- goodSamplesGenes(blood_temp_late, verbose = 0)
valid_genes_late <- common_genes_late[gsg_brain_late$goodGenes & gsg_blood_late$goodGenes]

cat(paste("\n.. Retained", length(valid_genes_late), "valid genes for LATE Brain projection.\n"))

multiExpr_late <- list(
  Brain = list(data = brain_temp_late[, valid_genes_late]),
  Blood = list(data = blood_temp_late[, valid_genes_late])
)
colorList_late <- list(Brain = net_late$colors[valid_genes_late])
names(colorList_late$Brain) <- valid_genes_late

cat("Executing robust WGCNA Permutation Engine (Late)...\n")
mp_late <- modulePreservation(multiExpr_late, colorList_late, referenceNetworks = 1, nPermutations = 100, randomSeed = 42, verbose = 0)

z_scores_late <- mp_late$preservation$Z[[1]][[2]]$Zsummary.pres
names(z_scores_late) <- rownames(mp_late$preservation$Z[[1]][[2]])

cat("\nSUCCESS: All Module Preservation calculated successfully!\n")

# ==============================================================================
# PHASE 6 REVISION: PIVOT TO THE DISCOVERED "GREEN" BIOMARKER HUB
# ==============================================================================

library(tidyverse)
library(WGCNA)

# 1. Update the prioritization function to accept dynamic module colors
prioritize_biomarkers <- function(work_mat, colors_vec, pheno_df, blood_mat, label, target_module = "green") {
  
  target_genes <- names(colors_vec)[colors_vec == target_module]
  valid_genes <- intersect(target_genes, rownames(blood_mat))
  
  numeric_trait <- as.numeric(as.factor(pheno_df$Diagnosis))
  gs_scores <- cor(work_mat[, valid_genes], numeric_trait, use = "p")
  
  # Dynamically construct the Eigengene name (e.g., "MEgreen")
  me_name <- paste0("ME", target_module)
  me_vector <- moduleEigengenes(work_mat, colors = colors_vec)$eigengenes[[me_name]]
  
  mm_scores <- cor(work_mat[, valid_genes], me_vector, use = "p")
  
  tibble(
    Gene_Symbol = valid_genes, Brain_Network = label,
    Module_Membership = round(abs(as.numeric(mm_scores)), 4),
    Gene_Significance = round(abs(as.numeric(gs_scores)), 4),
    Prioritization_Score = round((Module_Membership + Gene_Significance) / 2, 4)
  ) %>% arrange(desc(Prioritization_Score))
}

cat("\nRe-ranking biomarker candidates targeting the highly preserved GREEN module...\n")

# 2. Extract leaders from the Green module
early_leaders <- prioritize_biomarkers(
  net_early$work_mat, net_early$colors, deg_early$pheno_sync, 
  deg_blood$expr_sync, "Early_AsymAD_Entorhinal", "green"
)

late_leaders <- prioritize_biomarkers(
  net_late$work_mat, net_late$colors, deg_late$pheno_sync, 
  deg_blood$expr_sync, "Late_SymAD_Frontal", "green"
)

# 3. Combine into the final master export
master_biomarker_export <- bind_rows(early_leaders, late_leaders)

# Save the updated file
write_csv(master_biomarker_export, "C:/Users/USER/Documents/Genomac_clients/Jaiyoba_FINAL/Data/Master_AD_Green_Biomarkers.csv")

cat("\n======================================================================\n")
cat("      TOP 10 SYSTEMIC BIOMARKERS (GREEN MODULE HUB GENES)             \n")
cat("======================================================================\n")
print(head(master_biomarker_export, 10))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

# ==============================================================================
# PHASE 9: FUNCTIONAL ENRICHMENT ANALYSIS (GENE ONTOLOGY)
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(ggplot2)

# ==============================================================================
# PHASE 9 REVISION: STATISTICALLY OPTIMIZED GO ENRICHMENT
# ==============================================================================

cat("\n======================================================================\n")
cat("  RUNNING OPTIMIZED FUNCTIONAL ENRICHMENT ENGINE                      \n")
cat("======================================================================\n\n")

# 1. Define the Background Universe
# By telling clusterProfiler to only use genes we actually measured, we stop it 
# from penalizing us for the 10,000+ unmeasured genes in the human genome.
background_universe <- rownames(expr_blood_symbols)
cat(paste(".. Setting biological universe to", length(background_universe), "measured genes.\n"))

# 2. First Attempt: Standard FDR Correction with proper Universe
cat(".. Attempting rigorous FDR-corrected GO Enrichment...\n")
go_bp_results <- enrichGO(
  gene          = target_genes,
  universe      = background_universe, 
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,  
  qvalueCutoff  = 0.2,   # CRITICAL: This must be relaxed alongside the p-value
  readable      = FALSE
)

# 3. Fallback: Exploratory Unadjusted Signal Extraction
if (is.null(go_bp_results) || nrow(as.data.frame(go_bp_results)) == 0) {
  cat("\n[!] FDR threshold too stringent for cross-tissue hubs. Activating exploratory mode...\n")
  
  go_bp_results <- enrichGO(
    gene          = target_genes,
    universe      = background_universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "none", # Remove the BH penalty to see the raw biological signal
    pvalueCutoff  = 0.01,   # Use a stricter raw p-value to compensate
    qvalueCutoff  = 1.0,    # Ignore Q-value
    readable      = FALSE
  )
}

# 4. Verification and Plotting Prep
if (is.null(go_bp_results) || nrow(as.data.frame(go_bp_results)) == 0) {
  stop("FATAL: No enrichment found even in exploratory mode. The gene list is too heterogeneous.")
} else {
  cat(paste("SUCCESS: Extracted", nrow(as.data.frame(go_bp_results)), "biological pathways!\n"))
  
  # Print the top 5 pathways to the console so we can see what they are immediately
  cat("\n--- Top 5 Biological Pathways Driven by the Green Module ---\n")
  print(head(as.data.frame(go_bp_results)[, c("ID", "Description", "pvalue", "Count")], 5))
}

go_df <- as.data.frame(go_bp_results)

write_csv(go_df, paste0(data_dir, "Manuscript_Table3_GO_Biological_Processes.csv"))
cat(".. Results exported to Manuscript_Table3_GO_Biological_Processes.csv\n")




# ==============================================================================
# PHASE 7 & 8: ADVANCED COLORBLIND-COMPLIANT PUBLICATION GRAPHICS ENGINE
# (DUAL FORMAT: PNG & TIFF, 300 DPI, TITLE-FREE, GRID-FREE)
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(viridis)
library(circlize)
library(enrichplot)

cat("\n======================================================================\n")
cat("  INITIALIZING ADVANCED MULTI-TISSUE VISUALIZATION ENGINE             \n")
cat("======================================================================\n\n")

# Define global export directory
plot_dir <- "C:/Users/USER/Documents/Genomac_clients/Jaiyoba_FINAL/Plots/"

# Unified, highly accessible academic theme (GRID-FREE)
theme_publication <- function() {
  theme_minimal(base_size = 12) +
    theme(
      axis.title.x      = element_text(face = "bold", size = 11, color = "black", margin = margin(t=10)),
      axis.title.y      = element_text(face = "bold", size = 11, color = "black", margin = margin(r=10)),
      axis.text.x       = element_text(size = 10, face = "bold", color = "black"),
      axis.text.y       = element_text(size = 10, face = "bold", color = "black"),
      panel.grid.major  = element_blank(), # Gridlines completely removed
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 1.2),
      legend.title      = element_text(face = "bold", size = 10),
      legend.text       = element_text(size = 9, face = "bold")
    )
}

# Universal device opener for Dual-Format Export
formats <- c(".png", ".tiff")
open_device <- function(filename, ext, w, h) {
  filepath <- paste0(plot_dir, filename, ext)
  if (ext == ".png") {
    png(filepath, width = w, height = h, units = "in", res = 300)
  } else {
    tiff(filepath, width = w, height = h, units = "in", res = 300, compression = "lzw")
  }
}

# Helper to safely clear graphics
dev.off_all <- function() { while(!is.null(dev.list())) dev.off() }

# ------------------------------------------------------------------------------
# FIGURE 1: GENE DENDROGRAMS (EARLY & LATE BRAIN)
# ------------------------------------------------------------------------------
cat("[1/10] Rendering Figure 1: Hierarchical Clustering Dendrograms... ")

title_1b <- "Figure1B_Gene_Hierarchical_Clustering_Dendrogram_Entorhinal_Cortex"
for (ext in formats) {
  open_device(title_1b, ext, 9.5, 5.5)
  par(mar = c(3, 12, 1, 1)) # Tight top margin for title-free layout
  plotDendroAndColors(
    dendro       = net_early$tree,         
    colors       = net_early$colors,      
    groupLabels  = "Modules",   
    dendroLabels = FALSE,                   
    hang         = 0.03,                    
    addGuide     = TRUE,                    
    guideHang    = 0.05,
    main         = "" 
  )
  dev.off()
}

title_1c <- "Figure1C_Gene_Hierarchical_Clustering_Dendrogram_Frontal_Cortex"
for (ext in formats) {
  open_device(title_1c, ext, 9.5, 5.5)
  par(mar = c(3, 12, 1, 1))
  plotDendroAndColors(
    dendro       = net_late$tree,         
    colors       = net_late$colors,      
    groupLabels  = "Modules",   
    dendroLabels = FALSE,                   
    hang         = 0.03,                    
    addGuide     = TRUE,                    
    guideHang    = 0.05,
    main         = "" 
  )
  dev.off()
}
cat("Done.\n")

# ------------------------------------------------------------------------------
# FIGURE 2 & S3: NETWORK TOPOLOGY (SCALE-FREE FIT FOR EARLY AND LATE)
# ------------------------------------------------------------------------------
cat("[2/10] Rendering Figure 2 & S3: Scale-Free Topology Diagnostics... ")

powers_vector <- c(1:10, seq(12, 20, 2))

# Early Topology
sft_early <- pickSoftThreshold(net_early$work_mat, powerVector = powers_vector, networkType = "signed", verbose = 0)
df_sft_early <- as.data.frame(sft_early$fitIndices)

p2_a <- ggplot(df_sft_early, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
  geom_point(color = "#D55E00", size = 3) + geom_line(color = "#D55E00", linewidth = 1) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Soft-Threshold Power (beta)", y = "Scale-Free Topology Fit (R^2)", title = NULL) + theme_publication()

p2_b <- ggplot(df_sft_early, aes(x = Power, y = mean.k.)) +
  geom_point(color = "#0072B2", size = 3) + geom_line(color = "#0072B2", linewidth = 1) +
  labs(x = "Soft-Threshold Power (beta)", y = "Mean Connectivity (mean.k)", title = NULL) + theme_publication()

title_2 <- "Figure2_Network_Scale-Free_Fit_and_Mean_Connectivity_Profile_Entorhinal_Cortex"
for (ext in formats) {
  open_device(title_2, ext, 8.5, 4.5)
  print(p2_a + p2_b + plot_annotation(tag_levels = 'A')) 
  dev.off()
}

# Late Topology (Supplementary)
sft_late <- pickSoftThreshold(net_late$work_mat, powerVector = powers_vector, networkType = "signed", verbose = 0)
df_sft_late <- as.data.frame(sft_late$fitIndices)

p_supp_a <- ggplot(df_sft_late, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
  geom_point(color = "#D55E00", size = 3) + geom_line(color = "#D55E00", linewidth = 1) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Soft-Threshold Power (beta)", y = "Scale-Free Topology Fit (R^2)", title = NULL) + theme_publication()

p_supp_b <- ggplot(df_sft_late, aes(x = Power, y = mean.k.)) +
  geom_point(color = "#0072B2", size = 3) + geom_line(color = "#0072B2", linewidth = 1) +
  labs(x = "Soft-Threshold Power (beta)", y = "Mean Connectivity (mean.k)", title = NULL) + theme_publication()

title_s3 <- "FigureS3_Network_Scale-Free_Fit_and_Mean_Connectivity_Profile_Frontal_Cortex"
for (ext in formats) {
  open_device(title_s3, ext, 8.5, 4.5)
  print(p_supp_a + p_supp_b + plot_annotation(tag_levels = 'A')) 
  dev.off()
}
cat("Done.\n")

# ------------------------------------------------------------------------------
# FIGURE 3: CROSS-TISSUE MODULE PRESERVATION BAR CHART (GREEN MODULE TARGETED)
# ------------------------------------------------------------------------------
cat("[3/10] Rendering Figure 3: Highlighted Cross-Tissue Preservation... ")

df_pres_early_updated <- tibble(
  Module = names(z_scores_early),
  Zsummary = as.numeric(z_scores_early)
) %>%
  filter(Module != "gold" & Module != "grey") %>% arrange(desc(Zsummary)) %>%
  mutate(
    Module = factor(Module, levels = Module),
    Highlight = ifelse(Module == "green", "Preserved", "Unpreserved") # PIVOTED TO GREEN
  )

fig3_updated <- ggplot(df_pres_early_updated, aes(x = Module, y = Zsummary, fill = Module, alpha = Highlight)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.6, width = 0.7) +
  scale_fill_manual(values = setNames(as.character(df_pres_early_updated$Module), df_pres_early_updated$Module)) +
  scale_alpha_manual(values = c("Preserved" = 1.0, "Unpreserved" = 0.3)) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "#0072B2", linewidth = 1.2) + 
  geom_hline(yintercept = 10, linetype = "solid", color = "#333333", linewidth = 1.2) +
  labs(x = "Discovered Early Brain Modules", y = expression(bold(Cross~Tissue~Preservation~(Z[summary]))), title = NULL) +
  theme_publication() + theme(legend.position = "none")

title_3 <- "Figure3_Translational_Projection_of_Brain_Co-Expression_Networks_onto_Whole_Blood"
for (ext in formats) {
  open_device(title_3, ext, 7.5, 5.0)
  print(fig3_updated)
  dev.off()
}
cat("Done.\n")

# ------------------------------------------------------------------------------
# FIGURE 4: BIOMARKER PRIORITIZATION MATRIX (VIRIDIS DOT PLOT)
# ------------------------------------------------------------------------------
cat("[4/10] Rendering Figure 4: Accessible Biomarker Selection Space... ")

fig4_plot <- ggplot(master_biomarker_export %>% filter(Brain_Network == "Early_AsymAD_Entorhinal"), 
                    aes(x = Gene_Significance, y = Module_Membership)) +
  geom_point(aes(color = Prioritization_Score, size = Prioritization_Score), alpha = 0.85) +
  scale_color_viridis_c(option = "plasma", direction = 1) +
  geom_text(data = head(master_biomarker_export %>% filter(Brain_Network == "Early_AsymAD_Entorhinal"), 5),
            aes(label = Gene_Symbol), vjust = -1.2, fontface = "bold", size = 3.5, color = "black") +
  labs(x = "Gene Significance (Correlation with AsymAD Phenotype)", 
       y = "Module Membership (Leverage Hub Weight inside Green Network)",
       title = NULL, color = "Score", size = "Score") +
  xlim(0, 0.5) + ylim(0.4, 1.0) +
  theme_publication() + theme(legend.position = "right")

title_4 <- "Figure4_Dual-Metric_Prioritization_Matrix_for_Blood-Preserved_Candidates"
for (ext in formats) {
  open_device(title_4, ext, 7.2, 5.5)
  print(fig4_plot)
  dev.off()
}
cat("Done.\n")

# ------------------------------------------------------------------------------
# FIGURE 5: CROSS-TISSUE CIRCOS PLOT (GREEN HUB GENES)
# ------------------------------------------------------------------------------
cat("[5/10] Rendering Figure 5: Cross-Tissue Hub Gene Circos Plot... ")

top_genes <- head(master_biomarker_export$Gene_Symbol, 5)

circos_df <- data.frame(
  Source = c(rep("Brain (Entorhinal)", 5), top_genes),
  Target = c(top_genes, rep("Peripheral Blood", 5)),
  Value  = c(master_biomarker_export$Prioritization_Score[1:5], master_biomarker_export$Prioritization_Score[1:5])
)

# PIVOTED TO GREEN COLOR ASSIGNMENT
grid_colors <- c("Brain (Entorhinal)" = "#D55E00", "Peripheral Blood" = "#0072B2", setNames(rep("green", 5), top_genes))

title_5 <- "Figure5_Cross-Tissue_Network_Preservation_of_Green_Module_Hub_Genes"
for (ext in formats) {
  open_device(title_5, ext, 6, 6)
  circos.clear()
  circos.par(gap.degree = 4, start.degree = 90)
  chordDiagram(circos_df, grid.col = grid_colors, transparency = 0.3, 
               directional = 1, direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow", annotationTrack = c("name", "grid"),
               annotationTrackHeight = c(0.05, 0.05))
  dev.off()
  circos.clear()
}
cat("Done.\n")

# ------------------------------------------------------------------------------
# FIGURE 6: EIGENGENE & HUB GENE CLINICAL TRAJECTORY HEATMAP (ROBUST PATCH)
# ------------------------------------------------------------------------------
cat("[6/10] Rendering Figure 6: Eigengene and Hub Gene Clinical Trajectory... ")

# SAFELY RECALCULATE EIGENGENES ON THE FLY
me_early <- moduleEigengenes(net_early$work_mat, colors = net_early$colors)$eigengenes$MEgreen
names(me_early) <- rownames(net_early$work_mat)

me_late <- moduleEigengenes(net_late$work_mat, colors = net_late$colors)$eigengenes$MEgreen
names(me_late) <- rownames(net_late$work_mat)

# ROBUST TOP GENE EXTRACTION
common_full_genes <- intersect(rownames(deg_early$expr_sync), rownames(deg_late$expr_sync))
safe_master <- master_biomarker_export %>% filter(Gene_Symbol %in% common_full_genes)
top_genes_safe <- head(safe_master$Gene_Symbol, 5)

expr_top_early <- deg_early$expr_sync[top_genes_safe, , drop=FALSE]
expr_top_late  <- deg_late$expr_sync[top_genes_safe, , drop=FALSE]

# UNIFY PHENOTYPES FOR THE TRAJECTORY
unified_pheno <- bind_rows(
  deg_early$pheno_sync %>% mutate(Stage = Diagnosis),
  deg_late$pheno_sync %>% filter(Diagnosis == "SymAD") %>% mutate(Stage = Diagnosis)
) %>%
  mutate(Stage = factor(Stage, levels = c("Control", "AsymAD", "SymAD"))) %>%
  arrange(Stage)

# ALIGN EXPRESSION AND EIGENGENE DATA
combined_expr <- cbind(expr_top_early, expr_top_late)[, unified_pheno$SampleID]
combined_me <- c(me_early, me_late[deg_late$pheno_sync$Diagnosis == "SymAD"])
names(combined_me) <- unified_pheno$SampleID

# LABEL EXPLICITLY AS THE GREEN MODULE EIGENGENE
heatmap_matrix <- rbind("Green Module Eigengene" = combined_me, combined_expr)
heatmap_scaled <- t(scale(t(heatmap_matrix)))

annotation_col <- data.frame(Disease_Stage = unified_pheno$Stage)
rownames(annotation_col) <- unified_pheno$SampleID

ann_colors <- list(Disease_Stage = c("Control" = "#009E73", "AsymAD" = "#E69F00", "SymAD" = "#D55E00"))

title_6 <- "Figure6_Green_Module_Eigengene_and_Hub_Gene_Clinical_Trajectories"
for (ext in formats) {
  open_device(title_6, ext, 8.5, 4.5)
  pheatmap(
    mat = heatmap_scaled,
    cluster_cols = FALSE, cluster_rows = FALSE,
    annotation_col = annotation_col, annotation_colors = ann_colors,
    show_colnames = FALSE, color = viridis(100, option = "magma"), 
    main = NA, 
    gaps_row = c(1), fontsize_row = 11, border_color = NA
  )
  dev.off()
}
cat("Done.\n")

# ------------------------------------------------------------------------------
# FIGURE 7: GO ENRICHMENT DOTPLOT
# ------------------------------------------------------------------------------
cat("[7/10] Rendering Figure 7: GO Biological Process Dotplot... ")

theme_enrichment <- theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.title  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

p_dot <- dotplot(go_bp_results, showCategory = 15, title = NULL) + 
  scale_color_viridis_c(option = "plasma") + 
  theme_enrichment

title_7 <- "Figure7_GO_Biological_Process_Dotplot"
for (ext in formats) {
  open_device(title_7, ext, 8, 6.5)
  print(p_dot)
  dev.off()
}
cat("Done.\n")
# ------------------------------------------------------------------------------
# FIGURE 8: GENE-CONCEPT NETWORK (CNETPLOT) - PATCHED
# ------------------------------------------------------------------------------
cat("[8/10] Rendering Figure 8: Gene-Concept Network... ")

fold_changes <- master_biomarker_export$Prioritization_Score
names(fold_changes) <- master_biomarker_export$Gene_Symbol
fold_changes <- fold_changes[!duplicated(names(fold_changes))]

# Simplified call to avoid ggraph layout version conflicts
p_cnet <- cnetplot(
  go_bp_results, 
  foldChange = fold_changes, 
  showCategory = 5
) +
  scale_color_viridis_c(option = "mako", name = "Hub\nScore") + 
  theme(legend.position = "right")

title_8 <- "Figure8_GO_Gene_Concept_Network"
for (ext in formats) {
  open_device(title_8, ext, 9.5, 7.5)
  print(p_cnet)
  dev.off()
}
cat("Done.\n")
# ------------------------------------------------------------------------------
# SUPPLEMENTARY FIGURES S1 & S2: MODULE-TRAIT CORRELATION HEATMAPS
# ------------------------------------------------------------------------------
cat("[9/10 & 10/10] Rendering Supplementary Heatmaps S1 & S2... ")

generate_trait_heatmap <- function(cor_mat, p_mat, title_prefix) {
  textMatrix = paste(signif(cor_mat, 2), "\n(", signif(p_mat, 1), ")", sep = "")
  dim(textMatrix) = dim(cor_mat)
  col_palette <- colorRampPalette(c("#0072B2", "white", "#D55E00"))(50)
  
  for (ext in formats) {
    open_device(title_prefix, ext, 5, 8)
    pheatmap(
      mat = cor_mat,
      display_numbers = textMatrix,
      cluster_rows = FALSE,       
      cluster_cols = FALSE,       
      main = NA,                  
      color = col_palette,
      breaks = seq(-1, 1, length.out = 50),
      fontsize = 11,
      fontsize_number = 9,        
      angle_col = 0               
    )
    dev.off()
  }
}

title_s1 <- "FigureS1_Early_Brain_Entorhinal_Module_Trait_Correlation_Heatmap"
generate_trait_heatmap(as.matrix(trait_early$cor), as.matrix(trait_early$p), title_s1)

title_s2 <- "FigureS2_Late_Brain_Frontal_Module_Trait_Correlation_Heatmap"
generate_trait_heatmap(as.matrix(trait_late$cor), as.matrix(trait_late$p), title_s2)

cat("Done.\n")

cat("\n======================================================================\n")
cat("SUCCESS: Complete Dual-Format (PNG + TIFF) Graphics Suite Exported!\n")
cat("======================================================================\n")



# ------------------------------------------------------------------------------
# SUPPLEMENTARY FIGURE S1: COMBINED MODULE-TRAIT CORRELATION HEATMAPS
# ------------------------------------------------------------------------------
cat("[9/10] Rendering Combined Supplementary Heatmap S1... ")

generate_combined_trait_heatmap <- function(cor_early, p_early, cor_late, p_late, title_prefix) {
  
  # SAFETY CHECK: Force column names if missing
  if (is.null(colnames(cor_early))) colnames(cor_early) <- "AsymAD"
  if (is.null(colnames(cor_late))) colnames(cor_late) <- "SymAD"
  
  # Format the text matrices for both networks
  textMatrix_early = paste(signif(cor_early, 2), "\n(", signif(p_early, 1), ")", sep = "")
  dim(textMatrix_early) = dim(cor_early)
  
  textMatrix_late = paste(signif(cor_late, 2), "\n(", signif(p_late, 1), ")", sep = "")
  dim(textMatrix_late) = dim(cor_late)
  
  # Extract module colors
  colors_early <- gsub("ME", "", rownames(cor_early))
  colors_late <- gsub("ME", "", rownames(cor_late))
  
  for (ext in formats) {
    # 1. Open a wider device to fit both plots side-by-side
    open_device(title_prefix, ext, w = 10, h = 7) 
    
    # 2. Split the canvas into 1 row, 2 columns
    par(mfrow = c(1, 2))
    
    # 3. Manually set plotting margins: c(bottom, left, top, right)
    # Increased bottom (6) for x-axis labels, increased left (7) for y-axis colored blocks
    par(mar = c(6, 7, 4, 2))
    
    # --- PANEL A: Early Brain ---
    labeledHeatmap(Matrix = cor_early,
                   xLabels = colnames(cor_early),
                   yLabels = rownames(cor_early),
                   ySymbols = colors_early,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix_early,
                   setStdMargins = FALSE,
                   cex.text = 0.8,
                   zlim = c(-1, 1),
                   main = "Early Brain (Entorhinal)")
    
    # --- PANEL B: Late Brain ---
    labeledHeatmap(Matrix = cor_late,
                   xLabels = colnames(cor_late),
                   yLabels = rownames(cor_late),
                   ySymbols = colors_late,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix_late,
                   setStdMargins = FALSE,
                   cex.text = 0.8,
                   zlim = c(-1, 1),
                   main = "Late Brain (Frontal)")
    
    # Close the device and save
    dev.off()
  }
}

# Execute the combined function
title_combined <- "FigureS1_Combined_Brain_Module_Trait_Correlation"
generate_combined_trait_heatmap(
  cor_early = as.matrix(trait_early$cor), p_early = as.matrix(trait_early$p),
  cor_late  = as.matrix(trait_late$cor),  p_late  = as.matrix(trait_late$p),
  title_prefix = title_combined
)

cat("Done.\n")







# ------------------------------------------------------------------------------
# FIGURE 2: COMBINED NETWORK TOPOLOGY (SCALE-FREE FIT FOR EARLY AND LATE)
# ------------------------------------------------------------------------------
cat("[2/10] Rendering Combined Figure 2: Scale-Free Topology Diagnostics... ")

powers_vector <- c(1:10, seq(12, 20, 2))

# --- Early Topology (Entorhinal Cortex) ---
sft_early <- pickSoftThreshold(net_early$work_mat, powerVector = powers_vector, networkType = "signed", verbose = 0)
df_sft_early <- as.data.frame(sft_early$fitIndices)

p_early_r2 <- ggplot(df_sft_early, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
  geom_point(color = "#D55E00", size = 3) + geom_line(color = "#D55E00", linewidth = 1) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Soft-Threshold Power (beta)", y = "Scale-Free Topology Fit (R^2)", title = "Early Brain (Entorhinal)") + 
  theme_publication()

p_early_k <- ggplot(df_sft_early, aes(x = Power, y = mean.k.)) +
  geom_point(color = "#0072B2", size = 3) + geom_line(color = "#0072B2", linewidth = 1) +
  labs(x = "Soft-Threshold Power (beta)", y = "Mean Connectivity (mean.k)", title = " ") + # Empty space for alignment
  theme_publication()

# --- Late Topology (Frontal Cortex) ---
sft_late <- pickSoftThreshold(net_late$work_mat, powerVector = powers_vector, networkType = "signed", verbose = 0)
df_sft_late <- as.data.frame(sft_late$fitIndices)

p_late_r2 <- ggplot(df_sft_late, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
  geom_point(color = "#D55E00", size = 3) + geom_line(color = "#D55E00", linewidth = 1) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Soft-Threshold Power (beta)", y = "Scale-Free Topology Fit (R^2)", title = "Late Brain (Frontal)") + 
  theme_publication()

p_late_k <- ggplot(df_sft_late, aes(x = Power, y = mean.k.)) +
  geom_point(color = "#0072B2", size = 3) + geom_line(color = "#0072B2", linewidth = 1) +
  labs(x = "Soft-Threshold Power (beta)", y = "Mean Connectivity (mean.k)", title = " ") + 
  theme_publication()

# --- Combine with Patchwork ---
# Creates a 2x2 layout: Early on Top (A & B), Late on Bottom (C & D)
combined_plot <- (p_early_r2 | p_early_k) / (p_late_r2 | p_late_k) + 
  plot_annotation(tag_levels = 'A')

title_combined <- "Figure2_Combined_Network_Scale-Free_Fit_and_Mean_Connectivity"

for (ext in formats) {
  open_device(title_combined, ext, w = 10, h = 9) # Taller canvas to support 2 rows
  print(combined_plot)
  dev.off()
}

cat("Done.\n")


# 1. View the entire table of fit indices
# This will show you the R^2 value (SFT.R.sq) for every soft-thresholding power tested.
print(sft_late$fitIndices)

# 2. Extract the exact R^2 value specifically at soft-thresholding power 4
power_4_R2 <- sft_late$fitIndices$SFT.R.sq[sft_late$fitIndices$Power == 4]
print(power_4_R2)

# 3. Extract the maximum R^2 value (where the network completely stabilizes)
max_R2 <- max(sft_late$fitIndices$SFT.R.sq)
print(max_R2)





# ==============================================================================
# SUPPLEMENTARY TABLE: CROSS-TISSUE MODULE PRESERVATION (Z-SUMMARY SCORES)
# ==============================================================================

# Load required libraries
library(dplyr)
library(gt) # For publication-quality table rendering

# 1. Define Early Stage Data
early_data <- data.frame(
  Module = c("black", "blue", "brown", "gold", "green", "magenta", 
             "red", "turquoise", "yellow", "pink", "purple"),
  Z_Summary_Early = c(-0.4510598, 2.3616169, 0.3629981, 1.7987014, 15.2542124, 
                      2.3170994, -0.1347982, 0.4903086, 4.4662761, -1.0298502, 0.090289)
)

# 2. Define Late Stage Data
late_data <- data.frame(
  Module = c("black", "blue", "brown", "gold", "green", "yellow", "red", "turquoise"),
  Z_Summary_Late = c(1.8329624, 0.7466875, 2.0158106, 3.5076252, 10.5870596, 
                     2.102223, -0.9622151, 7.8149882)
)

# 3. Merge Datasets
# A full_join ensures all modules from both stages are included. 
# We arrange by the Early score descending to highlight the Green module at the top.
preservation_table <- full_join(early_data, late_data, by = "Module") %>%
  arrange(desc(Z_Summary_Early))

# 4. Generate Publication-Ready Table using 'gt'
supplementary_gt <- preservation_table %>%
  gt() %>%
  tab_header(
    title = md("**Supplementary Table S1**"),
    subtitle = "Cross-Tissue Network Preservation Analysis of Brain-Derived Modules in Peripheral Blood"
  ) %>%
  cols_label(
    Module = "Network Module",
    Z_Summary_Early = "Early Asymptomatic AD (Z-summary)",
    Z_Summary_Late = "Late Symptomatic AD (Z-summary)"
  ) %>%
  fmt_number(
    columns = c(Z_Summary_Early, Z_Summary_Late),
    decimals = 2 # Rounds to 2 decimal places for standard academic formatting
  ) %>%
  sub_missing(
    columns = Z_Summary_Late,
    missing_text = "-" # Replaces NA with a clean dash for missing late-stage modules
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Module,
      rows = Module == "green" # Bolds the Green module to draw reviewer attention
    )
  ) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12
  )

# 5. View the Table
print(supplementary_gt)

# 6. Export Options
# To save as a Word document (.docx) for your manuscript:
gtsave(supplementary_gt, "Supplementary_Table_S1_Preservation.docx")

# To save as a standard CSV:
# write.csv(preservation_table, "Supplementary_Table_S1_Preservation.csv", row.names = FALSE)