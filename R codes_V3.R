# ==============================================================================
# DESeq2 Analysis: Control A vs ABC, A vs AB, and A vs AC
# ==============================================================================

# 1. Install & Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("apeglm", quietly = TRUE)) BiocManager::install("apeglm")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

library(DESeq2)
library(ggplot2)
library(pheatmap)

# 2. Load Raw Counts File
counts_raw <- read.csv("raw_counts.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Set feature_id as row names and extract count matrix
rownames(counts_raw) <- counts_raw$feature_id
count_matrix <- as.matrix(counts_raw[, -1])

# Round fractional counts to integers (required for DESeq2's negative binomial model)
count_matrix_int <- round(count_matrix)

# 3. Build Sample Metadata (colData)
sample_names <- colnames(count_matrix_int)

# Extract group prefix from sample names (e.g., "A-1" -> "A", "ABC-2" -> "ABC")
group_labels <- gsub("-[0-9]+$", "", sample_names)

colData <- data.frame(
  sample_id = sample_names,
  group = factor(group_labels),
  row.names = sample_names
)

# Set 'A' as the control / reference group baseline
colData$group <- relevel(colData$group, ref = "A")

# 4. Construct DESeqDataSet and Run Pipeline
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_int,
  colData = colData,
  design = ~ group
)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# 5. Extract Differential Expression Results with LFC Shrinkage

# --- Contrast 1: ABC vs A (A as reference control) ---
res_ABC_vs_A <- results(dds, contrast = c("group", "ABC", "A"), alpha = 0.05)
res_ABC_vs_A_shunk <- lfcShrink(dds, contrast = c("group", "ABC", "A"), type = "ashr")

# --- Contrast 2: AB vs A (A as reference control) ---
res_AB_vs_A <- results(dds, contrast = c("group", "AB", "A"), alpha = 0.05)
res_AB_vs_A_shunk <- lfcShrink(dds, contrast = c("group", "AB", "A"), type = "ashr")

# --- Contrast 3: AC vs A (A as reference control) ---
res_AC_vs_A <- results(dds, contrast = c("group", "AC", "A"), alpha = 0.05)
res_AC_vs_A_shunk <- lfcShrink(dds, contrast = c("group", "AC", "A"), type = "ashr")

# 6. Save Results to CSV Files
write.csv(as.data.frame(res_ABC_vs_A_shunk), "DE_results_ABC_vs_A.csv")
write.csv(as.data.frame(res_AB_vs_A_shunk),  "DE_results_AB_vs_A.csv")
write.csv(as.data.frame(res_AC_vs_A_shunk),  "DE_results_AC_vs_A.csv")

cat("Differential expression results successfully exported!\n")

# ==============================================================================
# Visualizations
# ==============================================================================

# A. Principal Component Analysis (PCA)
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3) +
  labs(title = "PCA Plot of Samples",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

ggsave("PCA_plot.png", width = 7, height = 5)

# B. Sample Distance Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample Distance Clustering",
         filename = "Sample_Distance_Heatmap.png")

# C. MA Plots for all Comparisons
png("MA_Plots.png", width = 1200, height = 400)
par(mfrow = c(1, 3))
plotMA(res_ABC_vs_A_shunk, main = "ABC vs Control A", ylim = c(-5, 5))
plotMA(res_AB_vs_A_shunk,  main = "AB vs Control A",  ylim = c(-5, 5))
plotMA(res_AC_vs_A_shunk,  main = "AC vs Control A",  ylim = c(-5, 5))
dev.off()



















### with control gruop N

# ==============================================================================
# DESeq2 Analysis with Normalized Counts Export (Control = N)
# ==============================================================================

# 1. Load Required Packages
library(DESeq2)
library(ggplot2)
library(pheatmap)

# 2. Load and Prepare Raw Counts
counts_raw <- read.csv("raw_counts.csv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts_raw) <- counts_raw$feature_id
count_matrix <- as.matrix(counts_raw[, -1])

# Round fractional counts to integers for DESeq2 negative binomial modeling
count_matrix_int <- round(count_matrix)

# 3. Construct Metadata (colData)
sample_names <- colnames(count_matrix_int)
group_labels <- gsub("-[0-9]+$", "", sample_names)

colData <- data.frame(
  sample_id = sample_names,
  group = factor(group_labels),
  row.names = sample_names
)

# Set group 'N' as the reference control baseline
colData$group <- relevel(colData$group, ref = "N")

# 4. Construct DESeqDataSet and Run Pipeline
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_int,
  colData = colData,
  design = ~ group
)

# Run main DESeq2 workflow (calculates size factors, estimates dispersion, fits GLM)
dds <- DESeq(dds)

# ==============================================================================
# 5. Extract and Download Normalized Counts
# ==============================================================================

# A. Median-of-Ratios Normalized Counts (Primary Normalized Table)
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts), "DESeq2_normalized_counts.csv")

# B. Variance Stabilized Transformed (VST) Counts (Ideal for Heatmaps/PCA)
# ==============================================================================
# Fixed VST / Normalized Matrix Extraction for Small Datasets (<1,000 features)
# ==============================================================================

# Option 1: Use varianceStabilizingTransformation with fitType = "mean"
vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "mean")
vst_matrix <- assay(vsd)

# Export VST Normalized Counts to CSV
write.csv(as.data.frame(vst_matrix), "DESeq2_vst_normalized_counts.csv")

# ==============================================================================
# 6. Extract Differential Expression Results (LFC Shrinkage)
# ==============================================================================

res_ABC_vs_N <- lfcShrink(dds, contrast = c("group", "ABC", "N"), type = "ashr")
res_AB_vs_N  <- lfcShrink(dds, contrast = c("group", "AB",  "N"), type = "ashr")
res_AC_vs_N  <- lfcShrink(dds, contrast = c("group", "AC",  "N"), type = "ashr")
res_A_vs_N   <- lfcShrink(dds, contrast = c("group", "A",   "N"), type = "ashr")

# Export Differential Expression CSV Tables
write.csv(as.data.frame(res_ABC_vs_N), "DE_results_ABC_vs_N.csv")
write.csv(as.data.frame(res_AB_vs_N),  "DE_results_AB_vs_N.csv")
write.csv(as.data.frame(res_AC_vs_N),  "DE_results_AC_vs_N.csv")
write.csv(as.data.frame(res_A_vs_N),   "DE_results_A_vs_N.csv")

# ==============================================================================
# 7. Visualizations
# ==============================================================================

# PCA Plot
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3) +
  labs(title = "PCA Plot of Samples",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

ggsave("PCA_plot.png", width = 7, height = 5)

# Sample Distance Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample Distance Clustering",
         filename = "Sample_Distance_Heatmap.png")
















############## Correlation analysis #################


# ==============================================================================
# Script: Spearman Correlation between key lncRNAs and key target genes
# ==============================================================================

# 1. Load required libraries
required_packages <- c("tidyverse", "Hmisc", "pheatmap", "ggcorrplot", "reshape2", "svglite")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(Hmisc)
library(pheatmap)
library(ggcorrplot)
library(reshape2)
library(svglite)

# 2. Load Both Datasets
lnc_raw  <- read.csv("raw_counts_key_lncRNA.csv", check.names = FALSE)
gene_raw <- read.csv("raw_count_key_genes.csv", check.names = FALSE)

# 3. Format & Transpose (Samples as Rows, Features as Columns)
lnc_mat <- lnc_raw %>%
  column_to_rownames(var = "feature_id") %>%
  t() %>%
  as.data.frame()

gene_mat <- gene_raw %>%
  column_to_rownames(var = "gene_name") %>%
  t() %>%
  as.data.frame()

# Align sample order
common_samples <- intersect(rownames(lnc_mat), rownames(gene_mat))
lnc_mat  <- lnc_mat[common_samples, ]
gene_mat <- gene_mat[common_samples, ]

# Filter out zero-variance features
lnc_filtered  <- lnc_mat[, apply(lnc_mat, 2, var, na.rm = TRUE) > 0]
gene_filtered <- gene_mat[, apply(gene_mat, 2, var, na.rm = TRUE) > 0]

# 4. Compute Spearman Correlation (lncRNAs vs Key Genes)
lnc_names  <- colnames(lnc_filtered)
gene_names <- colnames(gene_filtered)

r_matrix <- matrix(NA, nrow = length(lnc_names), ncol = length(gene_names),
                   dimnames = list(lnc_names, gene_names))
p_matrix <- matrix(NA, nrow = length(lnc_names), ncol = length(gene_names),
                   dimnames = list(lnc_names, gene_names))

for (l in lnc_names) {
  for (g in gene_names) {
    test <- cor.test(lnc_filtered[[l]], gene_filtered[[g]], method = "spearman", exact = FALSE)
    r_matrix[l, g] <- test$estimate
    p_matrix[l, g] <- test$p.value
  }
}

# 5. Save Output CSV Files
write.csv(r_matrix, "lncRNA_vs_genes_spearman_r_values.csv")
write.csv(p_matrix, "lncRNA_vs_genes_spearman_p_values.csv")

# Create a clean pairwise summary table
melted_r <- melt(r_matrix, varnames = c("lncRNA", "Target_Gene"), value.name = "Spearman_r")
melted_p <- melt(p_matrix, varnames = c("lncRNA", "Target_Gene"), value.name = "p_value")

summary_table <- melted_r %>%
  mutate(
    p_value = melted_p$p_value,
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "NS"
    )
  ) %>%
  arrange(p_value)

write.csv(summary_table, "lncRNA_vs_genes_pairwise_summary.csv", row.names = FALSE)
cat("CSV files exported successfully!\n")

# 6. Plot 1: Cross-Correlation Heatmap (lncRNAs x Genes) with pheatmap
svg("lncRNA_vs_genes_heatmap.svg", width = 8, height = 6)
pheatmap(
  r_matrix,
  display_numbers = matrix(
    ifelse(p_matrix < 0.001, "***",
           ifelse(p_matrix < 0.01, "**",
                  ifelse(p_matrix < 0.05, "*", ""))),
    nrow = nrow(p_matrix)
  ),
  fontsize_number = 14,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  main = "lncRNA vs Target Gene Spearman Correlation (* p<0.05, ** p<0.01, *** p<0.001)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45
)
dev.off()

# 7. Plot 2: Scatter Plot with Linear Trendline and Correlation Stats (Example)
# Select top significant pair for demonstration
top_pair <- summary_table[1, ]

scatter_df <- data.frame(
  lncRNA_expr = lnc_filtered[[top_pair$lncRNA]],
  Gene_expr  = gene_filtered[[top_pair$Target_Gene]]
)

p_scatter <- ggplot(scatter_df, aes(x = lncRNA_expr, y = Gene_expr)) +
  geom_point(size = 3.5, shape = 21, fill = NA, color = "black", stroke = 1.2) +
  geom_smooth(method = "lm", color = "#B2182B", fill = "#B2182B", alpha = 0.15, linewidth = 1.1) +
  labs(
    title = paste("Correlation:", top_pair$lncRNA, "vs", top_pair$Target_Gene),
    subtitle = paste0("Spearman r = ", round(top_pair$Spearman_r, 3), 
                      ", p-value = ", format.pval(top_pair$p_value, digits = 3)),
    x = paste(top_pair$lncRNA, "Raw Count"),
    y = paste(top_pair$Target_Gene, "Raw Count")
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2)
  )

print(p_scatter)

ggsave(
  filename = "top_lncRNA_gene_scatterplot.svg", 
  plot = p_scatter, 
  device = svglite::svglite, 
  width = 6, 
  height = 5
)



######### For high quality Heatmap from the Spearman correlation ##########


# ==============================================================================
# Fixed Script: Heatmap Export without 0-byte issues
# ==============================================================================

library(tidyverse)
library(reshape2)
library(svglite)

# 1. Read files
lnc_raw  <- read.csv("raw_counts_key_lncRNA.csv", check.names = FALSE)
gene_raw <- read.csv("raw_count_key_genes.csv", check.names = FALSE)

# 2. Format & Transpose
lnc_mat <- lnc_raw %>% column_to_rownames(var = "feature_id") %>% t() %>% as.data.frame()
gene_mat <- gene_raw %>% column_to_rownames(var = "gene_name") %>% t() %>% as.data.frame()

# Align sample order
common_samples <- intersect(rownames(lnc_mat), rownames(gene_mat))
lnc_mat  <- lnc_mat[common_samples, ]
gene_mat <- gene_mat[common_samples, ]

# Remove zero variance features
lnc_filtered  <- lnc_mat[, apply(lnc_mat, 2, var, na.rm = TRUE) > 0]
gene_filtered <- gene_mat[, apply(gene_mat, 2, var, na.rm = TRUE) > 0]

# 3. Calculate Spearman Correlations & P-values
lnc_names  <- colnames(lnc_filtered)
gene_names <- colnames(gene_filtered)

r_matrix <- matrix(NA, nrow = length(lnc_names), ncol = length(gene_names),
                   dimnames = list(lnc_names, gene_names))
p_matrix <- matrix(NA, nrow = length(lnc_names), ncol = length(gene_names),
                   dimnames = list(lnc_names, gene_names))

for (l in lnc_names) {
  for (g in gene_names) {
    test <- cor.test(lnc_filtered[[l]], gene_filtered[[g]], method = "spearman", exact = FALSE)
    r_matrix[l, g] <- test$estimate
    p_matrix[l, g] <- test$p.value
  }
}

# 4. Save CSV files
write.csv(r_matrix, "lncRNA_vs_genes_spearman_r_values.csv")
write.csv(p_matrix, "lncRNA_vs_genes_spearman_p_values.csv")

# 5. Build ggplot2 Heatmap
df_r <- melt(r_matrix, varnames = c("lncRNA", "Gene"), value.name = "r_value")
df_p <- melt(p_matrix, varnames = c("lncRNA", "Gene"), value.name = "p_value")
df_plot <- cbind(df_r, p_value = df_p$p_value)

df_plot$stars <- ifelse(df_plot$p_value < 0.001, "***",
                        ifelse(df_plot$p_value < 0.01, "**",
                               ifelse(df_plot$p_value < 0.05, "*", "")))

p_heatmap <- ggplot(df_plot, aes(x = Gene, y = lncRNA, fill = r_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(r_value, 2), "\n", stars)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limit = c(-1, 1), name = "Spearman r") +
  labs(title = "Spearman Correlation: lncRNA vs Target Genes", x = "Target Gene", y = "lncRNA ID") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid = element_blank()
  )

# 6. Safe Save via ggsave (prevents zero-byte file creation)
ggsave(
  filename = "lncRNA_vs_genes_heatmap.svg",
  plot = p_heatmap,
  device = svglite::svglite,
  width = 8,
  height = 6,
  dpi = 300
)

# Optional PNG save for quick preview
ggsave(
  filename = "lncRNA_vs_genes_heatmap.png",
  plot = p_heatmap,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Successfully saved non-zero SVG and PNG heatmaps!\n")














#######################################################################################
New modified analysis for bar plot and spearman correlation analysis with CPM values
#######################################################################################



#group_colors <- c("ABC" = "#8E63CE", "AB"  = "#3C78D8", "AC"  = "#2AA198", "A"   = "#E69138", "N"   = "#7F7F7F")
#group_colors <- c("ABC" = "#C5A7EA", "AB"  = "#83C5F7", "AC"  = "#8ED9BA", "A"   = "#FBC078", "N"   = "#BFC3C7")
group_colors <- c("ABC" = "#C5A7EA", "AB"  = "#319EFC", "AC"  = "#4CBD8F", "A"   = "#FD9D30", "N"   = "#BFC3C7")


# ==============================================================================
# COMPLETE ANALYSIS
# Target ncRNAs and key target genes using CPM expression
#
# INPUT FILES:
#   1. target_ncRNAs_CPM.csv
#   2. Key_genes_CPM.csv
#
# ANALYSES:
#   A. Spearman correlation: ncRNAs vs key target genes
#   B. BH-FDR correction for multiple correlations
#   C. Correlation heatmap
#   D. Scatterplots for significant ncRNA-gene pairs
#   E. Bar/point plots of ncRNA expression among groups
#   F. ABC vs A statistical comparison
#
# IMPORTANT:
#   - Correlation is calculated using CPM values.
#   - log2(CPM + 1) is used for visualization.
#   - AB and AC have n = 1 and are NOT statistically tested.
# ==============================================================================











# ==============================================================================
# EXPRESSION PLOTS FOR CANDIDATE ncRNAs AND KEY TARGET GENES
# Groups: ABC, AB, AC, A, N
# Statistics: ABC vs A, Wilcoxon rank-sum test + BH correction
# ==============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(svglite)

# ==============================================================================
# 1. READ DATA
# ==============================================================================

lnc_raw <- read.csv("target_ncRNAs_CPM.csv", check.names = FALSE, stringsAsFactors = FALSE)
gene_raw <- read.csv("Key_genes_CPM.csv", check.names = FALSE, stringsAsFactors = FALSE)

# ==============================================================================
# 2. GROUPS AND COLORS
# ==============================================================================

plot_groups <- c("ABC", "AB", "AC", "A", "N")

group_colors <- c("ABC" = "#C5A7EA", "AB" = "#319EFC", "AC" = "#4CBD8F", "A" = "#FD9D30", "N" = "#BFC3C7")


# ==============================================================================
# PART A: CANDIDATE ncRNAs
# ==============================================================================

# ==============================================================================
# 3. CONVERT ncRNA DATA TO LONG FORMAT
# ==============================================================================

lnc_long <- lnc_raw %>%
  pivot_longer(cols = -feature_id, names_to = "Sample", values_to = "CPM") %>%
  mutate(Group = sub("-[0-9]+$", "", Sample), log2CPM = log2(CPM + 1))

lnc_plot_data <- lnc_long %>%
  filter(Group %in% plot_groups) %>%
  mutate(Group = factor(Group, levels = plot_groups))

cat("\nSample numbers used for ncRNA plots:\n")
print(lnc_plot_data %>% distinct(Sample, Group) %>% count(Group))

write.csv(lnc_plot_data, "ncRNA_sample_level_log2CPM_all_groups.csv", row.names = FALSE)


# ==============================================================================
# 4. ncRNA GROUP MEAN ± SEM
# ==============================================================================

lnc_summary <- lnc_plot_data %>%
  group_by(feature_id, Group) %>%
  summarise(n = n(),
            Mean = mean(log2CPM, na.rm = TRUE),
            SD = ifelse(n() > 1, sd(log2CPM, na.rm = TRUE), NA_real_),
            SEM = ifelse(n() > 1, SD / sqrt(n()), NA_real_),
            .groups = "drop")

write.csv(lnc_summary, "ncRNA_group_mean_SEM_all_groups.csv", row.names = FALSE)


# ==============================================================================
# 5. ncRNA STATISTICS: ABC vs A
# ==============================================================================

lnc_abc_vs_a <- lnc_plot_data %>%
  filter(Group %in% c("ABC", "A")) %>%
  group_by(feature_id) %>%
  group_modify(~ {
    x <- .x$log2CPM[.x$Group == "ABC"]
    y <- .x$log2CPM[.x$Group == "A"]
    if(length(unique(c(x, y))) <= 1) {
      tibble(p_value = NA_real_)
    } else {
      test <- wilcox.test(x, y, exact = FALSE)
      tibble(p_value = test$p.value)
    }
  }) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p_value, method = "BH"),
         significance = case_when(FDR < 0.001 ~ "***",
                                  FDR < 0.01 ~ "**",
                                  FDR < 0.05 ~ "*",
                                  is.na(FDR) ~ "NA",
                                  TRUE ~ "ns"))

write.csv(lnc_abc_vs_a, "ncRNA_ABC_vs_A_Wilcoxon_statistics.csv", row.names = FALSE)


# ==============================================================================
# 6. FUNCTION FOR INDIVIDUAL ncRNA PLOTS
# ==============================================================================

make_ncRNA_plot <- function(ncRNA_id) {
  
  dat <- lnc_plot_data %>% filter(feature_id == ncRNA_id)
  summ <- lnc_summary %>% filter(feature_id == ncRNA_id)
  stat <- lnc_abc_vs_a %>% filter(feature_id == ncRNA_id)
  
  max_point <- max(dat$log2CPM, na.rm = TRUE)
  max_errorbar <- max(summ$Mean + ifelse(is.na(summ$SEM), 0, summ$SEM), na.rm = TRUE)
  ymax <- max(max_point, max_errorbar, na.rm = TRUE)
  
  if(!is.finite(ymax) || ymax <= 0.001) {
    bracket_y <- 0.04
    upper_limit <- 0.05
  } else {
    bracket_y <- ymax * 1.08
    upper_limit <- ymax * 1.18
  }
  
  p <- ggplot() +
    
    # Mean bars without borders
    geom_col(data = summ, aes(x = Group, y = Mean, fill = Group),
             width = 0.65, color = NA, inherit.aes = FALSE) +
    
    # SEM
    geom_errorbar(data = summ %>% filter(n > 1),
                  aes(x = Group, ymin = Mean - SEM, ymax = Mean + SEM),
                  width = 0.15, linewidth = 0.9, color = "black",
                  inherit.aes = FALSE) +
    
    # Larger SOLID BLACK individual data points
    geom_jitter(data = dat, aes(x = Group, y = log2CPM),
                width = 0.08, height = 0, size = 3.8,
                shape = 16, color = "black",
                inherit.aes = FALSE) +
    
    scale_fill_manual(values = group_colors, drop = FALSE) +
    
    labs(title = ncRNA_id, x = NULL, y = expression(log[2]*"(CPM + 1)")) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    
    coord_cartesian(ylim = c(0, upper_limit), clip = "off") +
    
    theme_classic(base_size = 16) +
    
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title.y = element_text(face = "bold", color = "black", size = 16),
      axis.text.x = element_text(face = "bold", color = "black", size = 17),
      axis.text.y = element_text(face = "bold", color = "black", size = 16),
      axis.line = element_line(linewidth = 0.9, color = "black"),
      axis.ticks = element_line(linewidth = 0.9, color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
    )
  
  if(nrow(stat) == 1 && !is.na(stat$FDR)) {
    
    annotation_df <- data.frame(group1 = "ABC", group2 = "A",
                                y.position = bracket_y,
                                label = stat$significance)
    
    p <- p + ggpubr::stat_pvalue_manual(annotation_df,
                                        label = "label",
                                        xmin = "group1",
                                        xmax = "group2",
                                        y.position = "y.position",
                                        tip.length = 0.015,
                                        bracket.size = 0.8,
                                        size = 5.5,
                                        inherit.aes = FALSE)
  }
  
  return(p)
}


# ==============================================================================
# 7. GENERATE ncRNA PLOTS
# ==============================================================================

ncRNA_ids <- unique(lnc_plot_data$feature_id)
ncRNA_plots <- lapply(ncRNA_ids, make_ncRNA_plot)
names(ncRNA_plots) <- ncRNA_ids


# ==============================================================================
# 8. SAVE INDIVIDUAL ncRNA PLOTS
# ==============================================================================

dir.create("Individual_ncRNA_expression_plots_ALL_GROUPS", showWarnings = FALSE)

for(i in seq_along(ncRNA_plots)) {
  
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", names(ncRNA_plots)[i])
  
  ggsave(file.path("Individual_ncRNA_expression_plots_ALL_GROUPS",
                   paste0(safe_name, "_expression.png")),
         ncRNA_plots[[i]], width = 5.5, height = 5, dpi = 600)
  
  ggsave(file.path("Individual_ncRNA_expression_plots_ALL_GROUPS",
                   paste0(safe_name, "_expression.svg")),
         ncRNA_plots[[i]], device = svglite::svglite,
         width = 5.5, height = 5)
}


# ==============================================================================
# 9. COMBINE ALL ncRNA PLOTS
# ==============================================================================

combined_ncRNA_plot <- patchwork::wrap_plots(ncRNA_plots, ncol = 3) +
  patchwork::plot_annotation(
    title = "Expression of Candidate ncRNAs",
    subtitle = "Bars represent mean ± SEM; points represent individual samples",
    theme = theme(
      plot.title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 17)
    )
  )

print(combined_ncRNA_plot)

ggsave("All_ncRNA_expression_ALL_5_GROUPS.png", combined_ncRNA_plot, width = 16, height = 14, dpi = 600)
ggsave("All_ncRNA_expression_ALL_5_GROUPS.svg", combined_ncRNA_plot, device = svglite::svglite, width = 16, height = 14)
ggsave("All_ncRNA_expression_ALL_5_GROUPS.pdf", combined_ncRNA_plot, width = 16, height = 14)



# ==============================================================================
# PART B: KEY TARGET GENES
# ==============================================================================

# ==============================================================================
# 10. CONVERT GENE DATA TO LONG FORMAT
# ==============================================================================

gene_long <- gene_raw %>%
  pivot_longer(cols = -gene_name, names_to = "Sample", values_to = "CPM") %>%
  mutate(Group = sub("-[0-9]+$", "", Sample), log2CPM = log2(CPM + 1))

gene_plot_data <- gene_long %>%
  filter(Group %in% plot_groups) %>%
  mutate(Group = factor(Group, levels = plot_groups))

cat("\nSample numbers used for target-gene plots:\n")
print(gene_plot_data %>% distinct(Sample, Group) %>% count(Group))

write.csv(gene_plot_data, "Gene_sample_level_log2CPM_all_groups.csv", row.names = FALSE)


# ==============================================================================
# 11. GENE GROUP MEAN ± SEM
# ==============================================================================

gene_summary <- gene_plot_data %>%
  group_by(gene_name, Group) %>%
  summarise(n = n(),
            Mean = mean(log2CPM, na.rm = TRUE),
            SD = ifelse(n() > 1, sd(log2CPM, na.rm = TRUE), NA_real_),
            SEM = ifelse(n() > 1, SD / sqrt(n()), NA_real_),
            .groups = "drop")

write.csv(gene_summary, "Gene_group_mean_SEM_all_groups.csv", row.names = FALSE)


# ==============================================================================
# 12. GENE STATISTICS: ABC vs A
# ==============================================================================

gene_abc_vs_a <- gene_plot_data %>%
  filter(Group %in% c("ABC", "A")) %>%
  group_by(gene_name) %>%
  group_modify(~ {
    x <- .x$log2CPM[.x$Group == "ABC"]
    y <- .x$log2CPM[.x$Group == "A"]
    if(length(unique(c(x, y))) <= 1) {
      tibble(p_value = NA_real_)
    } else {
      test <- wilcox.test(x, y, exact = FALSE)
      tibble(p_value = test$p.value)
    }
  }) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p_value, method = "BH"),
         significance = case_when(FDR < 0.001 ~ "***",
                                  FDR < 0.01 ~ "**",
                                  FDR < 0.05 ~ "*",
                                  is.na(FDR) ~ "NA",
                                  TRUE ~ "ns"))

write.csv(gene_abc_vs_a, "Gene_ABC_vs_A_Wilcoxon_statistics.csv", row.names = FALSE)


# ==============================================================================
# 13. FUNCTION FOR INDIVIDUAL GENE PLOTS
# ==============================================================================

make_gene_plot <- function(gene_id) {
  
  dat <- gene_plot_data %>% filter(gene_name == gene_id)
  summ <- gene_summary %>% filter(gene_name == gene_id)
  stat <- gene_abc_vs_a %>% filter(gene_name == gene_id)
  
  max_point <- max(dat$log2CPM, na.rm = TRUE)
  max_errorbar <- max(summ$Mean + ifelse(is.na(summ$SEM), 0, summ$SEM), na.rm = TRUE)
  ymax <- max(max_point, max_errorbar, na.rm = TRUE)
  
  if(!is.finite(ymax) || ymax <= 0.001) {
    bracket_y <- 0.04
    upper_limit <- 0.05
  } else {
    bracket_y <- ymax * 1.08
    upper_limit <- ymax * 1.18
  }
  
  p <- ggplot() +
    
    # Mean bars without borders
    geom_col(data = summ, aes(x = Group, y = Mean, fill = Group),
             width = 0.65, color = NA, inherit.aes = FALSE) +
    
    # SEM
    geom_errorbar(data = summ %>% filter(n > 1),
                  aes(x = Group, ymin = Mean - SEM, ymax = Mean + SEM),
                  width = 0.15, linewidth = 0.9, color = "black",
                  inherit.aes = FALSE) +
    
    # Larger SOLID BLACK individual data points
    geom_jitter(data = dat, aes(x = Group, y = log2CPM),
                width = 0.08, height = 0, size = 3.8,
                shape = 16, color = "black",
                inherit.aes = FALSE) +
    
    scale_fill_manual(values = group_colors, drop = FALSE) +
    
    labs(title = gene_id, x = NULL, y = expression(log[2]*"(CPM + 1)")) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    
    coord_cartesian(ylim = c(0, upper_limit), clip = "off") +
    
    theme_classic(base_size = 16) +
    
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title.y = element_text(face = "bold", color = "black", size = 16),
      axis.text.x = element_text(face = "bold", color = "black", size = 17),
      axis.text.y = element_text(face = "bold", color = "black", size = 16),
      axis.line = element_line(linewidth = 0.9, color = "black"),
      axis.ticks = element_line(linewidth = 0.9, color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
    )
  
  if(nrow(stat) == 1 && !is.na(stat$FDR)) {
    
    annotation_df <- data.frame(group1 = "ABC", group2 = "A",
                                y.position = bracket_y,
                                label = stat$significance)
    
    p <- p + ggpubr::stat_pvalue_manual(annotation_df,
                                        label = "label",
                                        xmin = "group1",
                                        xmax = "group2",
                                        y.position = "y.position",
                                        tip.length = 0.015,
                                        bracket.size = 0.8,
                                        size = 5.5,
                                        inherit.aes = FALSE)
  }
  
  return(p)
}


# ==============================================================================
# 14. GENERATE GENE PLOTS
# ==============================================================================

gene_ids <- unique(gene_plot_data$gene_name)
gene_plots <- lapply(gene_ids, make_gene_plot)
names(gene_plots) <- gene_ids


# ==============================================================================
# 15. SAVE INDIVIDUAL GENE PLOTS
# ==============================================================================

dir.create("Individual_key_gene_expression_plots_ALL_GROUPS", showWarnings = FALSE)

for(i in seq_along(gene_plots)) {
  
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", names(gene_plots)[i])
  
  ggsave(file.path("Individual_key_gene_expression_plots_ALL_GROUPS",
                   paste0(safe_name, "_expression.png")),
         gene_plots[[i]], width = 5.5, height = 5, dpi = 600)
  
  ggsave(file.path("Individual_key_gene_expression_plots_ALL_GROUPS",
                   paste0(safe_name, "_expression.svg")),
         gene_plots[[i]], device = svglite::svglite,
         width = 5.5, height = 5)
}


# ==============================================================================
# 16. COMBINE ALL GENE PLOTS
# ==============================================================================

combined_gene_plot <- patchwork::wrap_plots(gene_plots, ncol = 3) +
  patchwork::plot_annotation(
    title = "Expression of Key Target Genes",
    subtitle = "Bars represent mean ± SEM; points represent individual samples",
    theme = theme(
      plot.title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 17)
    )
  )

print(combined_gene_plot)

ggsave("All_key_gene_expression_ALL_5_GROUPS.png", combined_gene_plot, width = 16, height = 13, dpi = 600)
ggsave("All_key_gene_expression_ALL_5_GROUPS.svg", combined_gene_plot, device = svglite::svglite, width = 16, height = 13)
ggsave("All_key_gene_expression_ALL_5_GROUPS.pdf", combined_gene_plot, width = 16, height = 13)


# ==============================================================================
# FINISHED
# ==============================================================================

cat("\n============================================================\n")
cat("ALL EXPRESSION PLOTS COMPLETED\n")
cat("Groups shown: ABC, AB, AC, A, N\n")
cat("Statistical comparison: ABC vs A\n")
cat("Test: Wilcoxon rank-sum + BH-FDR correction\n")
cat("Bars: no border\n")
cat("Individual samples: enlarged solid black dots\n")
cat("Panel titles and axis fonts enlarged\n")
cat("Y-axis automatically adjusted for each feature\n")
cat("============================================================\n")



















