
# ============================================================
# Alzheimer's Disease Gene Expression Analysis
# Dataset: GSE5281_EC (Entorhinal Cortex)
# ============================================================

# STEP 1: LOAD PACKAGES
library(KEGGdzPathwaysGEO)
library(Biobase)

# STEP 2: LOAD DATASET
data("GSE5281_EC", package = "KEGGdzPathwaysGEO")

eset <- GSE5281_EC
expr <- exprs(eset)

# STEP 3: DEFINE SAMPLE GROUPS
group <- factor(
  tolower(as.character(pData(eset)$Group)),
  levels = c("c", "d"),
  labels = c("Control", "Alzheimer")
)

# STEP 4: CHECK DATASET
cat("Total probes:", nrow(expr), "\n")
cat("Total samples:", ncol(expr), "\n")

print(table(group))
cat("Platform:", annotation(eset), "\n")

# STEP 5: CHECK MISSING VALUES
cat("Missing values:", sum(is.na(expr)), "\n")
cat("Infinite values:", sum(is.infinite(expr)), "\n")

# STEP 6: CREATE OUTPUT FOLDERS
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# STEP 7: GENE EXPRESSION BOXPLOT
png(
  "figures/01_Expression_Boxplot.png",
  width = 1100,
  height = 700
)

boxplot(
  expr,
  outline = FALSE,
  las = 2,
  col = "lightblue",
  main = "Gene Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Expression Values"
)

dev.off()



# Load packages for plotting and heatmap visualization
library(ggplot2)
library(pheatmap)


# ------------------------------------------------------------
#  PRINCIPAL COMPONENT ANALYSIS (PCA)
# ------------------------------------------------------------

# Calculate the variance of each gene probe across all samples
gene_variance <- apply(expr, 1, var)

# Select the 500 probes with the highest variance
top_genes <- order(gene_variance, decreasing = TRUE)[1:500]

# Perform PCA to identify patterns between samples
# Transpose the data so rows represent samples
pca <- prcomp(
  t(expr[top_genes, ]),
  center = TRUE,
  scale. = FALSE
)

# Calculate the percentage of variance explained by each PC
variance_percent <- round(
  100 * pca$sdev^2 / sum(pca$sdev^2), 2
)

# Create a data frame containing PCA coordinates and groups
pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = group
)

# Create a scatter plot to compare Control and Alzheimer samples
pca_plot <- ggplot(
  pca_data, aes(PC1, PC2, color = Group)
) +
  geom_point(size = 3) +
  labs(
    title = "PCA: Control vs Alzheimer's Disease",
    x = paste0("PC1 (", variance_percent[1], "%)"),
    y = paste0("PC2 (", variance_percent[2], "%)")
  ) +
  theme_minimal()

# Save the PCA plot as a PNG image
ggsave(
  "figures/02_PCA_Plot.png",
  plot = pca_plot,
  width = 8,
  height = 6
)

# Calculate Euclidean distances between all sample pairs
sample_distance <- as.matrix(dist(t(expr)))

# Assign sample names to the distance matrix
rownames(sample_distance) <- colnames(expr)
colnames(sample_distance) <- colnames(expr)

# Create group labels for heatmap annotations
sample_annotation <- data.frame(Group = group)
rownames(sample_annotation) <- colnames(expr)


# ------------------------------------------------------------
# SAMPLE DISTANCE HEATMAP
# ------------------------------------------------------------

# Create sample distance heatmap
heatmap_plot <- pheatmap(
  sample_distance,
  annotation_row = sample_annotation,
  annotation_col = sample_annotation,
  main = "Sample-to-Sample Distance Heatmap",
  border_color = NA,
  silent = TRUE
)

# Open PNG graphics device
png(
  "figures/03_Sample_Distance_Heatmap.png",
  width = 1100,
  height = 1000
)

# Explicitly draw the heatmap using grid graphics
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

# Close graphics device
dev.off()



# ============================================================
# DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ============================================================

# Load packages for statistical analysis and gene annotation
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)

# Create design matrix for Control and Alzheimer groups
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Define comparison: Alzheimer vs Control
contrast_matrix <- makeContrasts(
  Alzheimer - Control,
  levels = design
)

# Fit linear model and apply empirical Bayes statistics
fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)

# Extract differential expression results
# Adjust p-values using the Benjamini-Hochberg method
results <- topTable(
  fit,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Store probe identifiers as a separate column
results$Probe_ID <- rownames(results)

# Map probe IDs to gene symbols
results$Gene_Symbol <- mapIds(
  hgu133plus2.db,
  keys = results$Probe_ID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Add full gene names using the annotation database
results$Gene_Name <- mapIds(
  hgu133plus2.db,
  keys = results$Probe_ID,
  column = "GENENAME",
  keytype = "PROBEID",
  multiVals = "first"
)

# Select significant probes using adjusted p-value < 0.05
significant_results <- results[
  !is.na(results$adj.P.Val) &
    results$adj.P.Val < 0.05,
]

# Save complete and significant results as CSV files
write.csv(
  results,
  "results/All_Results.csv",
  row.names = FALSE
)

write.csv(
  significant_results,
  "results/Significant_Probes.csv",
  row.names = FALSE
)

# Display summary of differential expression results
cat("Total probes:", nrow(results), "\n")
cat("Significant probes:", nrow(significant_results), "\n")

# Count upregulated and downregulated significant probes
cat("Upregulated:",
    sum(significant_results$logFC > 0), "\n")
cat("Downregulated:",
    sum(significant_results$logFC < 0), "\n")

# Display the top 10 probes with their statistical results
print(
  head(
    results[, c(
      "Probe_ID", "Gene_Symbol",
      "logFC", "adj.P.Val"
    )],
    10
  )
)



# ============================================================
# DIFFERENTIAL EXPRESSION VISUALIZATION
# ============================================================

# Classify probes based on adjusted p-value and logFC
# Positive logFC = Upregulated, Negative logFC = Downregulated
results$Regulation <- ifelse(
  results$adj.P.Val < 0.05 & results$logFC > 0,
  "Upregulated",
  ifelse(
    results$adj.P.Val < 0.05 & results$logFC < 0,
    "Downregulated",
    "Not Significant"
  )
)

# Create a volcano plot to visualize differential expression
# X-axis shows log2 fold change
# Y-axis shows statistical significance (-log10 adjusted p-value)
volcano_plot <- ggplot(
  results,
  aes(
    x = logFC,
    y = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    color = Regulation
  )
) +
  geom_point(alpha = 0.5, size = 1) +
  
  # Assign colors based on probe regulation
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey"
    )
  ) +
  
  # Add a dashed line representing the significance threshold
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot: Alzheimer's vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

# Save the volcano plot as a PNG image
ggsave(
  "figures/04_Volcano_Plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6
)

# Select the top 30 statistically significant probes
top30 <- head(significant_results$Probe_ID, 30)

# Extract expression values for the selected probes
heatmap_expr <- expr[top30, , drop = FALSE]

# Retrieve gene symbols for the selected probes
top30_symbols <- results$Gene_Symbol[
  match(top30, results$Probe_ID)
]

# Use gene symbols as labels, with probe IDs for identification
# If a gene symbol is missing, use the probe ID instead
top30_labels <- ifelse(
  is.na(top30_symbols) | top30_symbols == "",
  top30,
  paste0(top30_symbols, " (", top30, ")")
)

rownames(heatmap_expr) <- top30_labels

# Generate a heatmap of the top 30 significant probes
png(
  "figures/05_Top30_Heatmap.png",
  width = 1100,
  height = 1000
)

# Standardize expression values by row for comparison
# Add sample group labels to identify Control and Alzheimer samples
pheatmap(
  heatmap_expr,
  scale = "row",
  annotation_col = sample_annotation,
  main = "Top 30 Differentially Expressed Probes",
  border_color = NA
)

# Close the graphics device and save the heatmap
dev.off()


expected_files <- c(
  "figures/01_Expression_Boxplot.png",
  "figures/02_PCA_Plot.png",
  "figures/03_Sample_Distance_Heatmap.png",
  "figures/04_Volcano_Plot.png",
  "figures/05_Top30_Heatmap.png",
  "results/All_Results.csv",
  "results/Significant_Probes.csv"
)

print(data.frame(
  File = expected_files,
  Exists = file.exists(expected_files)
))

stopifnot(all(file.exists(expected_files)))

cat("\nAll analysis outputs verified successfully.\n")





