
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

