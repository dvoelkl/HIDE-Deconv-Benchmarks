# This script reproduces the Morandini benchmark used in
# Dietrich, A., Merotto, L., Pelz, K. et al.
# omnideconv: a unifying framework for using and benchmarking single-cell-informed deconvolution of bulk RNA-seq data.
# Genome Biol 27, 6 (2026).
# https://doi.org/10.1186/s13059-026-03955-w


# Download Link Single Cells (haoSub_sce, deconvData):
# https://figshare.com/ndownloader/files/56063936

# Download Link Bulk (morandini, deconvData):
# https://figshare.com/ndownloader/articles/25347757/versions/3?folder_path=morandini

# The deconvolution is performed with CSsingle following Example 1
# https://wenjshen.github.io/CSsingle/
# Shen, W., Hu, Y., Lei, Y. et al.
# CSsingle: a unified tool for robust decomposition of bulk and spatial transcriptomic data across diverse single-cell references.
# Nucleic Acids Res 54, 8 (2026).
# https://doi.org/10.1093/nar/gkag410


# Assuming the files are in the located in the subfolder data
data_Hao_rds <- "./data/haoSub_sce.rds"
data_morandini_counts <- "./data/morandini_counts.rds"
data_morandini_facs <- "./data/morandini_facs.rds"



if (!"devtools" %in% rownames(installed.packages())) {
  install.packages("devtools")
}

if (!"CSsingle" %in% rownames(installed.packages())) {
  devtools::install_github("wenjshen/CSsingle")
}

library(CSsingle)
library(Seurat)
library(Biobase)
library(nnls)
library(parallel)
library(SingleCellExperiment)


## Data Conversions

# First convert the Hao_sce.rds file to a SeuratObject, as CSsingle requires the counts in the "RNA" assay
print("-> Converting Hao Single Cells")
sce_hao <- as(readRDS(data_Hao_rds), "SingleCellExperiment")

sc.eset <- CreateSeuratObject(
  counts = counts(sce_hao),
  meta.data = as.data.frame(colData(sce_hao))
)


# Second convert the morandini counts to an ExpressionSet
print("-> Convert Morandini counts to ExpressionSet")
morandini_counts <- as.matrix(readRDS(data_morandini_counts))

bulk.eset <- ExpressionSet(assayData = morandini_counts)


# Third load the morandini facs
print("-> Load Morandini facs")
morandini_facs <- as.matrix(readRDS(data_morandini_facs))


## Deconvolution Part

print("-> Running Deconvolution with CSsingle")

cellType <- "cell_type"

# Exclude sparse features and cells and trim bulk and single cells to the same features
prepData <- preprocessData(bulk.eset, sc.eset, cellType)
bulk.eset <- prepData$bulk.eset
sc.eset <- prepData$sc.eset

# Build the differentially expressed matrix
res <- build_diffexpMat(sc.eset, cellType)
ctDEGs <- res$ctDEGs
diffexpMat <- res$diffexpMat
designMat <- res$designMat


res.CSsingle_Comp <- CSsingle_optimize(designMat, diffexpMat, exprs(bulk.eset), ctDEGs, dampened = TRUE, tune.correction = TRUE, n_cores = 1)

C_est <- res.CSsingle_Comp$Est.prop


# Get original names again
celltype_names <- unique(as.character(sce_hao[[cellType]]))
names(celltype_names) <- gsub("[.]", "_", make.names(celltype_names))

colnames(C_est) <- celltype_names[colnames(C_est)]
C_est <- as.data.frame(C_est, check.names = FALSE)


## Aggregate the results according to the computeMetricsNF.R function from deconvBench
# Morandini specific aggregation from computeMetricsNF.R
C_est[["T cells CD4"]] <- C_est[["T cells CD4 conv"]]
C_est[["T cells"]] <- C_est[["T cells CD4"]] + C_est[["T cells CD8"]] + C_est[["Tregs"]]

if ("ILC" %in% colnames(C_est)) {
  C_est[["Lymphocytes"]] <- C_est[["T cells"]] + C_est[["B cells"]] + C_est[["ILC"]] + C_est[["NK cells"]]
} else {
  C_est[["Lymphocytes"]] <- C_est[["T cells"]] + C_est[["B cells"]] + C_est[["NK cells"]]
}

if ("Plasma cells" %in% colnames(C_est)) {
  C_est[["Lymphocytes"]] <- C_est[["Lymphocytes"]] + C_est[["Plasma cells"]]
}


## Calculate metrics according to deconvBench

C_est_long <- data.frame(
  sample = rep(rownames(C_est), times = ncol(C_est)),
  cell_type = rep(colnames(C_est), each = nrow(C_est)),
  estimated_frac = unlist(C_est, use.names = FALSE)
)

morandini_facs_long <- data.frame(
  cell_type = rep(rownames(morandini_facs), times = ncol(morandini_facs)),
  sample = rep(colnames(morandini_facs), each = nrow(morandini_facs)),
  true_frac = c(morandini_facs)
)

results_df <- merge(
  C_est_long,
  morandini_facs_long,
  by = c("sample", "cell_type")
)

results_df$estimated_frac[is.na(results_df$estimated_frac)] <- 0
results_df$true_frac[is.na(results_df$true_frac)] <- 0

x <- as.numeric(results_df$estimated_frac)
y <- as.numeric(results_df$true_frac)

pearson_r <- cor(x, y, method = "pearson")
rmse <- sqrt(mean((x - y)^2))


## Original results from deconvBench
benchmark_results <- data.frame(
  "Method" = c(
    "AutoGeneS",
    "BayesPrism",
    "Bisque",
    "CIBERSORTx",
    "DWLS",
    "MuSiC",
    "Scaden",
    "SCDC",
    "CSsingle"
  ),
  "Pearson R" = c(
    0.819,
    0.254,
    0.872,
    0.648,
    0.909,
    0.319,
    0.795,
    0.468,
    pearson_r
  ),
  "RMSE" = c(
    0.182,
    0.295,
    0.122,
    0.193,
    0.103,
    0.275,
    0.150,
    0.239,
    rmse
  ),
  check.names = FALSE
)

write.csv(benchmark_results, "morandini_results_cssingle.csv", row.names = FALSE, quote = FALSE)