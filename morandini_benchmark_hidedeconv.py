# This script reproduces the Morandini benchmark used in
# Dietrich, A., Merotto, L., Pelz, K. et al. 
# omnideconv: a unifying framework for using and benchmarking single-cell-informed deconvolution of bulk RNA-seq data. 
# Genome Biol 27, 6 (2026). 
# https://doi.org/10.1186/s13059-026-03955-w


# Download Link Single Cells (haoSub_sce, deconvData): 
# https://figshare.com/ndownloader/files/56063936

# Download Link Bulk (morandini, deconvData):
# https://figshare.com/ndownloader/articles/25347757/versions/3?folder_path=morandini


# Assuming the files are in the located in the subfolder data
data_Hao_rds = './data/haoSub_sce.rds'
data_morandini_counts = './data/morandini_counts.rds'
data_morandini_facs = './data/morandini_facs.rds'

## Imports
import scanpy as sc
import numpy as np
import pandas as pd
import anndata2ri
import anndata as ad
from scipy.stats import pearsonr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from hide_deconv import deconvolution


## Data Conversions

# First convert the Hao_sce.rds file to AnnData single cell
print("-> Converting Hao Single Cells")
r_hao = ro.r["readRDS"](data_Hao_rds)

ro.r(
    """
    library(SingleCellExperiment)
    """
)

r_hao = ro.r["as"](r_hao, "SingleCellExperiment")

with localconverter(anndata2ri.converter):
    adata_hao = ro.conversion.rpy2py(r_hao)

# Second convert the morandini counts to pandas
print("-> Convert Morandini counts to pandas")
r_morandini_counts = ro.r["readRDS"](data_morandini_counts)
with localconverter(pandas2ri.converter): 
    df_morandini_counts = ro.conversion.rpy2py(r_morandini_counts)

df_morandini_counts = pd.DataFrame(df_morandini_counts)

# Set indexes
r_morandini_counts_rownames = ro.r["rownames"](r_morandini_counts)

with localconverter(pandas2ri.converter): 
    index = ro.conversion.rpy2py(r_morandini_counts_rownames) 
    df_morandini_counts.index = pd.Index(index)

# Set column names
r_morandini_counts_colnames = ro.r["colnames"](r_morandini_counts)
with localconverter(pandas2ri.converter): 
    colnames = ro.conversion.rpy2py(r_morandini_counts_colnames) 
    df_morandini_counts.columns = pd.Index(colnames)


# Third convert the morandini facs to pandas
print("-> Convert Morandini facs to pandas")
r_morandini_facs = ro.r["readRDS"](data_morandini_facs)
with localconverter(pandas2ri.converter): 
    df_morandini_facs = ro.conversion.rpy2py(r_morandini_facs)

df_morandini_facs = pd.DataFrame(df_morandini_facs)

# Set indexes
r_morandini_facs_rownames = ro.r["rownames"](r_morandini_facs)

with localconverter(pandas2ri.converter): 
    index = ro.conversion.rpy2py(r_morandini_facs_rownames) 
    df_morandini_facs.index = pd.Index(index)

# Set column names
r_morandini_facs_colnames = ro.r["colnames"](r_morandini_facs)
with localconverter(pandas2ri.converter): 
    colnames = ro.conversion.rpy2py(r_morandini_facs_colnames) 
    df_morandini_facs.columns = pd.Index(colnames)


## Deconvolution Part
# Run the deconvolution with HIDE-Deconv with standard parameters
print("-> Running Deconvolution with HIDE-Deconv")
results = deconvolution(adata_hao, df_morandini_counts)

C_est = results[0].T # Transpose, as HIDE-deconv results have form cell_types x samples

## Aggregate the results according to the computeMetricsNF.R function from deconvBench
# Morandini specific aggregation from computeMetricsNF.R
C_est["T cells CD4"] = C_est["T cells CD4 conv"]
C_est["T cells"] = (
    C_est["T cells CD4"]
    + C_est["T cells CD8"]
    + C_est["Tregs"]
)

if "ILC" in C_est.columns:
    C_est["Lymphocytes"] = (
        C_est["T cells"]
        + C_est["B cells"]
        + C_est["ILC"]
        + C_est["NK cells"]
    )
else:
    C_est["Lymphocytes"] = (
        C_est["T cells"]
        + C_est["B cells"]
        + C_est["NK cells"]
    )

if "Plasma cells" in C_est.columns:
    C_est["Lymphocytes"] = (
        C_est["Lymphocytes"]
        + C_est["Plasma cells"]
    )


## Calculate metrics according to deconvBench

C_est_long = (
    C_est
    .rename_axis("sample")
    .reset_index()
    .melt(
        id_vars="sample",
        var_name="cell_type",
        value_name="estimated_frac",
    )
)

df_morandini_facs_long = (
    df_morandini_facs
    .rename_axis("cell_type")
    .reset_index()
    .melt(
        id_vars="cell_type",
        var_name="sample",
        value_name="true_frac",
    )
)

results_df = C_est_long.merge(
    df_morandini_facs_long,
    on=["sample", "cell_type"],
    how="inner",
)

results_df["estimated_frac"] = results_df["estimated_frac"].fillna(0)
results_df["true_frac"] = results_df["true_frac"].fillna(0)

x = results_df["estimated_frac"].to_numpy(dtype=float)
y = results_df["true_frac"].to_numpy(dtype=float)

pearson_r = pearsonr(x, y).statistic
rmse = np.sqrt(np.mean((x - y) ** 2))

## Original results from deconvBench
benchmark_results = pd.DataFrame({
    "Method": [
        "AutoGeneS",
        "BayesPrism",
        "Bisque",
        "CIBERSORTx",
        "DWLS",
        "MuSiC",
        "Scaden",
        "SCDC",
        "HIDE-Deconv",
    ],
    "Pearson R": [
        0.819,
        0.254,
        0.872,
        0.648,
        0.909,
        0.319,
        0.795,
        0.468,
        pearson_r,
    ],
    "RMSE": [
        0.182,
        0.295,
        0.122,
        0.193,
        0.103,
        0.275,
        0.150,
        0.239,
        rmse,
    ],
})

benchmark_results.to_csv("morandini_results.csv", index=False)
