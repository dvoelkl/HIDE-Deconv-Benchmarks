# This script reproduces the Morandini benchmark used in
# Dietrich, A., Merotto, L., Pelz, K. et al.
# omnideconv: a unifying framework for using and benchmarking single-cell-informed deconvolution of bulk RNA-seq data.
# Genome Biol 27, 6 (2026).
# https://doi.org/10.1186/s13059-026-03955-w


# Download Link Single Cells (haoSub_sce, deconvData):
# https://figshare.com/ndownloader/files/56063936

# Download Link Bulk (morandini, deconvData):
# https://figshare.com/ndownloader/articles/25347757/versions/3?folder_path=morandini

# This script was run with pip install hide-deconv==0.1.10, thus is using the same model as used in the preprint.
# As v0.1.10 does not provide the lazy deconvolution API, the backend pipelines of the
# command line interface (hide-deconv run) are called directly.

# Assuming the files are in the located in the subfolder data
data_Hao_rds = "./data/haoSub_sce.rds"
data_morandini_counts = "./data/morandini_counts.rds"
data_morandini_facs = "./data/morandini_facs.rds"

# HIDE-Deconv project folder (equivalent to hide-deconv run --path <project_dir>)
hidedeconv_path = "./hidedeconv_morandini_articleversion_wHierarchy"

## Imports
import numpy as np
import pandas as pd
import anndata as ad
import anndata2ri
from scipy.stats import pearsonr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from hide_deconv.config import hidedeconv_config
from hide_deconv.pipelines import (
    init_hidedeconv,
    preprocessing_pipeline,
    train_pipeline,
    deconvolve_hide_pipeline,
)


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
# Define the Morandini cell type hierarchy
adata_hao.obs["morandini_minor"] = adata_hao.obs["cell_type"]

adata_hao.obs.loc[
    adata_hao.obs["cell_type"].isin(
        [
            "T cells CD4 conv",
            "T cells CD8",
            "Tregs",
        ]
    ),
    "morandini_minor",
] = "T cells"

adata_hao.obs["morandini_major"] = adata_hao.obs["morandini_minor"]

adata_hao.obs.loc[
    adata_hao.obs["morandini_minor"].isin(
        [
            "T cells",
            "B cells",
            "ILC",
            "NK cells",
            "Plasma cells",
        ]
    ),
    "morandini_major",
] = "Lymphocytes"


# Run the deconvolution with HIDE-Deconv using the backend of the CLI (hide-deconv run)

print("-> Running Deconvolution with HIDE-Deconv")

# hide-deconv init: create project structure and default config.json
init_hidedeconv(hidedeconv_path)

# hide-deconv setup / config: the CLI reads the single cells and bulk from file
sc_file_name = hidedeconv_path + "/data/haoSub.h5ad"
bulk_file_name = hidedeconv_path + "/data/morandini_counts.csv"
ad.settings.allow_write_nullable_strings = True  # pandas string arrays from anndata2ri
adata_hao.write_h5ad(sc_file_name)
df_morandini_counts.to_csv(bulk_file_name)

hconf = hidedeconv_config.load(hidedeconv_path + "/config.json")
hconf.sc_file_name = sc_file_name
hconf.bulk_file_name = bulk_file_name
hconf.sub_ct_col = "cell_type"  # finest cell type annotation level
hconf.higher_ct_cols = ["morandini_minor", "morandini_major"]  # higher levels
hconf.n_genes = 5000
hconf.n_train_bulks = 10000
hconf.n_cells_per_bulk = 100
hconf.n_hide_iter = 1000
hconf.save(hidedeconv_path + "/config.json")

# hide-deconv preprocess (with domain transfer)
preprocessing_pipeline(hidedeconv_path, True)
hconf = hidedeconv_config.load(hidedeconv_path + "/config.json")
hconf.preprocessed = True
hconf.domainTransfer = True
hconf.save(hidedeconv_path + "/config.json")

# hide-deconv train
train_pipeline(hidedeconv_path)

# hide-deconv deconv
results = deconvolve_hide_pipeline(hidedeconv_path, None)


# HIDE-deconv results have form cell_types x samples

C_est = results[0].T
C_est = C_est.rename(
    columns={
        "T cells CD4 conv": "T cells CD4",
    }
)

C_est["T cells"] = results[1].T["T cells"]
C_est["Lymphocytes"] = results[2].T["Lymphocytes"]


## Calculate metrics according to deconvBench

C_est_long = (
    C_est.rename_axis("sample")
    .reset_index()
    .melt(
        id_vars="sample",
        var_name="cell_type",
        value_name="estimated_frac",
    )
)

df_morandini_facs_long = (
    df_morandini_facs.rename_axis("cell_type")
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
benchmark_results = pd.DataFrame(
    {
        "Method": [
            "AutoGeneS",
            "BayesPrism",
            "Bisque",
            "CIBERSORTx",
            "DWLS",
            "MuSiC",
            "Scaden",
            "SCDC",
            "HIDE-Deconv (v0.1.10, w. Hierarchy)",
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
    }
)

benchmark_results.to_csv("morandini_results_articleVersion_wHierarchy.csv", index=False)