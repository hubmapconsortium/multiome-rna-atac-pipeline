#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import scanpy as sc
from muon import prot as pt

from plot_utils import new_plot


def main(mudata_raw: Path):
    expr = mu.read(str(mudata_raw))
    expr.obs["num_genes_rna"] = (expr["rna"].X > 0).sum(axis=1)
    mu.pp.filter_obs(expr, "num_genes_rna", lambda x: x > 200)
    rna_expr = expr["rna"]
    rna_expr.X = rna_expr.layers["spliced"]
    print(rna_expr)
    atac_cbg_expr = expr["atac_cell_by_gene"]
    atac_cbb_expr = expr["atac_cell_by_bin"]
    print("Construct muon object with RNA and ATAC cell-by-gene data.")
    mdata_raw = mu.MuData({"rna": rna_expr, "atac_cbg": atac_cbg_expr})
    print(mdata_raw)

    mdata_rna = mdata_raw.mod["rna"]
    mdata_ataccbg = mdata_raw.mod["atac_cbg"]

    ## RNA normalization
    # Can filter out cells that do not pass QC
    print("Normalizing RNA data...")
    sc.pp.normalize_total(mdata_rna, target_sum=1e4)
    sc.pp.log1p(mdata_rna)

    # ATAC TFIDF normalization
    mu.atac.pp.tfidf(mdata_ataccbg, scale_factor=1e4)

    # find highly variable genes for ATAC
    sc.pp.highly_variable_genes(mdata_ataccbg, min_mean=0.02, max_mean=4, min_disp=0.5)
    print("Found", np.sum(mdata_ataccbg.var.highly_variable), "highly variable genes in ATAC-seq.")

    ## Downstream analysis for Atac
    print("Performing downstream analysis for ATAC...")
    sc.tl.pca(mdata_ataccbg)
    sc.pp.neighbors(mdata_ataccbg)

    for axis in [0, 1]:
        neighbor_counts = np.array((mdata_ataccbg.obsp["distances"] > 0).sum(axis=axis)).flatten()
        mdata_ataccbg.obs.loc[:, f"neighbor_counts_ax{axis}"] = neighbor_counts
        cells_with_neighbors = (neighbor_counts > 0).astype(int)
        mdata_ataccbg.obs.loc[:, f"has_neighbors_ax{axis}"] = cells_with_neighbors
    mdata_ataccbg.obs.loc[:, "has_neighbors"] = mdata_ataccbg.obs.loc[:, "has_neighbors_ax1"]
    mu.pp.filter_obs(mdata_ataccbg, "has_neighbors", lambda x: x > 0)

    sc.tl.leiden(mdata_ataccbg)
    sc.tl.umap(mdata_ataccbg)
    with new_plot():
        sc.pl.umap(mdata_ataccbg, color="leiden", legend_loc="on data")
        plt.savefig("leiden_cluster_atac.pdf")

    ## Downstream analysis for RNA
    print("Performing downstream analysis for RNA...")

    # find highly variable genes for RNA-seq
    sc.pp.highly_variable_genes(mdata_rna, min_mean=0.02, max_mean=4, min_disp=0.5)
    print("Found", np.sum(mdata_rna.var.highly_variable), "highly variable genes in RNA-seq.")

    # scaling the data
    # save log-normalised counts
    mdata_rna.raw = mdata_rna
    print("Scaling the log-normalised counts to zero mean and unit variance...")
    sc.pp.scale(mdata_rna, max_value=10)

    # pca and neighbourhood graph
    print("Performing PCA...")
    sc.tl.pca(mdata_rna)
    print("Constructing neighborhood graph...")
    sc.pp.neighbors(mdata_rna)

    for axis in [0, 1]:
        neighbor_counts = np.array((mdata_rna.obsp["distances"] > 0).sum(axis=axis)).flatten()
        mdata_rna.obs.loc[:, f"neighbor_counts_ax{axis}"] = neighbor_counts
        cells_with_neighbors = (neighbor_counts > 0).astype(int)
        mdata_rna.obs.loc[:, f"has_neighbors_ax{axis}"] = cells_with_neighbors
    mdata_rna.obs.loc[:, "has_neighbors"] = mdata_rna.obs.loc[:, "has_neighbors_ax1"]
    mu.pp.filter_obs(mdata_rna, "has_neighbors", lambda x: x > 0)

    # clustering
    print("Performing leiden clustering with the computed neighbourhood graph...")
    sc.tl.leiden(mdata_rna)
    sc.tl.umap(mdata_rna)
    with new_plot():
        sc.pl.umap(mdata_rna, color="leiden", legend_loc="on data")
        plt.savefig("leiden_cluster_rna.pdf")

    ## Multi-omics integration
    mdata_raw.update()
    # # if we filter the cells at the RNA QC step, subset them in the protein modality
    mu.pp.intersect_obs(mdata_raw)
    print(mdata_raw)
    mdata_raw.write("multiome_normalized.h5mu")
    ## Multi-omics factor analysis
    mu.tl.mofa(mdata_raw, outfile="multiome_mofa.hdf5", n_factors=30)

    # multiplex clustering
    sc.pp.neighbors(mdata_raw["rna"])
    sc.pp.neighbors(mdata_raw["atac_cbg"])
    print(mdata_raw)
    sc.pp.neighbors(mdata_raw, use_rep="X_mofa", key_added="mofa")
    sc.tl.umap(mdata_raw, neighbors_key="mofa")
    sc.tl.leiden(mdata_raw, resolution=1.0, neighbors_key="mofa", key_added="leiden_wnn")
    # print(mdata_raw)
    with new_plot():
        sc.pl.umap(mdata_raw, color="leiden_wnn", legend_loc="on data")
        plt.savefig("leiden_cluster_combined.pdf")

    # add the cellbybin data back for output
    atac_cbb_expr.var["highly_variable"] = False
    mdata_raw.mod["atac_cbb"] = atac_cbb_expr

    mdata_raw.write("secondary_analysis.h5mu")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--mudata_raw", type=Path)
    args = p.parse_args()

    main(args.mudata_raw)
