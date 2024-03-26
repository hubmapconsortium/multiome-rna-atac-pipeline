#!/usr/bin/env python3
from argparse import ArgumentParser
from os.path import exists
from pathlib import Path

import muon as mu
import pandas as pd
import scanpy as sc


def generate_barcode_dict():
    #colnames = ["transformed", "original"]
    #trans_info = pd.read_csv(directory, sep="\t", header=None, names=colnames)

    atac_barcodes_file = "/opt/atac_barcodes.txt"
    rna_barcodes_file = "/opt/cellranger_barcodes.txt"
    atac_barcodes = pd.read_csv(atac_barcodes_file, sep="\n", header=None)
    rna_barcodes = pd.read_csv(rna_barcodes_file, sep="\n", header=None)

    barcode_dict = dict(zip(atac_barcodes, rna_barcodes))
    print("Generated ", len(barcode_dict), " key-value pairs for RNA barcode transformation.")
    return barcode_dict


def transform_rna_barcode(barcode_dict: Path, rna_name: list):
    count = 0
    for i in range(len(rna_name)):
        barcode = rna_name[i]
        if barcode not in barcode_dict:
            count += 1
            continue
        else:
            rna_name[i] = barcode_dict[barcode]
    return rna_name, count


def main(
    rna_file: Path,
    atac_cell_by_bin: Path,
    atac_cell_by_gene: Path,
    assay: str,
):
    rna_expr = mu.read(str(rna_file))
    rna_name = list(rna_expr.obs_names)
    cbb = mu.read(str(atac_cell_by_bin))
    cbg = mu.read(str(atac_cell_by_gene))

    #get rid of the BAM_data# artifact from ArchR in the atac-seq data
    cbg_names = list(cbg.obs_names)
    cbg_names = [s.replace("BAM_data#",'') for s in cbg_names]
    cbg.obs.index = cbg_names

    cbb_names = list(cbb.obs_names)
    cbb_names = [s.replace("BAM_data#",'') for s in cbb_names]
    cbb.obs.index = cbb_names
    
    if(assay=="multiome_10x"):
        print("Performing transformation step of the cellular barcodes of RNA...")
        barcode_dict = generate_barcode_dict(trans_file_path)
        rna_expr.obs.index, count = transform_rna_barcode(barcode_dict, rna_name)
        rna_name = list(rna_expr.obs_names)
        print(
            "A total of",
            count,
            "out of",
            len(rna_name),
            "RNA cell barcodes were not found in RNA barcode transformation dictionary.",
        )

    # print("There are", len(common_cells), "common cells in RNA and Atac experiments.")
    mdata = mu.MuData({"rna": rna_expr, "atac_cell_by_bin": cbb, "atac_cell_by_gene": cbg})
    mu.pp.intersect_obs(mdata)
    print(
        "Saving MuData filtered by",
        len(mdata),
        "common cells in RNA and Atac experiments...",
    )

    mdata.write("mudata_raw.h5mu")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--rna_file", type=Path)
    p.add_argument("--atac_cell_by_bin", type=Path)
    p.add_argument("--atac_cell_by_gene", type=Path)
    p.add_argument("--assay", type=str)

    args = p.parse_args()

    main(
        args.rna_file,
        args.atac_cell_by_bin,
        args.atac_cell_by_gene,
        args.assay,
    )
