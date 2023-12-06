#!/usr/bin/env python3
from argparse import ArgumentParser
from os.path import exists
from pathlib import Path

import muon as mu
import pandas as pd
import scanpy as sc


def generate_barcode_dict(directory: Path):
    colnames = ["transformed", "original"]
    trans_info = pd.read_csv(directory, sep="\t", header=None, names=colnames)
    barcode_dict = dict(zip(trans_info.transformed, trans_info.original))
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
    trans_dir: Path,
    trans_filename: str,
):

    rna_expr = mu.read(str(rna_file))
    rna_name = list(rna_expr.obs_names)
    cbb = mu.read(str(atac_cell_by_bin))
    cbg = mu.read(str(atac_cell_by_gene))
    
    # if the transformation file of given name exist, perform transformation step
    if trans_dir != None:
        trans_file_path = trans_dir / trans_filename
        if exists(trans_file_path):
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
        else:
            raise ValueError(trans_filename, " is not found under given directory.")

    #print("There are", len(common_cells), "common cells in RNA and Atac experiments.")
    mdata = mu.MuData({"rna": rna_expr, "atac_cell_be_bin": cbb, "atac_cell_by_gene": cbg})
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
    p.add_argument("--trans_dir", type=Path)
    p.add_argument("--trans_filename", type=str)

    args = p.parse_args()

    main(args.rna_file, args.atac_cell_by_bin, args.atac_cell_by_gene, args.trans_dir, args.trans_filename)
