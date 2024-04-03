#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from pathlib import Path

import muon as mu
from anndata import AnnData

from utils import Assay

DATA_DIR = Path("/data")


def generate_barcode_dict() -> dict[str, str]:
    with open(DATA_DIR / "atac_barcodes_rev.txt") as f:
        atac_barcodes = [line.strip() for line in f]
    with open(DATA_DIR / "cellranger_barcodes.txt") as f:
        rna_barcodes = [line.strip() for line in f]
    assert len(atac_barcodes) == len(rna_barcodes)

    barcode_dict = dict(zip(atac_barcodes, rna_barcodes))
    print("Read", len(barcode_dict), "key-value pairs for ATAC-seq barcode transformation.")
    return barcode_dict


def drop_bam_data_obs_index_prefix(adata: AnnData):
    """
    Operates on adata in place. This is a separate function to reduce the
    possibility of copy/paste problems when applying this to multiple
    AnnData objects.

    :param adata:
    """
    adj_names = [s.removeprefix("BAM_data#") for s in adata.obs.index]
    adata.obs.index = adj_names


def map_atac_barcodes(adata: AnnData, barcode_mapping: dict[str, str]) -> AnnData:
    """
    This is a separate function to reduce the possibility of copy/paste
    problems when applying this to multiple AnnData objects. Does not
    operate in place since we might return a subset of `adata` if not
    all barcodes are present in `barcode_mapping`.

    :param adata:
    :param barcode_mapping:
    """
    selection = [b in barcode_mapping for b in adata.obs.index]
    adata = adata[selection, :]
    mapped_barcodes = [barcode_mapping[b] for b in adata.obs.index]
    adata.obs.index = mapped_barcodes
    return adata.copy()


def main(
    rna_file: Path,
    atac_cell_by_bin: Path,
    atac_cell_by_gene: Path,
    rna_genome_build_path: Path,
    atac_genome_build_path: Path,
    assay_atac: Assay,
):
    rna_expr = mu.read(str(rna_file))
    cbb = mu.read(str(atac_cell_by_bin))
    cbg = mu.read(str(atac_cell_by_gene))

    with open(rna_genome_build_path) as f:
        rna_genome_build_info = json.load(f)

    with open(atac_genome_build_path) as f:
        atac_genome_build_info = json.load(f)

    rna_expr.uns["genome_build"] = rna_genome_build_info
    cbg.uns["genome_build"] = atac_genome_build_info
    cbb.uns["genome_build"] = atac_genome_build_info

    # get rid of the BAM_data# artifact from ArchR in the atac-seq data
    drop_bam_data_obs_index_prefix(cbg)
    drop_bam_data_obs_index_prefix(cbb)

    if assay_atac == Assay.MULTIOME_10X:
        print("Performing transformation step of the cellular barcodes of ATAC-seq")
        barcode_dict = generate_barcode_dict()
        cbg = map_atac_barcodes(cbg, barcode_dict)
        cbb = map_atac_barcodes(cbb, barcode_dict)

    # print("There are", len(common_cells), "common cells in RNA and Atac experiments.")
    mdata = mu.MuData({"rna": rna_expr, "atac_cell_by_bin": cbb, "atac_cell_by_gene": cbg})
    mu.pp.intersect_obs(mdata)
    print(
        "Saving MuData filtered by",
        len(mdata),
        "common cells in RNA and ATAC experiments",
    )

    mdata.write("mudata_raw.h5mu")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--rna_file", type=Path)
    p.add_argument("--atac_cell_by_bin", type=Path)
    p.add_argument("--atac_cell_by_gene", type=Path)
    p.add_argument("--rna_genome_build_path", type=Path)
    p.add_argument("--atac_genome_build_path", type=Path)
    p.add_argument("--assay_atac", choices=list(Assay), type=Assay)

    args = p.parse_args()

    main(
        args.rna_file,
        args.atac_cell_by_bin,
        args.atac_cell_by_gene,
        args.rna_genome_build_path,
        args.atac_genome_build_path,
        args.assay_atac,
    )
