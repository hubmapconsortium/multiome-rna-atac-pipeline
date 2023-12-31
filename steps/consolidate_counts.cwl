cwlVersion: v1.1
class: CommandLineTool
label: Consolidate RNA, ATAC-seq
requirements:
  DockerRequirement:
      dockerPull: hubmap/multiome_analysis:latest
baseCommand: /opt/consolidate_counts.py

inputs:
  count_matrix_h5ad_rna:
    type: File
    inputBinding:
      position: 0
      prefix: "--rna_file"
  cell_by_bin_matrix_h5ad_atac:
    type: File
    inputBinding:
      position: 1
      prefix: "--atac_cell_by_bin"
  cell_by_gene_matrix_h5ad_atac:
    type: File
    inputBinding:
      position: 2
      prefix: "--atac_cell_by_gene"
  transformation_dir:
    type: Directory?
    inputBinding:
      position: 3
      prefix: "--trans_dir"
  transformation_filename:
    type: string?
    inputBinding:
      position: 4
      prefix: "--trans_filename"
outputs:
  muon_dir:
    type: File
    outputBinding:
      glob: "mudata_raw.h5mu"