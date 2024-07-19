cwlVersion: v1.1
class: CommandLineTool
label: Downstream analysis for RNA and ATAC
requirements:
  DockerRequirement:
      dockerPull: hubmap/multiome_analysis:1.1.6
baseCommand: /opt/downstream.py

inputs:
  muon_dir:
    type: File
    inputBinding:
      position: 0
      prefix: "--muon_dir"
outputs:
  muon_processed:
    type: File
    outputBinding:
      glob: "secondary_analysis.h5mu"
  mofa_out:
    type: File
    outputBinding:
      glob: "multiome_mofa.hdf5"
  joint_embedding:
    type: File
    outputBinding:
      glob: "leiden_cluster_combined.pdf"
  rna_embedding:
    type: File
    outputBinding:
      glob: "leiden_cluster_rna.pdf"
  atac_embedding:
    type: File
    outputBinding:
      glob: "leiden_cluster_atac.pdf"

