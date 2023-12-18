cwlVersion: v1.1
class: CommandLineTool
label: Downstream analysis for RNA and ATAC 
requirements:
  DockerRequirement:
      dockerPull: hubmap/multiome:latest
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
      glob: "multiome_normalized.h5mu"
  mofa_out:
    type: File
    outputBinding:
      glob: "multiome_mofa.hdf5"
  joint_embedding:
    type: File
    outputBinding:
      glob: "leiden_cluster_combined.pdf"
  