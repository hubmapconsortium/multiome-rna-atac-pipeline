#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
label: multiome pipeline using Salmon and Alevin (HuBMAP scRNA-seq pipeline) and HuBMAP scATAC-seq pipeline
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
inputs:
  fastq_dir_rna:
    label: "Directory containing RNA-seq FASTQ files"
    type: Directory[]
  fastq_dir_atac:
    label: "Directory containing ATAC-seq FASTQ files"
    type: Directory[]
  assay_rna:
    label: "scRNA-seq assay"
    type: string
  assay_atac:
    label: "scATAC-seq assay"
    type: string
  threads_rna:
    label: "Number of threads for Salmon"
    type: int
    default: 1
  threads_atac:
    label: "Number of threads for scATAC-seq"
    type: int?
  expected_cell_count:
    type: int?
  keep_all_barcodes:
    type: boolean?
  trans_dir:
    label: "Directory of barcode transformation mapping file for feature barcoding protocol use TotalSeq B or C"
    type: Directory?
  trans_filename:
    label: "Filename of barcode transformation mapping file for feature barcoding protocol use TotalSeq B or C"
    type: string?
  exclude_bam:
    type: boolean?
outputs:
  muon_original_h5mu:
    outputSource: consolidate_counts/muon_dir
    type: File
    label: "Consolidated expression cell-by-gene, cell-by-bin"
  muon_processed_h5mu:
    outputSource: downstream_analysis/muon_processed
    type: File
    label: "Processed version of raw expression for each modality"
  mofa_model:
    outputSource: downstream_analysis/mofa_out
    type: File
    label: "Multi-omics factor analysis model"
  rna_embedding_result:
    outputSource: downstream_analysis/rna_embedding
    type: File
    label: "Leiden clustering result on rna modality"
  joint_embedding_result:
    outputSource: downstream_analysis/joint_embedding
    type: File
    label: "Leiden clustering result on joint modality"
steps:
  rna_quantification:
    in:
      fastq_dir:
        source: fastq_dir_rna
      assay:
        source: assay_rna
      threads:
        source: threads_rna
      expected_cell_count:
        source: expected_cell_count
      keep_all_barcodes:
        source: keep_all_barcodes
    out:
      - salmon_output
      - count_matrix_h5ad
      - raw_count_matrix
      - genome_build_json
    run: salmon-rnaseq/steps/salmon-quantification.cwl
  atac_quantification:
    in:
      sequence_directory:
        source: fastq_dir_atac
      assay:
        source: assay_atac
      threads:
        source: threads_atac
      exclude_bam:
        source: exclude_bam
    out:
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad
    run: sc-atac-seq-pipeline/sc_atac_seq_prep_process_analyze.cwl
  consolidate_counts:
    in:
      count_matrix_h5ad_rna:
        source:
          rna_quantification/count_matrix_h5ad
      cell_by_bin_matrix_h5ad_atac:
        source:
          atac_quantification/cell_by_bin_h5ad
      cell_by_gene_matrix_h5ad_atac:
        source:
          atac_quantification/cell_by_gene_h5ad
      transformation_dir:
        source:
          trans_dir
      transformation_filename:
        source:
          trans_filename
    out: [muon_dir]
    run: steps/consolidate_counts.cwl
  downstream_analysis:
    in:
      muon_dir:
        source:
          consolidate_counts/muon_dir
    out: 
      - muon_processed
      - mofa_out
      - joint_embedding
    run: steps/downstream.cwl
