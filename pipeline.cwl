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
     label: "RNA-seq assay"
     type: string
  assay_atac:
     label: "ATAC-seq assay"
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
  exclude_bam:
    type: boolean?
  atac_metadata_file:
    type: File?
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
  atac_embedding_result:
    outputSource: downstream_analysis/atac_embedding
    type: File
    label: "Leiden clustering result on atac modality"
  joint_embedding_result:
    outputSource: downstream_analysis/joint_embedding
    type: File
    label: "Leiden clustering result on joint modality"
  scanpy_qc_results:
    outputSource: rna_qc/scanpy_qc_results
    type: File
    label: "Quality control metrics from Scanpy"
  rna_qc_report:
    outputSource: rna_qc/qc_metrics
    type: File
    label: "Quality control report in JSON format"
  atac_qc_report:
    outputSource: atac_qc/qc_report
    type: File
    label "Quality control report in JSON format"
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
      metadata_file:
        source: atac_metadata_file
    out:
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad
      - genome_build_json
      - bam_file
      - bam_index
      - image_file
      - archR_project
    run: sc-atac-seq-pipeline/steps/sc_atac_seq_prep_process_init.cwl

  analyze_with_ArchR:
    run: sc-atac-seq-pipeline/steps/sc_atac_seq_analyze_steps/archr_clustering.cwl
    in:
      image_file: sc_atac_seq_prep_process_init/image_file
      archr_project: sc_atac_seq_prep_process_init/archr_project
    out:
      - peaks_bed

  atac_qc:
    run: sc-atac-seq-pipeline/steps/qc_measures.cwl
    in:
      bam_file: atac_quantification/bam_file
      bam_index: atac_quantification/bam_index
      peak_file: analyze_with_ArchR/peaks_bed
      cell_by_bin_h5ad: atac_quantification/cell_by_bin_h5ad
    out:
      - qc_report

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
      rna_genome_build_path:
        source:
          rna_quantification/genome_build_json
      atac_genome_build_path:
        source:
          atac_quantification/genome_build_json
      assay_atac:
        source:
          assay_atac
      atac_metadata_file:
        source:
          atac_metadata_file
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
      - rna_embedding
      - atac_embedding
    run: steps/downstream.cwl

  rna_qc:
    in:
      assay:
        source: assay_rna
      primary_matrix_path:
        source: rna_quantification/count_matrix_h5ad
      secondary_matrix_path:
        source: downstream_analysis/muon_processed
      salmon_dir:
        source: rna_quantification/salmon_output
    out:
      - scanpy_qc_results
      - qc_metrics
    run: salmon-rnaseq/steps/compute-qc-metrics.cwl
    label: "Compute QC metrics"