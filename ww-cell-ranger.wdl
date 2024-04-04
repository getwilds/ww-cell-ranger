version 1.0
# Generates single cell feature counts via Cell Ranger
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count

#### STRUCT DEFINITIONS

struct CellRangerSample {
    # sample specific parameters
    String sample_name
    Int? cell_number
    Array[File] r1_fastqs
    Array[File] r2_fastqs
}

#### WORKFLOW DEFINITION

workflow cell_ranger {
  input {
      Array[CellRangerSample] all_samples
      File cr_ref_tar
  }

  # scatter over each of the samples in the array all_samples
  scatter (sample in all_samples) {
      # run cell ranger count on the GEX library 
      call CellRangerCount {
        input:
          r1_fastqsGEX = sample.r1_fastqs,
          r2_fastqsGEX = sample.r2_fastqs,
          cr_ref_tar = cr_ref_tar,
          cr_sample = sample.sample_name,
          threads = 6,
          cell_number = sample.cell_number
      }
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] cr_bam = CellRangerCount.bam
    Array[File] cr_bai = CellRangerCount.bai
    Array[File] cr_barcodes = CellRangerCount.barcodes
    Array[File] cr_filtered_barcodes = CellRangerCount.filtered_barcodes
    Array[File] cr_features = CellRangerCount.features
    Array[File] cr_matrix = CellRangerCount.matrix
  }

  parameter_meta {
    all_samples: "list of structs containing details about each sample including names and fastq locations"
    cr_ref_tar: "tar file containing the reference files needed for Cell Ranger"

    cr_bam: "array of aligned bam files for each sample"
    cr_bai: "array of corresponding index files for each aligned bam file"
    cr_barcodes: "array of tsv files containing barcode sequences and GEM wells for each sample"
    cr_filtered_barcodes: "array of tsv files containing a filtered list of barcode sequences for each sample"
    cr_features: "array of tsv files containing a list of features being measured"
    cr_matrix: "array of text files containing Cell Ranger's feature-barcode matrices for each sample"
  }
} # End workflow

#### TASK DEFINITIONS

task CellRangerCount {
  input {
    Array[File] r1_fastqsGEX
    Array[File] r2_fastqsGEX
    File cr_ref_tar
    String cr_sample
    Int threads
    Int? cell_number
  }

  String cr_ref_dir = basename(cr_ref_tar, ".tar.gz")

  command <<<
    set -eo pipefail
    mkdir fastqs
    cp ~{sep=' ' r1_fastqsGEX} fastqs
    cp ~{sep=' ' r2_fastqsGEX} fastqs
    tar -xvf ~{cr_ref_tar}
    cellranger count \
      --id=countrun \
      --fastqs=fastqs \
      --sample="~{cr_sample}" \
      --transcriptome="~{cr_ref_dir}" \
      --localcores=~{threads} \
      ~{"--expect-cells=" + cell_number}
    mv countrun/outs/possorted_genome_bam.bam "~{cr_sample}.possorted_genome_bam.bam"
    mv countrun/outs/possorted_genome_bam.bam.bai "~{cr_sample}.possorted_genome_bam.bam.bai"
    mv countrun/outs/raw_feature_bc_matrix/barcodes.tsv.gz "~{cr_sample}.barcodes.tsv.gz"
    mv countrun/outs/filtered_feature_bc_matrix/barcodes.tsv.gz "~{cr_sample}.filtered.barcodes.tsv.gz"
    mv countrun/outs/filtered_feature_bc_matrix/features.tsv.gz "~{cr_sample}.features.tsv.gz"
    mv countrun/outs/filtered_feature_bc_matrix/matrix.mtx.gz "~{cr_sample}.matrix.mtx.gz"
  >>>

  output {
    File bam = "~{cr_sample}.possorted_genome_bam.bam"
    File bai = "~{cr_sample}.possorted_genome_bam.bam.bai"
    File barcodes = "~{cr_sample}.barcodes.tsv.gz"
    File filtered_barcodes = "~{cr_sample}.filtered.barcodes.tsv.gz"
    File features = "~{cr_sample}.features.tsv.gz"
    File matrix = "~{cr_sample}.matrix.mtx.gz"
  }

  runtime {
    docker: "ghcr.io/getwilds/cellranger:6.0.2"
    cpu: threads
  }

  parameter_meta {
    r1_fastqsGEX: "location of the R1 fastq files for the sample in question"
    r2_fastqsGEX: "location of the R2 fastq files for the sample in question"
    cr_ref_tar: "tar file containing the reference files needed for Cell Ranger"
    cr_sample: "name of the sample in question, used for output file names"
    threads: "number of threads to use when performing Cell Ranger analysis"
    cell_number: "expected number of recovered cells"

    bam: "aligned bam file for the sample in question"
    bai: "corresponding index file for the aligned bam file"
    barcodes: "tsv file containing barcode sequences and GEM wells for the sample in question"
    filtered_barcodes: "tsv file containing a filtered list of barcode sequences for the sample in question"
    features: "tsv file containing a list of features being measured for the sample in question"
    matrix: "text files containing Cell Ranger's feature-barcode matrix"
  }
}

