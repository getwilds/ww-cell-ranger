version 1.0

#### STRUCT DEFINITIONS

struct cellRangerSample {
    # sample specific parameters
    String sampleName
    Int? cellNumber
    Array[File] R1Fastqs
    Array[File] R2Fastqs
}

#### WORKFLOW DEFINITION

workflow cell_ranger {
  input {
      Array[cellRangerSample] allSamples
      String cellRangerReferenceTar
  }

  # scatter over each of the samples in the array allSamples
  scatter ( sample in allSamples ) {
      # run cell ranger count on the GEX library 
      call cellRangerCount {
        input:
          R1FastqsGEX = sample.R1Fastqs,
          R2FastqsGEX = sample.R2Fastqs,
          cellRangerReferenceTar = cellRangerReferenceTar,
          cellrangerSample = sample.sampleName,
          threads = 6,
          cellNumber = sample.cellNumber
      }
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] cellrangerbam = cellRangerCount.bam
    Array[File] cellrangerBai = cellRangerCount.bai
    Array[File] cellrangerBarcodes = cellRangerCount.barcodes
    Array[File] cellrangerfilteredBarcodes = cellRangerCount.filteredBarcodes
    Array[File] cellrangerfeatures = cellRangerCount.features
    Array[File] cellrangerMatrix = cellRangerCount.matrix
    Array[Array[File]] cellrangerGlob = cellRangerCount.outputDir
  }
} # End workflow

#### TASK DEFINITIONS

task cellRangerCount {
  input {
    Array[File] R1FastqsGEX
    Array[File] R2FastqsGEX
    File cellRangerReferenceTar
    String cellrangerSample
    Int threads
    Int? cellNumber
  }

  String cellRangerReferenceDir = basename(cellRangerReferenceTar, ".tar.gz")

  command <<<
    set -eo pipefail
    mkdir fastqs
    cp ~{sep=' ' R1FastqsGEX} fastqs
    cp ~{sep=' ' R2FastqsGEX} fastqs
    tar -xvf ~{cellRangerReferenceTar}
    cellranger count \
      --id=countrun \
      --fastqs=fastqs \
      --sample="~{cellrangerSample}" \
      --transcriptome="~{cellRangerReferenceDir}" \
      --localcores=~{threads} \
      ~{"--expect-cells=" + cellNumber}
    mv countrun/outs/possorted_genome_bam.bam "~{cellrangerSample}.possorted_genome_bam.bam"
    mv countrun/outs/possorted_genome_bam.bam.bai "~{cellrangerSample}.possorted_genome_bam.bam.bai"
    mv countrun/outs/raw_feature_bc_matrix/barcodes.tsv.gz "~{cellrangerSample}.barcodes.tsv.gz"
    mv countrun/outs/filtered_feature_bc_matrix/barcodes.tsv.gz "~{cellrangerSample}.filtered.barcodes.tsv.gz"
    mv countrun/outs/filtered_feature_bc_matrix/features.tsv.gz "~{cellrangerSample}.features.tsv.gz"
    mv countrun/outs/filtered_feature_bc_matrix/matrix.mtx.gz "~{cellrangerSample}.matrix.mtx.gz"
  >>>

  output {
    File bam = "~{cellrangerSample}.possorted_genome_bam.bam"
    File bai = "~{cellrangerSample}.possorted_genome_bam.bam.bai"
    File barcodes = "~{cellrangerSample}.barcodes.tsv.gz"
    File filteredBarcodes = "~{cellrangerSample}.filtered.barcodes.tsv.gz"
    File features = "~{cellrangerSample}.features.tsv.gz"
    File matrix = "~{cellrangerSample}.matrix.mtx.gz"
    Array[File] outputDir = glob("./countrun/outs/*")
  }

  runtime {
    docker: "ghcr.io/tefirman/cellranger:latest"
    cpu: threads
  }
}

