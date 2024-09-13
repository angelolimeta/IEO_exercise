# IEO post-doc exercise

## Overview

This repository contains the complete analysis pipeline and results for both **ChIP-seq** and **RNA-seq** experiments, focusing on differential binding and expression analysis. The analysis was performed on data from an experiment targeting **H3K27ac** for ChIP-seq and gene expression profiling for RNA-seq, with the aim of identifying changes following a gene knockdown.

## Contents

### Nextflow Pipeline

The core pipeline for processing raw data is implemented in **Nextflow**. The following tools are integrated into the workflow:

- **Fastp**: Quality filtering and trimming of raw sequencing reads.
- **Bowtie2**: Alignment of ChIP-seq reads to the reference genome.
- **Kallisto**: Quantification of RNA-seq transcript abundance using pseudo-alignment.
- **MACS2**: Peak calling on ChIP-seq data for identifying enriched regions of H3K27ac.
- **Samtools**: Sorting, indexing, and other manipulations of BAM files.

The Nextflow pipeline processes the data from raw FASTQ files to final peak calls and transcript quantifications. All intermediate files (e.g., BAM, sorted BAM, and peak files) and results can be found in the `results` directory after running the pipeline.

### RMarkdown Analysis

Once the primary processing is complete, further analysis is conducted in R using **RMarkdown** scripts.

- **ChIP-seq differential peak analysis** and **RNA-seq differential expression analysis** are performed in `scripts/analysis/`.
- The **knitted HTML reports** of the main analysis can be found in the `scripts/analysis` directory. These reports provide the full analysis pipeline, including:
  - Quality control.
  - Differential expression analysis using **EdgeR**.
  - Peak annotation using **ChIPseeker** and gene enrichment analysis.

## How to Use

1. **Running the Nextflow Pipeline**:
   - Make sure **Nextflow** and the required tools are installed.
   - Clone this repository and set the paths for your input files (FASTQ and reference files) in the Nextflow configuration file.
   - Run the pipeline with:
     ```
     nextflow run main.nf
     ```

2. **Exploring the RMarkdown Analysis**:
   - Go to the `scripts/analysis` directory.
   - Open the RMarkdown files for each analysis type (e.g., ChIP-seq or RNA-seq).
   - The knitted HTML reports are also available for immediate viewing.

## Dependencies

- **Nextflow**
- **R** with the following packages:
  - **tidyverse**
  - **edgeR**
  - **ChIPseeker**
  - **GenomicRanges**
  - **TxDb.Hsapiens.UCSC.hg38.knownGene**
  - **AnnotationDbi**
  - **biomaRt**
  - **org.Hs.eg.db**
  - **ggrepel**
  - **tidyHeatmap**
  - **piano**

## Notes

- The **ChIP-seq differential peak analysis** was limited by the small sample size, which impacts the statistical power. For better results, future analyses should include more replicates and use dedicated tools like **DiffBind**. A lot of differential peaks were lost due to missing values and similar in the statistical analysis for now.
- The **RNA-seq analysis** was performed with **EdgeR** to identify differentially expressed genes, with normalization using TMM. Gene set analysis was performed with Piano.

## License

This repository is released under [MIT License](LICENSE).
