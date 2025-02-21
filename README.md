Below is a revised example of a README you could use for your **CommonBean_DE_RNASeq_Pipeline**. It describes the overall workflow for analyzing RNA-seq data using Kallisto and Sleuth, focusing on two common bean genotypes under control and phosphorus stress conditions. Feel free to adjust any sections to align with your specific lab requirements, software versions, or naming conventions.

---

# CommonBean_DE_RNASeq_Pipeline

This repository contains a pipeline for conducting RNA-seq differential expression analysis on *Phaseolus vulgaris* (common bean) using [Kallisto](https://pachterlab.github.io/kallisto/) for quantification and [Sleuth](https://pachterlab.github.io/sleuth/) for differential expression analysis. It is designed to handle multiple genotypes and treatments.

## Table of Contents
- [Background](#background)
- [Dataset Overview](#dataset-overview)
- [Pipeline Overview](#pipeline-overview)
  - [Step 1: Data Download](#step-1-data-download)
  - [Step 2: Quality Control](#step-2-quality-control)
  - [Step 3: Quantification with Kallisto](#step-3-quantification-with-kallisto)
  - [Step 4: Differential Expression with Sleuth](#step-4-differential-expression-with-sleuth)
  - [Step 5: Visualization & Results](#step-5-visualization--results)
- [Dependencies & Installation](#dependencies--installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [References & Tutorials](#references--tutorials)
- [License](#license)

---

## Background

This pipeline was initially conceived to support an RNA-seq study investigating gene expression changes under drought stress in *Phaseolus vulgaris*. To prototype the workflow, we used a publicly available dataset from [Cichy et al., 2019](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0210428) (BioProject [PRJNA498535](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA498535)), which investigated two genotypes of common bean under control and phosphorus (P) stress conditions.

## Dataset Overview

- **Genotypes**: DOR364 and IAC Imperador  
- **Treatments**: Control (no phosphorus stress) vs. phosphorus stress  
- **Replicates**: 3 biological replicates per treatment, per genotype  
- **Total conditions**: 4 (DOR364-Control, DOR364-P-stress, IACImperador-Control, IACImperador-P-stress)

**SRA Run IDs** (for quick reference):
- **DOR364**  
  - Control: SRR8113212, SRR8113218, SRR8113219  
  - P-stress: SRR8113217, SRR8113210, SRR8113213  
- **IAC Imperador**  
  - Control: SRR8113220, SRR8113221, SRR8113211  
  - P-stress: SRR8113214, SRR8113215, SRR8113216

**Reference Genome**: [*Phaseolus vulgaris* genome](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Pvulgaris) from JGI.

---

## Pipeline Overview

### Step 1: Data Download
1. **RNA-seq Data**  
   - Download raw fastq files from the [NCBI BioProject PRJNA498535](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA498535) using `fastq-dump` (SRA Toolkit) or similar tools.  
2. **Reference Genome**  
   - Download from the [JGI portal](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Pvulgaris).  

*Example shell script:*  
```
bash download_data.sh
```
*(Adjust to your own script names and usage.)*

### Step 2: Quality Control
- **QC Tools**: While this dataset was already trimmed, you may perform additional QC with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [MultiQC](https://multiqc.info/) to check read quality.  
- **Trimming**: (Optional) If needed, use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for adapter removal.

### Step 3: Quantification with Kallisto
1. **Indexing**  
   - Create a Kallisto index using the *Phaseolus vulgaris* reference transcriptome FASTA.  
   ```
   kallisto index -i pvulgaris_index.idx <reference_transcriptome>.fasta
   ```
2. **Quantification**  
   - For each sample (paired-end example):
   ```
   kallisto quant -i pvulgaris_index.idx -o <output_directory> \
       -b 100 <R1.fastq.gz> <R2.fastq.gz>
   ```
   - This step generates abundance estimates (`abundance.h5`, `abundance.tsv`).

### Step 4: Differential Expression with Sleuth
1. **Prepare Sample Table**  
   - Construct a metadata table (`sample_id`, `condition`, `genotype`, `path_to_quantification`, etc.).
2. **R Script**  
   - Use the [Sleuth R package](https://pachterlab.github.io/sleuth/) to read the quantification data and perform differential expression analysis.  
   - An example script is provided in `slueth_commonbean.R`.
3. **Modeling**  
   - Set up the appropriate model for your experimental design, e.g.  
   ```r
   so <- sleuth_prep(sample_table, full_model = ~ condition + genotype, ...)
   so <- sleuth_fit(so, ~ condition + genotype, 'full')
   so <- sleuth_fit(so, ~ genotype, 'reduced')
   so <- sleuth_lrt(so, 'reduced', 'full')
   ```
   - Adjust your model formula according to the factors you want to test.

### Step 5: Visualization & Results
- **Volcano Plot / MA Plot**  
  - Use base R, ggplot2, or Sleuth’s plotting functions to visualize significantly differentially expressed genes/transcripts.
- **Outputs**  
  - A CSV of significantly differentially expressed genes (`sleuth_significant.csv`)  
  - Additional diagnostic plots and PCA

---

## Dependencies & Installation
- [Kallisto](https://pachterlab.github.io/kallisto/) (Version >= 0.46.0)
- [R](https://www.r-project.org/) (Version >= 4.0)  
  - [Sleuth](https://pachterlab.github.io/sleuth/)  
  - [tidyverse](https://www.tidyverse.org/)  
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (optional)  
- [SRA Toolkit](https://github.com/ncbi/sra-tools) for data download (optional)

Install R packages via CRAN or Bioconductor as needed:
```r
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
install.packages("tidyverse")
```
*(Check for additional dependencies in your scripts.)*

---

## Usage
1. **Clone this repository**  
   ```
   git clone https://github.com/yourusername/CommonBean_DE_RNASeq_Pipeline.git
   ```
2. **Prepare directories** for raw data, reference genome, and output:
   ```
   mkdir data/ reference/ results/ logs/
   ```
3. **Download data** (Step 1)  
   ```
   bash download_data.sh
   ```
4. **Perform QC** (Step 2; optional if already done)
   ```
   bash run_fastqc.sh
   ```
5. **Run Kallisto** (Step 3)
   ```
   bash run_kallisto.sh
   ```
6. **Execute Sleuth** (Step 4)
   ```r
   Rscript slueth_commonbean.R
   ```
7. **Review and visualize** (Step 5)  
   - Check the `results/` folder for significant transcripts, plots, and logs.

---

## Project Structure
```
CommonBean_DE_RNASeq_Pipeline/
├── data/
│   └── fastq/                  # Raw FASTQ files (downloaded)
├── reference/
│   ├── pvulgaris.fasta         # Reference transcriptome
│   └── pvulgaris_index.idx     # Kallisto index
├── scripts/
│   ├── download_data.sh
│   ├── run_fastqc.sh
│   ├── run_kallisto.sh
│   └── ...
├── R/
│   └── slueth_commonbean.R     # Main Sleuth script
├── results/
│   ├── kallisto_output/
│   ├── sleuth_significant.csv
│   ├── figures/
│   └── logs/
└── README.md
```

---

## References & Tutorials
- [Cichy et al., 2019, *PLOS ONE*](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0210428)  
- [Harvard Informatics RNA-seq Workshop](https://informatics.fas.harvard.edu/workshops/HarvardInformatics_DEworkshop_Fall2017.html)  
- [Sleuth Walkthroughs](https://pachterlab.github.io/sleuth_walkthroughs/)  
- [DESeq2 Workshop Materials](https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)  
- [Formula Notation in R](https://faculty.chicagobooth.edu/richard.hahn/teaching/formulanotation.pdf)

---

## License
This project is publicly available under an open-source license (e.g., MIT). Feel free to modify as needed for your research purposes.

---

### Notes
- Ensure you have the appropriate permissions and credentials (if required) when downloading the reference genome from JGI.  
- Adjust the number of bootstrap samples (`-b` option in Kallisto) or the model formula in Sleuth depending on your experimental design.  
- For larger-scale or production runs, consider containerizing the environment with Docker or Singularity to maintain reproducibility.

---

**Happy analyzing!** If you have any questions or suggestions, please open an issue or contact the project maintainer.
