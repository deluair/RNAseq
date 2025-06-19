# RNA-seq Analysis Platform

This repository contains a comprehensive platform for RNA-seq data analysis, from raw data retrieval to differential expression and splicing analysis. It features a user-friendly web interface and a powerful command-line pipeline, making it suitable for both biologists and bioinformaticians. The platform is designed to be robust, with features like automatic retry mechanisms for API calls and fallback mock data generation.

![Web App Screenshot](https://i.imgur.com/your-screenshot.png) <!-- Replace with an actual screenshot -->

## Table of Contents

- [Features](#features)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
  - [Web Application](#web-application)
  - [Demonstration Script](#demonstration-script)
  - [Full Analysis Pipeline](#full-analysis-pipeline)
  - [Manuscript Generation](#manuscript-generation)
- [Modules](#modules)
  - [Data Retrieval](#data-retrieval)
  - [Quality Control](#quality-control)
  - [Alignment](#alignment)
  - [Differential Expression Analysis](#differential-expression-analysis)
  - [Splicing Analysis](#splicing-analysis)
  - [Visualization](#visualization)
- [Contributing](#contributing)
- [License](#license)

## Features

-   **Automated Data Retrieval**: Download data directly from NCBI's Sequence Read Archive (SRA).
-   **Robust Error Handling**: Includes exponential backoff and retry mechanisms for NCBI API calls to handle rate limiting.
-   **Mock Data Generation**: Can generate realistic mock FASTQ files if downloads fail, ensuring the pipeline can run end-to-end for demonstration and testing.
-   **Quality Control**: Integrated FastQC for assessing raw read quality.
-   **Alignment**: Uses the STAR aligner for mapping reads to a reference genome.
-   **Differential Expression Analysis**:
    -   Primary analysis using R/DESeq2.
    -   A pure Python-based implementation as a fallback for systems without R.
    -   Generates volcano plots and MA plots.
-   **Splicing Analysis**: A mock module to identify and visualize alternative splicing events (SE, MXE, A3SS, A5SS, RI).
-   **Web Interface**: A Flask-based web application to run and monitor analysis steps.
-   **Comprehensive Reporting**: Generates detailed HTML reports for each analysis stage.
-   **Scientific Manuscript Generation**: Includes a LaTeX template and `Makefile` to automatically generate a publication-quality manuscript from the analysis results.

## Project Structure

```
rnaseq/
├── app.py                      # Flask web application
├── pipeline.py                 # Main command-line pipeline script
├── demo.py                     # Script to run a quick demonstration
├── run_full_analysis.py        # Script to run the full analysis pipeline with mock data
├── requirements.txt            # Python package dependencies
├── setup.py                    # Project setup script
├── config/                     # Configuration files for analysis runs
├── data/                       # Raw and reference data
├── results/                    # Output directory for analysis results
├── logs/                       # Log files
├── manuscript/                 # LaTeX source and figures for the scientific manuscript
├── src/                        # Main source code
│   ├── data_retrieval/         # Data downloaders (e.g., from NCBI)
│   ├── preprocessing/          # QC and read trimming tools
│   ├── alignment/              # STAR aligner wrapper
│   ├── analysis/               # DE and splicing analysis modules
│   ├── quantification/         # Gene quantification tools
│   └── visualization/          # Plot generation modules
├── static/                     # CSS and JavaScript for the web app
├── templates/                  # HTML templates for the web app
└── README.md                   # This file
```

## Getting Started

### Prerequisites

Before you begin, ensure you have the following dependencies installed:

1.  **Python 3.8+**
2.  **Miniconda or Anaconda**: Recommended for managing environments.
3.  **Bioinformatics Tools**:
    -   **SRA Toolkit**: For downloading data from NCBI.
    -   **FastQC**: For quality control of raw sequencing data.
    -   **STAR**: For aligning reads to the reference genome.
    -   **R and DESeq2** (Optional but recommended): For differential expression analysis.
    -   **LaTeX Distribution** (Optional): For compiling the manuscript (e.g., MacTeX on macOS, TeX Live on Linux).

**On macOS, you can install most of these with Homebrew:**
```bash
brew install sratoolkit
brew install fastqc
brew install star
brew install R
```

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/deluair/RNAseq.git
    cd RNAseq
    ```

2.  **Create and activate a conda environment:**
    ```bash
    conda create -n rnaseq_env python=3.9
    conda activate rnaseq_env
    ```

3.  **Install Python dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

4.  **Install DESeq2 in R** (if you installed R):
    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("DESeq2")
    ```

5.  **Download Reference Genome and Annotation**:
    Place your reference genome (`.fa`) and gene annotation (`.gtf`) files in the `data/reference/` directory. Example files for soybean are provided.

## Usage

### Web Application

To start the web application, run:
```bash
python app.py --port 5001
```
The application will be accessible at `http://127.0.0.1:5001`. The web interface allows you to:
-   Search for and download datasets from NCBI.
-   Run quality control.
-   Perform differential expression and splicing analysis.
-   View results and reports.

### Demonstration Script

The `demo.py` script runs a small, pre-configured analysis on sample data.
```bash
python demo.py
```
This will generate sample reports in the `demo_results/` directory.

### Full Analysis Pipeline

The `run_full_analysis.py` script executes the complete end-to-end pipeline using mock data. This is useful for testing the integration of all modules.
```bash
python run_full_analysis.py
```
This script simulates a full analysis run, including DE and splicing analysis, and generates a comprehensive HTML report in `results/full_analysis_example/`.

### Manuscript Generation

A LaTeX manuscript template is provided in the `manuscript/` directory. After running an analysis, you can compile the manuscript to a PDF.
```bash
cd manuscript/
make
```
This will produce a `soybean_drought_rnaseq.pdf` file.

## Modules

### Data Retrieval (`src/data_retrieval`)
Handles the download of public datasets from NCBI SRA. It is designed to be resilient to network issues and API limits, using a retry strategy with exponential backoff. If a download ultimately fails, it can create a mock FASTQ file to allow downstream steps to proceed.

### Quality Control (`src/preprocessing`)
Uses FastQC to generate quality reports for each FASTQ file, which are essential for identifying issues with raw sequencing data.

### Alignment (`src/alignment`)
A wrapper for the STAR aligner. It takes FASTQ files and a reference genome index to produce BAM alignment files.

### Differential Expression Analysis (`src/analysis`)
-   `deseq2_analyzer.py`: Interfaces with R/DESeq2 for robust DE analysis.
-   `python_de_analyzer.py`: A pure Python implementation providing an alternative when R is not available. It performs normalization and statistical tests to find differentially expressed genes.

### Splicing Analysis (`src/analysis`)
The `splicing_analyzer.py` module provides a framework for identifying alternative splicing events. The current implementation is a mock that generates a realistic report, but it can be extended with tools like rMATS.

### Visualization (`src/visualization`)
Generates various plots for interpreting the analysis results, including volcano plots, MA plots, heatmaps, and splicing event diagrams.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue to discuss proposed changes.

1.  Fork the repository.
2.  Create a new branch (`git checkout -b feature/your-feature-name`).
3.  Make your changes.
4.  Commit your changes (`git commit -m 'Add some feature'`).
5.  Push to the branch (`git push origin feature/your-feature-name`).
6.  Open a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details. 