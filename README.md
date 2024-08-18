
# Proteomics Data Processing Pipeline

<img src="/SagharT/MasterThesis/raw/main/Images/Workflow.png" alt="Workflow" width="70%">

This repository hosts a collection of scripts and a Snakemake workflow designed to streamline the analysis of proteomics data, accommodating both Data-Independent Acquisition (DIA) and Data-Dependent Acquisition (DDA) methods. The pipeline handles raw data conversion, feature extraction, peptide identification, quality control and generation of detailed analytical reports.

## Components of the Pipeline

The following scripts and configuration files are included in this repository:

- **`Compare.py`**: Compares features across samples to identify significant trends and outliers.
- **`DDA-Report.py`**: Generates reports for DDA analyzed samples, including data visualizations.
- **`DIA-Report.py`**: Similar to `DDA-Report.py`, but for DIA data.
- **`DIA-or-DDA.py`**: Detects whether the samples are DDA or DIA and based on this, the samples will proceed under DIA or DDA protocols .
- **`MZid-Parser.py`**: Parses .mzid files to extract peptide identification data.
- **`cleanup_empty_files.py`**: Removes empty files from output directories by the end of process.
- **`config.yaml`**: A configuration file that users must edit to specify paths and parameters relevant to their data and system environment.
- **`mzML-Parser.py`**: Processes .mzML files to extract essential data for downstream analysis.
- **`snakefile`**: Coordinates the workflow defined for Snakemake, ensuring each step of the process is executed in the correct sequence.

## Workflow Description

This Snakemake pipeline processes proteomics data through several stages:

- **Raw Data Conversion**: Converts raw mass spectrometry data files into .mzML format using the MSConvert tool.
- **Feature Extraction**: Utilizes algorithms like Dinosaur to detect and quantify features from .mzML files.
- **Quality Contro**l: Generates summary reports on the quality of the processed data.
- **Peptide Identification**: Employs MS-GF+ for DDA data and DIANN for DIA data to identify peptides.
- **Report Generation**: Produces detailed reports, including visualizations of the analyzed data.

## Getting Started
### Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/yourusername/proteomics-pipeline.git
cd proteomics-pipeline
```

### Prerequisites


Before you begin setting up and running this workflow, ensure the following software and tools are installed on your system:

##### Software and Tools

- **Python 3.x**: Required for running the scripts. Install the following Python packages:
  - numpy
  - pandas
  - matplotlib
  - reportlab
  - pyteomics
  - lxml

- **Java**: Needed to run Java-based applications such as MSGFPlus and Dinosaur.

- **Snakemake**: Essential for managing and executing the workflow. Install via Python's package manager or conda.

- **Docker**: Used for running the `msconvert` tool within a Docker container. Docker ensures consistency across different computing environments.

##### Proteomics Software

- **MSGFPlus**: A powerful tool for peptide identification. Download the latest version [here](https://github.com/MSGFPlus/msgfplus) (version 2021.03.22).

- **Dinosaur**: An open-source tool for feature detection in mass spectrometry data. Available on [GitHub](https://github.com/fickludd/dinosaur) (version 1.1.4).
- **QC_LCMS.jar**:  Contains FeaturesQCSummary.jar and MzMLSummary.jar. This scripts were developed previously in house and not published yet.

## Configuration

Before running the workflow, ensure that all paths and parameters in the `config.yaml` file are correctly set up according to your system's directory structure and processing needs.

### Running the Workflow
To execute the workflow, navigate to the directory containing the snakefile and run the following command:

```bash
snakemake
```
This will process all .raw files specified in the directory, performing all steps from data conversion to report generation automatically.

