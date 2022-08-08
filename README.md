# Nextflow De Novo Assembly Pipeline

This pipeline automates the process of performing de novo assembly from Illumina NGS Data.

## Installation

To install this pipeline, enter the following command:
```
# Clone the Repository
git clone https://github.com/rchapman2000/de-novo-assembly-pipeline.git

# Create a conda environment using the provided environment.yml file
conda env create -f environment.yml

# Activate the conda environment
conda activate DeNovoAssembly
```

## Usage

To run the pipeline, use the following command:
```
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --ref REFERENCE_FASTA --adapters ADAPTER_FASTA [--host_fasta HOST_FASTA | --host_bt2_index INDEX_DIR]
```

### Options
Documentation to be added. To view the list of options, use the following command:
```
nextflow rn main.nf --help
```