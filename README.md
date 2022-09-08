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

### Update
If you already have the pipeline installed, you can update it using the following commands:
```
# Navigate to your installation directory
cd de-novo-assembly-pipeline

# Use git pull to get the lastest update
git pull

# Activate the conda environment and use the environment.yml file to download updates
conda activate DeNovoAssembly
conda env update --file environment.yml --prune
```

## Usage

To run the pipeline, use the following command:
```
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR [--host_fasta HOST_FASTA | --host_bt2_index INDEX_DIR]
```

### Options
The pipeline also supports the following optional arguments:

| Option | Type | Description |
|---|---|---|
| --unicycler | *None* | This pipeline also supports assembly using the tool unicycler. Supplying this option will use unicycler as opposed to SPADES (Default). |
| --ref | *Fasta File* | If provided a reference fasta file, this pipeline will align the generated contigs/scaffolds to this file using minimap2 |
| --minCov | *int* | The minimum depth of coverage, below which a a position will be masked. [Default = 20] |
| --minLen | *int* | The minimum length of a read (in base pairs) to keep post trimming. [Default = 75] |
| --threads | *int* | The number of threads that can be use to run pipeline tools in parallel. [Default = 1] |

To view the list of options from the command line, use the following command:
```
nextflow rn main.nf --help
```