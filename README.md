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
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR
```

### Options
The pipeline also supports the following optional arguments:

| Option | Type | Description |
|---|---|---|
| --unicycler | *None* | This pipeline also supports assembly using the tool unicycler. Supplying this option will use unicycler as opposed to SPADES (Default). |
| --host_fasta | *Fasta File* | A host fasta file that the reads will be aligned to to remove host contamination. If genome is large, creating this index will slow the pipeline. |
| --host_bt2_index | *Directory* | A directory containing an existing bowtie2 index to be used for host read removal. Must be in its own directory. If using a large genome, this option will greatly improve the pipeline runtime. |
| --ref | *Fasta File* | If provided a reference fasta file, this pipeline will align the generated contigs/scaffolds to this file using minimap2 |
| --minLen | *int* | The minimum length of a read (in base pairs) to keep post trimming. [Default = 75] |
| --minTrimQual | *int* | Sets the average basecall quality threshold below which to trim a read. During trimming, trimmomatic performs a sliding window checking the average base quality, and removing the rest of the read if it drops below this treshold. [Default = 20] |
| --threads | *int* | The number of threads that can be use to run pipeline tools in parallel. [Default = 1] |

To view the list of options from the command line, use the following command:
```
nextflow rn main.nf --help
```