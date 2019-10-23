# FR-AgENCODE RNA-seq pipeline
This is the pipeline to process RNA-seq data from raw reads to reference and novel gene and transcript expression used in the FR-AgENCODE project http://www.fragencode.org/.
It can be used in the context of other projects but data need to be paired-end and stranded (type mate2_sense).
Also, provided pairs of fastq need to correspond to bioreplicates and not to technical replicates or runs.

Follow the instructions below to use this pipeline on genologin.

## Prerequisites

Conda and Snakemake are required to install and run this pipeline.

1. Load modules
    ```
    module purge
    module load system/Anaconda3-5.2.0
    module load bioinfo/snakemake-4.8.0
    ```

2. Enable the conda activate command
    ```
    source /usr/local/bioinfo/src/Anaconda/Anaconda3-5.2.0/etc/profile.d/conda.sh
    ```

## Installation

1. Download the code
    ```
    git clone https://github.com/sdjebali/fragencode-rnaseq-mapquantassemble.git
    ```

2. (Recommended) Specify the conda environments and packages paths if you want them elsewhere than ~/.conda
    ```
    conda config --add envs_dirs /path/to/conda/envs
    conda config --add pkgs_dirs /path/to/conda/pkgs
    ```

3. Create the environment in the conda environments path
    ```
    cd fragencode-rnaseq-mapquantassemble
    conda env create --file environment.yaml
    ```

## Usage

1. Configure the samples.tsv, reads.tsv and config.yaml files

2. Activate the conda environment
    ```
    conda activate mapquantassemble
    ```

3. Export the DRMAA library path
   ```
   export DRMAA_LIBRARY_PATH="/tools/libraries/slurm-drmaa/slurm-drmaa-1.0.7/lib/libdrmaa.so"
   ```

4. Run the pipeline
   ```
   cd fragencode-rnaseq-mapquantassemble
   snakemake --use-conda --conda-prefix /path/to/conda/envs --debug-dag --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" --configfile config.yaml -p
   ```
