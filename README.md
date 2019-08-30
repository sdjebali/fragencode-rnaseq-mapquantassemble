# FR-AgENCODE RNA-seq pipeline
This is the pipeline to process RNA-seq data from raw reads to reference and novel gene and transcript expression used in the FR-AgENCODE project http://www.fragencode.org/.
It can be used in the context of other projects but data need to be paired-end and stranded (type mate2_sense).
Also, provided pairs of fastq need to correspond to bioreplicates and not to technical replicates or runs.

## Download the code
<pre>
git clone https://github.com/sdjebali/fragencode-rnaseq-mapquantassemble.git
</pre>

## Installation
Use conda
<pre>
conda env create --file environment.yaml
</pre>

If everything went fine activate the environment
<pre>
conda activate mapquantassemble
conda activate base
</pre>


## Usage
**************************************
First make sure that the files samples.tsv, reads.tsv and config.yaml are correct


## Running snakemake (on genologin here)
<pre>
export DRMAA_LIBRARY_PATH="/tools/libraries/slurm-drmaa/slurm-drmaa-1.0.7/lib/libdrmaa.so"
snakemake --use-conda --debug-dag --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" --configfile config.yaml -n -p
</pre>
