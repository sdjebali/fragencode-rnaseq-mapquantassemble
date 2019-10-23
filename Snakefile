import pandas as pd
#from snakemake.utils import validate

include: "rules/common.smk"

##--------------------------------------
##  Global variables
##--------------------------------------
GENOME = config["genome"]
GTF = config["gtf"]

WORKDIR = config["workdir"]
message("The current working directory is "+WORKDIR)


# paths to various files
snakedir = os.path.dirname(workflow.snakefile)
genomedir = os.path.dirname(GENOME)
gtfdir = os.path.dirname(GTF)


##### load config and read file #####
#configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["sample_file"]).set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

reads = pd.read_table(config["read_file"]).set_index("sample", drop=False)
#validate(reads, schema="schemas/reads.schema.yaml")


wildcard_constraints:
    sample="\w+"


###############
# main rule   #
###############
rule all:
    input:
        expand("{annot}/quantifying/{tiss}/{anim}/{sample}/Quant.isoforms.results", annot=["ref", "new"], sample=samples.index, tiss=samples["tissue"], anim=samples["animal"]),
        expand("{annot}/quantifying/{tiss}/{anim}/{sample}/Quant.genes.results", annot=["ref", "new"], sample=samples.index, tiss=samples["tissue"], anim=samples["animal"])


#### load rules #####
include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/quantifying.smk"
include: "rules/assembling.smk"
include: "rules/lncrnas.smk"
