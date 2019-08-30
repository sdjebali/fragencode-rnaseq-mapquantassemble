

#################################
# Trim reads to remove adaptors #
#################################
rule trimming:
    input:
        read1=get_fastq1,
        read2=get_fastq2
    output:
        read1="trimming/{tiss}/{anim}/{sample}/{sample}_R1_val_1.fq.gz",
        read2="trimming/{tiss}/{anim}/{sample}/{sample}_R2_val_2.fq.gz"
    params:
        readlength=config["read_length"],
        seqfstread=config["adapt_3p_first_read"],
        seqsndread=config["adapt_3p_second_read"]
    log:
        "logs/trimming/{tiss}/{anim}/{sample}.log"
    shell:
        "minLength=`expr {params.readlength} / 3`; "
        "cutadapt -a {params.seqfstread} -A {params.seqsndread} "
        "--minimum-length $minLength --output {output.read1} "
        "--paired-output {output.read2} {input.read1} {input.read2} 2> {log}"


#####################################################
# Remove unstranded transcripts from the gtf file   #
#####################################################
rule cleangtf:
    input:
        get_gtf
    output:
        "{annot}/assembling/clean.annot.gtf"
    params:
        path=snakedir
    log:
        "logs/cleaning/clean.{annot}.log"
    shell:
        "awk -f {params.path}/scripts/keep.stranded.awk {input} > {output} 2> {log}"
