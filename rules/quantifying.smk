
##############################################
# Prepare indices for quantifying the gtf    #
##############################################
rule idxquant:
    input:
        gtf = os.path.join(WORKDIR, "{annot}/assembling/clean.annot.gtf"),
        genome = GENOME
    output:
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.grp"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.ti"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.chrlist"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.transcripts.fa"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.seq"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.idx.fa"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.n2g.idx.fa")
    params:
        rsem_prefix = os.path.join(WORKDIR, "{annot}/assembling/RSEM")
    threads:
        4
    log:
        os.path.join(WORKDIR, "logs/indexing/transcriptome/idx.transcriptome.{annot}.log")
    conda:
        "../envs/rsem.yaml"
    shell:
        "export PERL5LIB=\"\"; "
        "rsem-prepare-reference -p {threads} --gtf {input.gtf} {input.genome} "
        "{params.rsem_prefix} 2> {log}"


################################
# Quantification of the gtf    #
################################
rule quant:
    input:
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.grp"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.ti"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.chrlist"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.transcripts.fa"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.seq"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.idx.fa"),
        os.path.join(WORKDIR, "{annot}/assembling/RSEM.n2g.idx.fa"),
        bam=os.path.join(WORKDIR, "{annot}/mapping/{tiss}/{anim}/{sample}/Aligned.toTranscriptome.out.bam")
    output:
        tr=os.path.join(WORKDIR, "{annot}/quantifying/{tiss}/{anim}/{sample}/Quant.isoforms.results"),
        gn=os.path.join(WORKDIR, "{annot}/quantifying/{tiss}/{anim}/{sample}/Quant.genes.results")
    params:
        rsem_prefix = os.path.join(WORKDIR, "{annot}/assembling/RSEM"),
        out_prefix = os.path.join(WORKDIR, "{annot}/quantifying/{tiss}/{anim}/{sample}/Quant")
    threads:
        4
    log:
        os.path.join(WORKDIR, "logs/quantification/{tiss}/{anim}/{sample}.{annot}.quant.log")
    conda:
        "../envs/rsem.yaml"
    shell:
        "export PERL5LIB=\"\"; "
        "rsem-calculate-expression --estimate-rspd --calc-ci --no-bam-output "
        "--seed 12345  --ci-memory 3000 --strand-specific --forward-prob 0 "
        "-p {threads} --alignments --paired-end  {input.bam} "
        "{params.rsem_prefix} {params.out_prefix} 2> {log} "
