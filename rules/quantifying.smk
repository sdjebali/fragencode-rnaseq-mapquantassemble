
##############################################
# Prepare indices for quantifying the gtf    #
##############################################
rule idxquant:
    input:
        gtf = "{annot}/assembling/clean.annot.gtf",
        genome = GENOME
    output:
        "{annot}/assembling/RSEM.grp",
        "{annot}/assembling/RSEM.ti",
        "{annot}/assembling/RSEM.chrlist",
        "{annot}/assembling/RSEM.transcripts.fa",
        "{annot}/assembling/RSEM.seq",
        "{annot}/assembling/RSEM.idx.fa",
        "{annot}/assembling/RSEM.n2g.idx.fa"
    params:
        rsem_prefix = "{annot}/assembling/RSEM"
    threads:
        4
    log:
        "logs/indexing/transcriptome/idx.transcriptome.{annot}.log"
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
        "{annot}/assembling/RSEM.grp",
        "{annot}/assembling/RSEM.ti",
        "{annot}/assembling/RSEM.chrlist",
        "{annot}/assembling/RSEM.transcripts.fa",
        "{annot}/assembling/RSEM.seq",
        "{annot}/assembling/RSEM.idx.fa",
        "{annot}/assembling/RSEM.n2g.idx.fa",
        bam="{annot}/mapping/{tiss}/{anim}/{sample}/Aligned.toTranscriptome.out.bam"
    output:
        tr="{annot}/quantifying/{tiss}/{anim}/{sample}/Quant.isoforms.results",
        gn="{annot}/quantifying/{tiss}/{anim}/{sample}/Quant.genes.results"
    params:
        rsem_prefix = "{annot}/assembling/RSEM",
        out_prefix = "{annot}/quantifying/{tiss}/{anim}/{sample}/Quant"
    threads:
        4
    log:
        "logs/quantification/{tiss}/{anim}/{sample}.{annot}.quant.log"
    conda:
        "../envs/rsem.yaml"
    shell:
        "export PERL5LIB=\"\"; "
        "rsem-calculate-expression --estimate-rspd --calc-ci --no-bam-output "
        "--seed 12345  --ci-memory 3000 --strand-specific --forward-prob 0 "
        "-p {threads} --alignments --paired-end  {input.bam} "
        "{params.rsem_prefix} {params.out_prefix} 2> {log} "
