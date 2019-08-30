
#####################################
# Assemble transcripts per sample   #
#####################################
rule assemble:
    input:
        "ref/mapping/{tiss}/{anim}/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        tr="new/assembling/{tiss}/{anim}/{sample}/transcripts.gtf",
        gn="new/assembling/{tiss}/{anim}/{sample}/genes.fpkm_tracking",
        iso="new/assembling/{tiss}/{anim}/{sample}/isoforms.fpkm_tracking",
        skp="new/assembling/{tiss}/{anim}/{sample}/skipped.gtf"
    params:
        gtf=GTF,
        outdir="new/assembling/{tiss}/{anim}/{sample}"
    threads:
        4
    log:
        "logs/assembling/{tiss}/{anim}/{sample}.assemble.log"
    shell:
        "cufflinks --max-intron-length 100000 --num-threads {threads} "
        "--overlap-radius 5 --library-type fr-firststrand "
        "--GTF-guide {params.gtf} -o {params.outdir} {input} 2> {log}"


###################################################################
# Merge the transcripts of all samples into a single new model    #
###################################################################
rule merge:
    input:
        expand("new/assembling/{tiss}/{anim}/{sample}/transcripts.gtf", sample=samples.index, tiss=samples["tissue"], anim=samples["animal"])
    output:
        "new/assembling/merged.gtf"
    params:
        gtf=GTF,
        genome=GENOME,
        outdir="new/assembling"
    threads:
        4
    log:
        "logs/merging/merge.log"
    conda:
        "../envs/cuffmerge.yaml"
    shell:
        "ls -1 {input} > {params.outdir}/gtf_list.txt; "
        "cuffmerge --num-threads {threads} --ref-gtf {params.gtf} "
        "--ref-sequence {params.genome} -o {params.outdir} "
        "{params.outdir}/gtf_list.txt 2> {log}"
