
#####################################
# Assemble transcripts per sample   #
#####################################
rule assemble:
    input:
        os.path.join(WORKDIR, "ref/mapping/{sample}/Aligned.sortedByCoord.out.bam")
    output:
        tr=os.path.join(WORKDIR, "new/assembling/{sample}/transcripts.gtf"),
        gn=os.path.join(WORKDIR, "new/assembling/{sample}/genes.fpkm_tracking"),
        iso=os.path.join(WORKDIR, "new/assembling/{sample}/isoforms.fpkm_tracking"),
        skp=os.path.join(WORKDIR, "new/assembling/{sample}/skipped.gtf")
    params:
        gtf=GTF,
        outdir=os.path.join(WORKDIR, "new/assembling/{sample}")
    threads:
        4
    log:
        os.path.join(WORKDIR, "logs/assembling/{sample}.assemble.log")
    shell:
        "cufflinks --max-intron-length 100000 --num-threads {threads} "
        "--overlap-radius 5 --library-type fr-firststrand "
        "--GTF-guide {params.gtf} -o {params.outdir} {input} 2> {log}"


###################################################################
# Merge the transcripts of all samples into a single new model    #
###################################################################
rule merge:
    input:
        expand(os.path.join(WORKDIR, "new/assembling/{sample}/transcripts.gtf"), sample=samples.index, tiss=samples["tissue"], anim=samples["animal"])
    output:
        os.path.join(WORKDIR, "new/assembling/merged.gtf")
    params:
        gtf=GTF,
        genome=GENOME,
        outdir=os.path.join(WORKDIR, "new/assembling")
    threads:
        4
    log:
        os.path.join(WORKDIR, "logs/merging/merge.log")
    shell:
        "ls -1 {input} > {params.outdir}/gtf_list.txt; "
        "cuffmerge --num-threads {threads} --ref-gtf {params.gtf} "
        "--ref-sequence {params.genome} -o {params.outdir} "
        "{params.outdir}/gtf_list.txt 2> {log}"
