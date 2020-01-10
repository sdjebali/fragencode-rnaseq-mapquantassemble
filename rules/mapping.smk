
##############################################################################
# Index the genome and gene model for mapping to genome and reference model  #
##############################################################################
rule idxmap:
    input:
        genome = GENOME,
        gtf = os.path.join(WORKDIR, "{annot}/assembling/clean.annot.gtf")
    output:
        os.path.join(WORKDIR, "{annot}/assembling/genomeParameters.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrStart.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrNameLength.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrName.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrLength.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/geneInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/exonGeTrInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/transcriptInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/exonInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/sjdbList.fromGTF.out.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/sjdbList.out.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/sjdbInfo.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/SAindex"),
        os.path.join(WORKDIR, "{annot}/assembling/SA"),
        os.path.join(WORKDIR, "{annot}/assembling/Genome")
    params:
        readlength=config["read_length"],
        stardir=os.path.join(WORKDIR, "{annot}/assembling")
    threads:
        4
    log:
        os.path.join(WORKDIR, "logs/indexing/genome/idx.genome.{annot}.log")
    shell:
        "overhang=`expr {params.readlength} - 1`; "
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {params.stardir} --genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} --sjdbOverhang $overhang "
        "--outFileNamePrefix . --limitGenomeGenerateRAM 25000000000 2> {log}"

##########################################
# Map trimmed reads to genome and to gtf #
##########################################
rule mapref:
    input:
        os.path.join(WORKDIR, "{annot}/assembling/genomeParameters.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrStart.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrNameLength.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrName.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/chrLength.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/geneInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/exonGeTrInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/transcriptInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/exonInfo.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/sjdbList.fromGTF.out.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/sjdbList.out.tab"),
        os.path.join(WORKDIR, "{annot}/assembling/sjdbInfo.txt"),
        os.path.join(WORKDIR, "{annot}/assembling/SAindex"),
        os.path.join(WORKDIR, "{annot}/assembling/SA"),
        os.path.join(WORKDIR, "{annot}/assembling/Genome"),
        read1=os.path.join(WORKDIR, "trimming/{sample}/{sample}_R1_val_1.fq.gz"),
        read2=os.path.join(WORKDIR, "trimming/{sample}/{sample}_R2_val_2.fq.gz")
    output:
        os.path.join(WORKDIR, "{annot}/mapping/{sample}/Log.final.out"),
        os.path.join(WORKDIR, "{annot}/mapping/{sample}/Aligned.sortedByCoord.out.bam"),
        os.path.join(WORKDIR, "{annot}/mapping/{sample}/Aligned.toTranscriptome.out.bam")
    params:
        stardir=os.path.join(WORKDIR, "{annot}/assembling"),
        prefix=os.path.join(WORKDIR, "{annot}/mapping/{sample}/")
    threads:
        4
    log:
        os.path.join(WORKDIR, "logs/mapping/{sample}.map.{annot}.log")
    shell:
        # difference with GB version is from --outSAMattributes to --alignSJDBoverhangMin but could be from dcc site
        # that I took from raffine but these could also be defaults parameters so it should not change much the output
        "STAR --genomeDir {params.stardir} --readFilesIn {input.read1} {input.read2} "
        "--readFilesCommand zcat --outFilterType BySJout --outSAMunmapped Within "
        "--outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical " # for cufflinks
        "--alignSoftClipAtReferenceEnds Yes " # for cufflinks and no need for post-processing with custom script then
        "--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 "
        "--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 "
        "--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 "
        "--runThreadN {threads} "
        "--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM "
        "--outWigType bedGraph --outWigStrand Stranded "
        "--outFileNamePrefix {params.prefix} 2> {log}"
