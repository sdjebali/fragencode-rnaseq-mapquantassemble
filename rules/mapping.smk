
##############################################################################
# Index the genome and gene model for mapping to genome and reference model  #
##############################################################################
rule idxmap:
    input:
        genome = GENOME,
        gtf = "{annot}/assembling/clean.annot.gtf"
    output:
        "{annot}/assembling/genomeParameters.txt",
        "{annot}/assembling/chrStart.txt",
        "{annot}/assembling/chrNameLength.txt",
        "{annot}/assembling/chrName.txt",
        "{annot}/assembling/chrLength.txt",
        "{annot}/assembling/geneInfo.tab",
        "{annot}/assembling/exonGeTrInfo.tab",
        "{annot}/assembling/transcriptInfo.tab",
        "{annot}/assembling/exonInfo.tab",
        "{annot}/assembling/sjdbList.fromGTF.out.tab",
        "{annot}/assembling/sjdbList.out.tab",
        "{annot}/assembling/sjdbInfo.txt",
        "{annot}/assembling/SAindex",
        "{annot}/assembling/SA",
        "{annot}/assembling/Genome"
    params:
        readlength=config["read_length"],
        stardir="{annot}/assembling"
    threads:
        4
    log:
        "logs/indexing/genome/idx.genome.{annot}.log"
    shell:
        "overhang=`expr {params.readlength} - 1`; "
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {params.stardir} --genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} --sjdbOverhang $overhang "
        "--outFileNamePrefix . --limitGenomeGenerateRAM 120000000000 2> {log}"

##########################################
# Map trimmed reads to genome and to gtf #
##########################################
rule mapref:
    input:
        "{annot}/assembling/genomeParameters.txt",
        "{annot}/assembling/chrStart.txt",
        "{annot}/assembling/chrNameLength.txt",
        "{annot}/assembling/chrName.txt",
        "{annot}/assembling/chrLength.txt",
        "{annot}/assembling/geneInfo.tab",
        "{annot}/assembling/exonGeTrInfo.tab",
        "{annot}/assembling/transcriptInfo.tab",
        "{annot}/assembling/exonInfo.tab",
        "{annot}/assembling/sjdbList.fromGTF.out.tab",
        "{annot}/assembling/sjdbList.out.tab",
        "{annot}/assembling/sjdbInfo.txt",
        "{annot}/assembling/SAindex",
        "{annot}/assembling/SA",
        "{annot}/assembling/Genome",
        read1="trimming/{tiss}/{anim}/{sample}/{sample}_R1_val_1.fq.gz",
        read2="trimming/{tiss}/{anim}/{sample}/{sample}_R2_val_2.fq.gz"
    output:
        "{annot}/mapping/{tiss}/{anim}/{sample}/Log.final.out",
        "{annot}/mapping/{tiss}/{anim}/{sample}/Aligned.sortedByCoord.out.bam",
        "{annot}/mapping/{tiss}/{anim}/{sample}/Aligned.toTranscriptome.out.bam"
    params:
        stardir="{annot}/assembling",
        prefix="{annot}/mapping/{tiss}/{anim}/{sample}/"
    threads:
        4
    log:
        "logs/mapping/{tiss}/{anim}/{sample}.map.{annot}.log"
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
