
FEELnc_filter.pl -i merged.gtf -b transcript_biotype=protein_coding,pseudogene --size=200 --monoex=1 --biex=25 --proc=12 --outlog=out.capra_filter.log -a capra_hircus_training_prot_corrected_wo_tRNA_exonOnly_woTTN.gtf > candidate_capra_lncRNA.gtf

FEELnc_codpot.pl -i candidate_capra_lncRNA.gtf -a  capra_hircus_training_prot_corrected_wo_tRNA_exonOnly_woTTN.gtf -g capra_hircus.fa --mode=shuffle -k "1,2,3,6,9,12" --outname="lst_lnc.capra.shuffle_wTrainingProt_98" --outdir="." --spethres=0.98,0.98 --seed=201

cat *.noORF.gtf *.lncRNA.gtf > lst_lnc.capra.shuffle_wTrainingProt_98.lncRNA_noORF.gtf

FEELnc_classifier.pl -i lst_lnc.capra.shuffle_wTrainingProt_98.lncRNA_noORF.gtf -a  capra_hircus_training_prot_corrected_wo_tRNA_exonOnly_woTTN.gtf --window=10000 --maxwindow=1000000 > lncRNA_noORF_classes.txt



