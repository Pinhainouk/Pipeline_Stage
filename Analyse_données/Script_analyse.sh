#!/bin/bash -i

list="gold_12x_on_data_elodie gold_24x_on_data_elodie gold_40x_on_data_elodie"

for sample in $list
  do

path="/home/elodie/Documents"
aln_mem_sort="${path}/Alignment/${sample}_aln_mem_sort.bam"
rmduplicates="${path}/Alignment/${sample}_aln_mem_sort_rmduplicates.bam"
aln_mem_sort_bamcov="${path}/Analysis/bam_coverage/${sample}_aln_mem_sort_bamcov.bw"
apply_bqsr="${path}/Alignment/${sample}_aln_mem_sort_rmduplicates_apply_bqsr.bam"
apply_bqsr_bamcov="${path}/Analysis/bam_coverage/${sample}_aaply_bqsr_bamcov.bw"
ref_bam="${path}/Elodie/Datas_ismael_ref/gold_40x_on_data_elodie.bam"
ref_bamcov="${path}/Elodie/Datas_ismael_ref/bam_coverage/gold_40x_on_data_elodie_bamcov.bw"
vcf_ref="${path}/Elodie/Datas_ismael_ref/gold_on_data_elodie_ismael_decomposed_normalize_uniq_chr.vcf"
vcf_uniq="${path}/Analysis/${sample}_filtervarianttranches_2D_decomposed_normalized_uniq.vcf"
bed="${path}/Elodie/01_data_elodie_V2.bed"
reference="${path}/Genome/hg19.p13.plusMT.no_alt_analysis_set.fa.gz"

################################################################################
# BAMCOVERAGE SUR LES BAM AVANT ET APRES ENLEVEMENT DES DUPLICATES bin=10
################################################################################
#BamCoverage sur les échantillons

if [ ! -f "${aln_mem_sort_bamcov}" ] ;then
  echo "Les bamcoverage sur les bam n'existent pas";
bamCoverage -b ${aln_mem_sort} --binSize 10 -o ${aln_mem_sort_bamcov}
else
  echo "Les bamcoverage sur les bam existent déjà"
fi

if [ ! -f "${apply_bqsr_bamcov}" ] ;then
  echo "Les bamcoverage sur les bam après rmduplicates n'existent pas";
bamCoverage -b ${apply_bqsr} --binSize 10 -o ${apply_bqsr_bamcov}
else
  echo "Les bamcoverage sur les bam après rmduplicates existent déjà"
fi
done

#BamCoverage sur la référence
if [ ! -f "${ref_bamcov}" ] ;then
  echo "Le bamcoverage sur la référence n'existe pas";
bamCoverage -b ${ref_bam} --binSize 10 -o ${ref_bamcov}
else
  echo "Le bamcoverage sur la référence existe déjà"
fi

################################################################################
# HAP.PY commande lancée sur le cluster de MOABI
################################################################################

#export HGREF=/srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz

# Commandes lancées sur le cluster
# Comparaison 12x-ref
#hap.py /srv/gold_on_data_elodie_ismael_decomposed_normalize_uniq_copie_chr.vcf.gz \
#/srv2/gold_12x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_uniq.vcf.gz \
#-r /srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
#-o /srv/ref-12x \
#-f /srv4/01_data_elodie_V2.bed

# Comparaison 24x-ref
#hap.py /srv/gold_on_data_elodie_ismael_decomposed_normalize_uniq_copie_chr.vcf.gz \
#/srv2/gold_24x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_uniq.vcf.gz \
#-r /srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
#-o /srv/ref-24x \
#-f /srv4/01_data_elodie_V2.bed

# Comparaison 40x-ref
#hap.py /srv/gold_on_data_elodie_ismael_decomposed_normalize_uniq_copie_chr.vcf.gz \
#/srv2/gold_40x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_uniq.vcf.gz \
#-r /srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
#-o /srv/ref-40x \
#-f /srv4/01_data_elodie_V2.bed

# Commande pour lancer en local (essayer d'installer hap.py dans un docker)
#hap.py ${vcf_ref} \
#${vcf_uniq} \
#-r ${reference} \
#-o ${sample}_ref \
#-f ${bed}

################################################################################
# Lancement du Rmarkdown
################################################################################

Rscript -e "rmarkdown::render('/home/elodie/Documents/Elodie/Pipeline_Stage/Analyse_données/Analyse.Rmd',
output_file='/home/elodie/Documents/Elodie/Pipeline_Stage/Analyse_données/Analyse.html')"
