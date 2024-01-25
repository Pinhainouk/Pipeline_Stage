#!/bin/bash -i

while getopts "o:s:" opt
do
  case "$opt" in
    o) path=${OPTARG};;
    s) sample=${OPTARG};;
  esac
done

if [ -z "$path" ] || [ -z "$sample" ]; then
  echo "Usage: $0 -o path -s sample"
  echo "La commande $opt nécessite une option"
  exit 1
fi
echo "options fournies:"
echo "path : $path"
echo "sample : $sample"

#shopt -s expand_aliases - non utilisé car utlilsation dans le head du -i pour les alias
#source ~/.bashrc - non utilisé

#outils ajoutés dans bashrc en alias

# Création d'une liste avec les noms des samples
# Boucle sur les samples avec les outils : fastqc - trimming - fastqc_trimming - alignement bwa - samtools stats
# Concaténation des chemins des samples avec création de variables
# echo = print

# Création des dossiers pour ranger les fichiers générés par le pipeline

#list_directory="QC_Raw Trimming QC_Trimming Alignment QC_Alignment"
#for directory in $list_directory
#do
#if [ ! -d "$directory" ]; then
#  echo "Les dossiers n'existent pas"
  #-p indique au mkdir pour créer d'abord le répertoire principal s'il n'existe pas déjà
  mkdir -p "${path}/{QC_Raw,Trimming,QC_Trimming,Alignment,QC_Alignment,Analysis}"
#  echo "$directory"
#else
#  echo "Les dossiers ont déjà été créés"
#fi
#done
#activation de l'environnement GATK
#conda activate gatk

end=".fastq.gz"
R1_raw="$path/Raw/${sample}_1${end}"
R2_raw="${path}/Raw/${sample}_2${end}"
log_qc_raw="${path}/Raw/${sample}_log_qc_raw"
R1_trim="${path}/Trimming/${sample}_1_trimming${end}"
R2_trim="${path}/Trimming/${sample}_2_trimming${end}"
stats_trim="${path}/Trimming/${sample}_trimming_stats"
log_trim="${path}/Trimming/${sample}_log_trim"
log_qc_trim="${path}/Trimming/${sample}_log_qc_trim"
aln_mem="${path}/Alignment/${sample}_aln_mem.sam"
log_alignment="${path}/QC_Alignment/${sample}_aln_log"
aln_mem_sort="${path}/Alignment/${sample}_aln_mem_sort.bam"
flagstat="${path}/QC_Alignment/${sample}_aln_mem_sort_flagstat"
idxstats="${path}/QC_Alignment/${sample}_aln_mem_sort_idxstats"
stats="${path}/QC_Alignment/${sample}_aln_mem_sort_stats"
rmduplicates="${path}/Alignment/${sample}_aln_mem_sort_rmduplicates.bam"
rmduplicates_metrics="${path}/QC_Alignment/${sample}_aln_mem_sort_rmduplicates_metrics"
log_rmduplicates="${path}/QC_Alignment/${sample}_aln_mem_sort_rmduplicates_log"
flagstat_rmduplicates="${path}/QC_Alignment/${sample}_aln_mem_sort_rmduplicates_flagstat"
idxstats_rmduplicates="${path}/QC_Alignment/${sample}_aln_mem_sort_rmduplicates_idxstats"
stats_rmduplicates="${path}/QC_Alignment/${sample}_aln_mem_sort_rmduplicates_stats"
read_groups="${path}/Alignment/${sample}_aln_mem_sort_rmduplicates_read_groups.bam"
vcf_indel="${path}/Genome/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
vcf_snp="${path}/Genome/dbsnp_138.hg19.vcf"
reference="${path}/Genome/hg19.p13.plusMT.no_alt_analysis_set.fa.gz"
before_base_recalibrator="${path}/Alignment/${sample}_before_bqsr.report"
apply_bqsr="${path}/Alignment/${sample}_aln_mem_sort_rmduplicates_apply_bqsr.bam"
after_base_recalibrator="${path}/Alignment/${sample}_after_bqsr.report"
before_after_bqsr_plot="${path}/Alignment/${sample}_analyze_covariates_plot.pdf"
vcf="${path}/Analysis/${sample}.vcf.gz"
bamout="${path}/Alignment/${sample}_bamout.bam"
bed="${path}/Elodie/01_data_elodie_V2.bed"
vcf_cnnscorevariants="${path}/Analysis/${sample}_cnnscorevariants_2D.vcf.gz"
vcf_filtervarianttranches="${path}/Analysis/${sample}_filtervarianttranches_2D.vcf"
vcf_decomposed="${path}/Analysis/${sample}_filtervarianttranches_2D_decomposed.vcf"
vcf_normalized="${path}/Analysis/${sample}_filtervarianttranches_2D_decomposed_normalized.vcf"
vcf_funcotated="${path}/Analysis/${sample}_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf"

echo $R1_raw
echo $R2_raw
echo $R1_trim
echo $R2_trim

################################################################################
# FASTQC SUR LES RAW
################################################################################
# condition si fichier de sortie non trouvé alors on exécute la commande

if [ ! -f "${path}/QC_Raw/${sample}_1_fastqc.html" ] && [ ! -f "${path}/QC_Raw/${sample}_2_fastqc.html" ];then
  echo "Les fichiers fastqc n'existent pas";
  fastqc $R1_raw $R2_raw -o /home/elodie/Documents/QC_Raw/ 2>&1 | tee -a $log_qc_raw
else
  echo "Les fichiers fastqc existent déjà"
fi
################################################################################
# TRIMMING avec bbduk
################################################################################

if [ ! -f "$R1_trim" ] && [ ! -f "$R2_trim" ] && [ ! -f "$log_trim" ]; then
  echo "Les fichiers trimming et log trimming n'existent pas"
  bbduk in1=$R1_raw in2=$R2_raw k=27 mink=11 ktrim=r ref=/home/elodie/Documents/Elodie/adapters.fa qtrim=rl trimq=20 out1=$R1_trim out2=$R2_trim \
  stats=$stats_trim tpe 2>&1 | tee -a $log_trim
else
  echo "Les fichiers trimming et log trimming existent déjà"
fi

################################################################################
  # FASTQC SUR LES RAW_TRIM
################################################################################

if [ ! -f "${path}/QC_Trimming/${sample}_1_trimming_fastqc.html" ] && [ ! -f "${path}/QC_Trimming/${sample}_2_trimming_fastqc.html" ]; then
  echo "Les fichiers fastqc_trim n'existent pas"
  fastqc $R1_trim $R2_trim -o /home/elodie/Documents/QC_Trimming/ 2>&1 | tee -a $log_qc_trim
else
  echo "Les fichiers fastqc_trim existent déjà"
fi

################################################################################
# INDEX DU GENOME lancé à part de ce pipeline
################################################################################

################################################################################
# ALIGNEMENT AVEC COMMANDE BWA MEM : fichiers sam générés
################################################################################

if [ ! -f "$aln_mem" ];then
  echo "Les fichiers sam n'existent pas"
  bwa mem /home/elodie/Documents/Genome/hg19_genome_index $R1_trim $R2_trim -o $aln_mem 2>&1 | tee -a $log_alignment
else
  echo "Les fichiers sam existent déjà"
fi

# génération des bam à partir des sam, tri des bams, génération des index des bams triés
if [ ! -f "$aln_mem_sort" ];then
  echo "Les fichiers bam triés et indexés n'existent pas"
  samtools view -b $aln_mem  | samtools sort -o $aln_mem_sort
  samtools index $aln_mem_sort
else
  echo "Les fichiers bam triés et indexés existent déjà"
fi

# stats à partir des bam triés
if [ ! -f "$flagstat" ];then
  echo "Les fichiers flagstat n'existent pas"
samtools flagstat $aln_mem_sort > $flagstat
else
  echo "Les fichiers flagstat existent déjà"
fi

if [ ! -f "$idxstats" ];then
  echo "Les fichiers idxstats n'existent pas"
samtools idxstats $aln_mem_sort > $idxstats
else
  echo "Les fichiers idxstats existent déjà"
fi

if [ ! -f "$stats" ];then
  echo "Les fichiers stats n'existent pas"
samtools stats $aln_mem_sort > $stats
else
  echo "Les fichiers stats existent déjà"
fi

################################################################################
# ELIMINATION DES DUPLICATES DES BAM
################################################################################

# génération des bam en enlevant les duplicats PCR et optiques
if  [ ! -f "$rmduplicates" ];then
  echo "Les fichiers bam sans les duplicats n'existent pas"
java -jar ${path}/Tools/picard.jar MarkDuplicates --INPUT $aln_mem_sort \
--OUTPUT $rmduplicates \
--METRICS_FILE $rmduplicates_metrics \
--REMOVE_DUPLICATES true \
--CREATE_INDEX true \
--TAGGING_POLICY All 2>&1 | tee -a $log_rmduplicates
else
  echo "Les fichiers bam sans les duplicats et les fichiers de métriques existent déjà"
fi

# Stats à partir des bam sans duplicats
if [ ! -f "$flagstat_rmduplicates" ];then
  echo "Les fichiers flagstat des bam sans duplicats n'existent pas"
samtools flagstat $rmduplicates > $flagstat_rmduplicates
else
  echo "Les fichiers flagstat des bam sans duplicats existent déjà"
fi

if [ ! -f "$idxstats_rmduplicates" ];then
  echo "Les fichiers idxstats des bam sans duplicats n'existent pas"
samtools idxstats $rmduplicates > $idxstats_rmduplicates
else
  echo "Les fichiers idxstats existent déjà"
fi

if [ ! -f "$stats_rmduplicates" ];then
  echo "Les fichiers stats des bam sans duplicats n'existent pas"
samtools stats $rmduplicates > $stats_rmduplicates
else
  echo "Les fichiers stats des bam sans duplicats existent déjà"
fi

################################################################################
# RECALIBRATION DU SCORE DE QUALITE DES BASES BQSR
################################################################################
# 3 outils impliqués: BaseRecalibrator, Apply Recalibration, AnalyseCovariates
# BASERECALIBRATOR

if  [ ! -f "$read_groups" ];then
  echo "Les fichiers bam_read_groups n'existent pas"
java -jar ${path}/Tools/picard.jar AddOrReplaceReadGroups I=$rmduplicates O=$read_groups RGLB=capture RGPL=illumina RGPU=stage RGSM=${sample}
else
  echo "Les fichiers bam_read_groups existent déjà"
fi

if  [ ! -f "$before_base_recalibrator" ];then
  echo "Les reports before_bqsr n'existent pas"
gatk BaseRecalibrator --input $read_groups --known-sites $vcf_snp --known-sites $vcf_indel --reference $reference --output $before_base_recalibrator
else
  echo "Les reports before_bqsr existent déjà"
fi

# APPLY RECALIBRATION
if  [ ! -f "$apply_bqsr" ];then
  echo "Les fichiers bam_apply_bqsr n'existent pas"
gatk ApplyBQSR --input $read_groups --reference $reference --bqsr-recal-file $before_base_recalibrator --output $apply_bqsr
else
  echo "Les fichiers bam_apply_bqsr existent déjà"
fi

# BASERECALIBRATOR
if  [ ! -f "$after_base_recalibrator" ];then
  echo "Les reports after_bqsr n'existent pas"
gatk BaseRecalibrator --input $apply_bqsr --known-sites $vcf_snp --known-sites $vcf_indel --reference $reference --output $after_base_recalibrator
else
  echo "Les reports after_bqsr existent déjà"
fi

# ANALYSE COVARIATES
if  [ ! -f "$before_after_bqsr_plot" ];then
  echo "Les fichiers plot_analyze_covariates n'existent pas"
gatk AnalyzeCovariates --before-report-file $before_base_recalibrator --after-report-file $after_base_recalibrator --plots-report-file $before_after_bqsr_plot
else
  echo "Les fichiers plot_analyze_covariates existent déjà"
fi

################################################################################
# HAPLOTYPE CALLER : base quality score = 20
################################################################################
if [ ! -f "$vcf" ];then
  echo "Les fichiers vcf n'existent pas"
gatk HaplotypeCaller --input $apply_bqsr --output $vcf --reference $reference --dbsnp $vcf_snp --bam-output $bamout --min-base-quality-score 20 --intervals $bed
else
  echo "Les fichiers vcf existent déjà"
fi

################################################################################
# CNNScoreVariants 2D et FilterVariantTranches 2D
################################################################################
# CNNSCOREVARIANTS
if [ ! -f "$vcf_cnnscorevariants" ];then
  echo "Les fichiers vcf_cnnscorevariants n'existent pas"
gatk CNNScoreVariants --variant $vcf --reference $reference --output $vcf_cnnscorevariants --intervals $bed --input $apply_bqsr --tensor-type read_tensor
else
  echo "Les fichiers vcf_cnnscorevariants existent déjà"
fi

# FILTERVARIANTTRANCHES
if  [ ! -f "$vcf_filtervarianttranches" ];then
  echo "Les fichiers vcf_filtervarianttranches n'existent pas"
gatk FilterVariantTranches --variant $vcf_cnnscorevariants --resource $vcf_snp --resource $vcf_indel --info-key CNN_2D --snp-tranche 99.95 --indel-tranche 99.4 --output $vcf_filtervarianttranches
else
  echo "Les fichiers vcf_filtervarianttranches existent déjà"
fi

################################################################################
# VT DECOMPOSE et VT NORMALIZE
################################################################################
# VT DECOMPOSE
if [ ! -f "$vcf_decomposed" ];then
  echo "Les fichiers vcf_decomposed n'existent pas"
vt decompose $vcf_filtervarianttranches -o $vcf_decomposed
else
  echo "Les fichiers vcf_decomposed existent déjà"
fi

# VT NORMALIZE
if  [ ! -f "$vcf_normalized" ];then
  echo "Les fichiers vcf_normalized n'existent pas"
vt normalize $vcf_decomposed -q -m -r $reference -o $vcf_normalized
else
  echo "Les fichiers vcf_normalized existent déjà"
fi

################################################################################
# FUNCOTATOR = annotation des variants avec Gencode
################################################################################
if  [ ! -f "$vcf_funcotated" ];then
  echo "Les fichiers vcf_funcotated n'existent pas"
gatk Funcotator --variant $vcf_normalized --reference $reference --ref-version hg19 --data-sources-path $path/Analysis/funcotator_dataSources.v1.4.20180615/gencode \
--output $vcf_funcotated --output-file-format VCF
else
  echo "Les fichiers vcf_funcotated existent déjà"
fi
