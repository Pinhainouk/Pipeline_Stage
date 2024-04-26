#!/bin/bash -i
while getopts "o:s:r:n:i:f:a:v:" opt
do
  case "$opt" in
    o) path=${OPTARG};;
  esac
  done
  if [ -z "${path}" ]; then
    echo "Usage: $0 -o path"
    echo "La commande $opt nécessite une option"
    exit 1
  fi
  echo "options fournies:"
  echo "path : ${path}"

################################################################################
# MULTIQC SUR LES RAW
################################################################################

if [ ! -d "${path}/QC_Raw/MultiQC_Raw" ] ;then
  echo "Le multiqc sur les fastqc_raw n'existe pas";
  multiqc ${path}/QC_Raw --module fastqc --filename multiQC_Raw --outdir ${path}/QC_Raw/MultiQC_Raw
else
  echo "Le multiqc sur les fastqc_raw existe déjà"
fi

################################################################################
# MULTIQC SUR LES RAW_TRIM
################################################################################

if [ ! -d "${path}/QC_Trimming/MultiQC_Trimming" ] ;then
  echo "Le multiqc sur les fastqc_trimming n'existe pas";
  multiqc ${path}/QC_Trimming --module fastqc --filename multiQC_Trimming --outdir ${path}/QC_Trimming/MultiQC_Trimming
else
  echo "Le multiqc sur les fastqc_trimming existe déjà"
fi

################################################################################
# MULTIQC SUR STATS FLAGSTAT ET IDXSTATS
################################################################################

if [ ! -d "${path}/QC_Alignment/MultiQC_Alignment" ] ;then
  echo "Le multiqc sur les stats alignement n'existe pas";
  multiqc ${path}/QC_Alignment --module samtools --filename multiQC_Alignement --outdir ${path}/QC_Alignment/MultiQC_Alignment
else
  echo "Le multiqc sur les stats alignement existe déjà"
fi
