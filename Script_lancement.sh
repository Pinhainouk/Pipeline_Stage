#!/bin/bash

fichier_liste="/home/elodie/Documents/Elodie/Pipeline_Stage/Config.txt" # Fichier contenant la liste des echantillons
script="/home/elodie/Documents/Elodie/Pipeline_Stage/Script.sh"

while IFS= read -r ligne
do
  # Vérifier si la ligne n'est pas vide
    if [ -n "$ligne" ]; then
        bash "$script" -o /home/elodie/Documents/ -r /home/elodie/Documents/Genome/hg19.p13.plusMT.no_alt_analysis_set.fa.gz -n /home/elodie/Documents/Genome/dbsnp_138.hg19.vcf \
        -i /home/elodie/Documents/Genome/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -f /home/elodie/Documents/Analysis/funcotator_dataSources.v1.4.20180615 -s "$ligne"
    fi
done < "$fichier_liste"
