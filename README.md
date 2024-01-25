
# **README : Pipeline d'analyse de données de séquençage à haut débit pour la détection de variants de type SNP.**

## **Description**
Ce document détaille les différentes étapes de la réalisation d'un pipeline d'analyse de données de séquençage haut débit pour la détection de variants de type SNP. 

## **Les différentes étapes:**
* L'évaluation de la qualité des données issues du séquençage avec l'outil **FastQC** permet de déterminer le nombre de séquences pour chaque échantillon, la qualité des séquences, leur longueur, une éventuelle contamination, la présence résiduelle d'adaptateurs de séquençage.

* Le trimming avec l'outil **BBDuk** permet l'élimination des adaptateurs et des bases de mauvaise qualité.

* Le contrôle de la qualité des fichiers de séquences nettoyées avec **FastQC**.

* L'indexation du génome de référence (hg19) avec **BWA** qui permet de préparer une structure de données qui facilite la recherche rapide des positions des séquences dans le génome à l'étape suivante d'alignement.

* L'alignement avec **BWA-MEM** consiste à trouver des correspondances exactes maximales de sous-séquences dans le génome de référence indexé et de procéder à une extension de ces sous-séquences avec l'algorithme d'alignement local de Smith-Waterman pour obtenir des alignements plus précis. **BWA-MEM** permet d'aligner des séquences requêtes de 70 à 1Mbp.

* L'élimination des duplicats de PCR liés aux étapes de PCR durant la préparation des librairies et des duplicats optiques liés au séquençage avec l'outil **MarkDuplicates de Picard**. 

* La recalibration du niveau de qualité des bases avec l'outil **BaseQualityScoreRecalibration (BQSR) de GATK**. Il consiste à corriger tout biais systématiques (sur-estimations ou sous-estimations des scores de qualité) afin d'améliorer la précision des scores de qualité de base. Cette étape à pour but de réduire les erreurs de séquençage et d'amélorer la fiabilité des appels de variants ultérieurs. 

* L'appel des variants avec l'outil **Haplotype Caller de GATK**. L'algorithme détermine des régions actives (celles présentant des variations), identifie les haplotypes possibles, ré-aligne chaque haplotype par rapport à l'haplotype de référence à l'aide de Smith-Waterman afin de détecter de potentielles variations.

* L'annotation du fichier .vcf avec **CNNScoreVariants de GATK** pour attribuer des scores de confiance aux variants génétiques en se basant sur l'apprentissage profond fourni par un réseau de neurones convolutionnel (CNN). Ces scores aident à évaluer la probabilité qu'un variant soit réel par opposition à un faux positif.

* L'application d'un filtre par tranche avec **FilterVariantTranches de GATK** aus fichiers VCF, en fonction des scores provenant de l'étape précédente **CNNScoreVariants de GATK**, permet de limiter les faux variants.

* **Vt decompose** est utilisé pour décomposer un variant complexe représenté par plusieurs variants simples (une insertion suivie d'une substitution décomposée en deux variants distinces sur deux lignes du vcf).

* **Vt normalize** est utilisé pour normaliser les variants (pas toujours représenter de la même façon dans le vcf selon l'outil d'appel de variants.) La normalisation ajuste les positions de début et de fin de la représentation des variants.

* La dernière étape est celle de l'annotation fonctionnelle des variants à l'aide de **Funcotator de GATK**. Il utilise une ou plusieurs bases de données mises à jour régulièrement qui contiennent des informations sur les régions génomiques, les transcrits, les protéines et les effets fonctionnels des variants sur les gènes associés.  

## **Pré-requis à l'utilisation du pipeline:**
Installation de Python (v3.11.5) et java (v17.0.9) pour le fonctionnement des outils.

Installation de l'environnement conda (v23.11.0) pour l'outil GATK.

Avant d'exécuter le Script, toujours être dans l'environnement conda:

```
conda activate gatk
```
Des alias des outils ont été créés
```
vim ~/.bashrc
```

## **Guide d'installation:** 
| Outils            | Liens d'installation
|--------           |-----------------------------------------------------------------------
|**FastQC**         |https://www.bioinformatics.babraham.ac.uk/projects/download.html#FastQC
|**BBDuk**          |https://sourceforge.net/projects/bbmap/
|**BWA**            |https://github.com/lh3/bwa/releases
|**Samtools**       |https://sourceforge.net/projects/samtools/files/samtools/1.3.1/
|**Picard**         |https://broadinstitute.github.io/Picard/
|**GATK**           |https://github.com/broadinstitute/gatk/releases
|**Vt**             |https://genome.sph.umich.edu/wiki/Vt#Installation
|**Bcftools**       |https://github.com/samtools/bcftools/releases?page=2
|**MultiQC**        |https://multiqc.info/docs/getting_started/installation/

## **Guide d'utilisation - étapes et paramètres utilisés**

### FastQC sur les données brutes:
```
fastqc R1_file R2_file -o path_output 2>&1 | tee -a log_qc_raw
```

### BBDuk : 
```
bbduk in1=R1_file in2=R2_file k=27 mink=11 ktrim=r qtrim=rl trimq=20 out1=R1_file_trim out2=R2_file_trim stats=stats_trim tpe 2>&1 | tee -a log_trim
```

### FastQC sur les données nettoyées:
```
fastqc R1_file_trim R2_file_trim -o path_output 2>&1 | tee -a log_qc_trim
```

### Création de l'index du génome: env. 2-3h
Génome hg19_genome_index path_hg19.p13.plusMT.no_alt_analysis_set.fa.gz téléchargé sur Genome Browser : https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/
```
./bwa index -p hg19_genome_index hg19_genome_file.fa.gz
```

### Alignement :
#### 1- Génération de fichiers sam:
 ```
 bwa mem hg19_genome_index R1_file_trim R2_file_trim -o file_alignment 2>&1 | tee -a log_alignment
```
#### 2- Génération de fichiers bam, tri des bam et création d'index des bam:
 ```
 samtools view -b file_alignment
 samtools sort -o file_alignment_sort
samtools index file_alignment_sort
```

### Retrait des duplicats de PCR et optiques:
```
java -jar picard.jar MarkDuplicates --INPUT file_alignment_sort \
--OUTPUT rmduplicates \
--METRICS_FILE file_rmduplicates_metrics \
--REMOVE_DUPLICATES true \
--CREATE_INDEX true \
--TAGGING_POLICY All 2>&1 | tee -a log_rmduplicates
```

### Pré-requis nécessaires à la commande BQSR:
#### 1- Création d'un index .fai du génome avec Samtools
Le genome en gzip et non en bgzip requis par samtools
```
samtools faidx hg19_genome_file.fa.gz
```
#### 2- Création d'un dictionnaire .dict avec Picard
```
java -jar picard.jar CreateSequenceDictionary R=/home/elodie/Documents/Genome/hg19.p13.plusMT.no_alt_analysis_set.fa.gz O=/home/elodie/Documents/Genome/hg19.p13.plusMT.no_alt_analysis_set.dict
```

Dans les fichiers **dbsnp_138.b37.vcf** et **Mills_and_1000G_gold_standard.indels.b37.vcf**, le chromosome est noté "1" et dans la référence du génome le chromosome 1 est noté "chr1" donc on remplace dans les fichiers vcf le "1" par "chr1" : 
Tous les débuts de ligne commençant par  "#" (header) sont ignorées puis il remplace le "1" par "chr1" des lignes suivantes.
```
sed -e '/^#/!s/^/chr/' /path/dbsnp_138.b37.vcf > /path/dbsnp_138.b37_chr.vcf
sed -e '/^#/!s/^/chr/' /path/Mills_and_1000G_gold_standard.indels.b37.vcf > /path/Mills_and_1000G_gold_standard.indels.b37_chr.vcf
```
#### 3- Création de fichiers bam avec readsGroup
```
java -jar picard.jar AddOrReplaceReadGroups I=file_alignment_sort_rmduplicates.bam O=file_alignment_sort_rmduplicates_readGroup.bam RGID=4 RGLB=capture RGPL=illumina RGPU=psl RGSM=gold_12x
```

### Base Quality Score Recalibration:
#### 1- Rapport avant application de BQSR
```
gatk BaseRecalibrator --input file_alignment_sort_rmduplicates_readGroup.bam  --known-sites dbsnp_138.b37.vcf --known-sites Mills_and_1000G_gold_standard.indels.b37.vcf --reference hg19_genome_file.fa.gz --output file_before_base_recalibrator
```
#### 2- Application de BQSR : obtention de bam recalibrés
```
gatk ApplyBQSR --input file_alignment_sort_rmduplicates_readGroup.bam --reference hg19_genome_file.fa.gz --bqsr-recal-file file_before_base_recalibrator --output apply_bqsr.bam
```
#### 3- Rapport après application de BQSR
```
gatk BaseRecalibrator --input apply_bqsr.bam --known-sites dbsnp_138.b37.vcf --known-sites Mills_and_1000G_gold_standard.indels.b37.vcf --reference hg19_genome_file.fa.gz --output file_after_base_recalibrator
```
#### 4- Rapport d'analyse de la covariation avant-après BQSR
```
gatk AnalyzeCovariates --before-report-file file_before_base_recalibrator --after-report-file file_after_base_recalibrator --plots-report-file file_before_after_bqsr_plot
```

### Appel des variants Haplotype Caller
```
gatk HaplotypeCaller --input apply_bqsr.bam --output file_vcf --reference hg19_genome_file.fa.gz --dbsnp dbsnp_138.b37.vcf --bam-output bamout --min-base-quality-score 20 --intervals file_bed
```

### CNNScoreVariants - Model 2D
```
gatk CNNScoreVariants --variant file_vcf --reference hg19_genome_file.fa.gz --output file_vcf_cnnscorevariants --intervals file_bed --input apply_bqsr.bam --tensor-type read_tensor
```

### FilterVariantTranches - Model 2D
```
gatk FilterVariantTranches --variant file_vcf_cnnscorevariants --resource dbsnp_138.b37.vcf --resource Mills_and_1000G_gold_standard.indels.b37.vcf --info-key CNN_2D --snp-tranche 99.95 --indel-tranche 99.4 --output file_vcf_filtervarianttranches
```

### Vt decompose
```
vt decompose file_vcf_filtervarianttranches -o file_vcf_filtervarianttranches_decomposed
```

### Vt normalize
```
vt normalize file_vcf_filtervarianttranches_decomposed -q -m -r hg19_genome_file.fa.gz -o file_vcf_filtervarianttranches_decomposed_normalized
```

### Funcotator
```
gatk Funcotator --variant file_vcf_filtervarianttranches_decomposed_normalized --reference hg19_genome_file.fa.gz --ref-version hg19 --data-sources-path $path/funcotator_dataSources.v1.4.20180615/gencode --output file_vcf_funcotated --output-file-format VCF
```

## **Problèmes rencontrés:** 
### GATK v4.0 
La nouvelle version 4.5.0.0 de GATK résoud les erreurs suivantes rencontrées à l'exécution de CNNScoreVariants avec les versions v4.0 antérieures:
```
java.lang.RuntimeException: A required Python package ("gatktool") could not be imported into the Python environment. This tool requires that the GATK Python environment is properly established and activated. Please refer to GATK [README.md](http://readme.md/) file for instructions on setting up the GATK Python environment.
```
```
org.broadinstitute.hellbender.utils.python.PythonScriptExecutorException: A nack was received from the Python process (most likely caused by a raised exception caused by): nck received
```

### Vt v0.5
Vt installé sans le git clone renvoie une erreur à l'exécution de Vt normalize mais pas à l'exécution de Vt decompose:
```
normalize v0.5

options:     input VCF file 
         [o] output VCF file                             
         [w] sorting window size                   
         [n] no fail on reference inconsistency for non SNPs
         [q] quiet
         [d] debug
         [r] reference FASTA file    

Exception en point flottant (core dumped)
```
Paramètres requis à la bonne exécution de Vt normalize:
```
normalize v0.5

description : normalizes variants in a VCF file.

usage : vt normalize [options] <in.vcf>

options : -o  output VCF file [-]
          -d  debug [false]
          -q  do not print options and summary [false]
          -m  warns but does not exit when REF is inconsistent
              with masked reference sequence for non SNPs.
              This overides the -n option [false]
          -n  warns but does not exit when REF is inconsistent
              with reference sequence for non SNPs [false]
          -f  filter expression []
          -w  window size for local sorting of variants [10000]
          -I  file containing list of intervals []
          -i  intervals []
          -r  reference sequence fasta file []
          -?  displays help
 ```         

## **Versions des outils utilisés:**
| Outils       | Versions
|--------      |----------
|**FastQC**    |v0.12.1
|**BBDuk**     |v39.03
|**BWA**       |v0.7.17-r1188
|**Samtools**  |v1.3.1
|**Picard**    |v3.1.1 
|**GATK**      |v4.5.0.0 
|**Vt**        |v0.5  
|**Bcftools**  |v1.10.2  
|**MultiQC**   |v1.19

## **Statistiques**
Stats à partir des bam triés
Stats à partir des bam sans duplicats

## **Glossaire:**
SNP: Single Nucleotide Polymorphism, variation d'une seule paire de base entre individus d'une même espèce

BQSR: Base Quality Score Recalibration

read: fragment de plusieurs bases générés et lu par un séquenceur appelé aussi une lecture, une séquence

NGS: Next Generation Sequencing

BBDuk: DUK pour Decontamination Using Kmers (de la suite d'outils BBTools).

hg19: human genome GRCh37

BWA: Burrows-Wheeler Aligner

MEM: Maximal Exact Matches

PCR: Polymerase Chain Reaction

duplicats de PCR: paires de lectures qui ont le même début et la même fin d'alignement.

GATK: Genome Analysis ToolKit

Réseau de neurones convolutionnel: système dont la conception est à l'origine schématiquement inspirée du fonctionnement des neurones biologiques, et qui par la suite s'est rapproché des méthodes statistiques. Lors de l'entraînement, un ensemble de données est utilisé pour permettre au modèle d'apprendre et de reconnaître un objet ou quelque chose en particulier.

fichiers VCF: Variant Call Format

Funcotator: FUNCtional annOTATOR

Vt: Variant Tool






