
# **README : Pipeline d'analyse de données de séquençage à haut débit pour la détection de variants de type SNP.**

## **Description**
Ce document détaille les différentes étapes de la réalisation d'un pipeline d'analyse de données de séquençage haut débit pour la détection de variants de type SNP. 

## **Les différentes étapes :**
* L'évaluation de la qualité des données issues du séquençage avec l'outil **FastQC** permet de déterminer le nombre de séquences pour chaque échantillon, la qualité des séquences, leur longueur, une éventuelle contamination, la présence résiduelle d'adaptateurs de séquençage.

* Le trimming avec l'outil **BBDuk** permet l'élimination des adaptateurs et des bases de mauvaise qualité. k=27: longueur des sous séquences pour la recherche des adaptateurs. mink=11: longueur minimal pour la recherche des adaptateurs. ktrim=r: pour retirer les adaptateurs à l'extrémité 3'. qtrim=rl: pour retirer les bases de mauvaises qualité des 2 extrémités.  trimq=20: seuil de qualité minimal pour le trimming des bases de mauvaise qualité.

* **Statistique de trimming** : un fichier de statistiques détaillées sur le trimming et le filtrage des séquences.

* Le contrôle de la qualité des fichiers de séquences nettoyées avec **FastQC**.

* L'alignement avec **BWA-MEM** consiste à trouver des correspondances exactes maximales de sous-séquences dans le génome de référence indexé et de procéder à une extension de ces sous-séquences avec l'algorithme d'alignement local de Smith-Waterman pour obtenir des alignements plus précis. **BWA-MEM** permet d'aligner des séquences requêtes de 70 à 1Mbp.

* Réalisation de statistiques à partir des bam avec **Samtools**.

* L'élimination des duplicats de PCR liés aux étapes de PCR durant la préparation des librairies et des duplicats optiques liés au séquençage avec l'outil **MarkDuplicates de Picard**. 
* Fichier de métriques après enlèvement des duplicats de PCR et optiques.

* Réalisation de statistiques à partir des bam sans duplicats de PCR avec **Samtools**.

* La recalibration du niveau de qualité des bases avec l'outil **BaseQualityScoreRecalibration (BQSR) de GATK**. Il consiste à corriger tout biais systématiques (sur-estimations ou sous-estimations des scores de qualité) afin d'améliorer la précision des scores de qualité de base. Cette étape à pour but de réduire les erreurs de séquençage et d'amélorer la fiabilité des appels de variants ultérieurs. 

* L'appel des variants avec l'outil **Haplotype Caller de GATK**. L'algorithme détermine des régions actives (celles présentant des variations), identifie les haplotypes possibles, ré-aligne chaque haplotype par rapport à l'haplotype de référence à l'aide de Smith-Waterman afin de détecter de potentielles variations.

* L'annotation du fichier .vcf avec **CNNScoreVariants de GATK Model 2D** pour attribuer des scores de confiance aux variants génétiques en se basant sur l'apprentissage profond fourni par un réseau de neurones convolutionnel (CNN). Ces scores aident à évaluer la probabilité qu'un variant soit réel par opposition à un faux positif.

* L'application d'un filtre par tranche avec **FilterVariantTranches de GATK Model 2D** aux fichiers vcf, en fonction des scores provenant de l'étape précédente **CNNScoreVariants de GATK**, permet de limiter les faux variants. --snp-tranche 99.95: filtre par tranche 99.95% des snp
--indel-tranche 99.4 : filtre par tranche 99.4% des indel en fonction des scores provenant de l'annotation dans le champ INFO de l'étape précédente.

* **Vt decompose** est utilisé pour décomposer un variant complexe représenté par plusieurs variants simples (une insertion suivie d'une substitution décomposée en deux variants distinces sur deux lignes du vcf). -s divise les champs INFO et GENOTYPE de façon appropriée.

* **Vt normalize** est utilisé pour normaliser les variants (qui ne sont pas toujours représenter de la même façon dans le vcf selon l'outil d'appel de variants.) La normalisation ajuste les positions de début et de fin de la représentation des variants. -q n'imprime pas les ptions ni le résumé -m donne des avertissements mais ne quitte pas lorsque REF est incompatible avec la séquence de référence masquée pour les non-SNP. 

* **Vt uniq** est utilisé pour supprimer les éventuels variants en double dans le vcf. 

* La dernière étape est celle de l'annotation fonctionnelle des variants à l'aide de **Funcotator de GATK**. Il utilise une ou plusieurs bases de données mises à jour régulièrement qui contiennent des informations sur les régions génomiques, les transcrits, les protéines et les effets fonctionnels des variants sur les gènes associés. Nous avons utilisé la base Gencode (funcotator_dataSources.v1.4.20180615/gencode)

## **Pré-requis à l'utilisation du pipeline :**

* **Installation de Python (v3.11.5) et java (v17.0.9) pour le fonctionnement des outils.**


* **Création de l'index du génome hg19: env. 2-3h**

L'indexation du génome de référence avec **BWA** permet de préparer une structure de données qui facilite la recherche rapide des positions des séquences dans le génome à l'étape d'alignement.
Il s'agit du génome version hg19 sans les chromosomes alternatifs et avec les mitochrondries.
**hg19.p13.plusMT.no_alt_analysis_set.fa.gz** téléchargé sur Genome Browser (dernière modification 2020-03-09 10:21) : 
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/

**Quand on créé l'index, il est important qu'il ait le même nom que la référence.**

```
./bwa index -p hg19.p13.plusMT.no_alt_analysis_set.fa {path}/hg19.p13.plusMT.no_alt_analysis_set.fa.gz
```

* **Installation de l'environnement conda (v23.11.0) pour l'outil GATK.**

* **Création des alias des outils :**
```
vim ~/.bashrc
```

| Outils | Noms alias |
|----- |------
| **FastQC** | fastqc
| **BBDuk** | bbduk
| **BWA** | bwa
| **Samtools** | samtools
| **Picard** | picard
| **GATK** | gatk
| **Vt** | vt

## **Guide d'installation et versions des outils :** 
| Outils | Liens d'installation | Versions
|----- |----- |-----
| **FastQC** | https://www.bioinformatics.babraham.ac.uk/projects/download.html#FastQC | **v0.12.1**
| **BBDuk**          | https://sourceforge.net/projects/bbmap/ | **v39.03**
| **BWA**            | https://github.com/lh3/bwa/releases | **v0.7.17-r1188**
| **Samtools**       | https://sourceforge.net/projects/samtools/files/samtools/1.3.1/ | **v1.3.1**
| **Picard**         | https://broadinstitute.github.io/Picard/ | **v3.1.1**
| **GATK**           | S'assurer d'avoir la version java 17.<BR>Installer GATK : https://github.com/broadinstitute/gatk/releases<BR>Dans le dossier gatk-4.5.0.0, pour créer un environnement GATK avec toutes les dépendances, lancer la commande :<BR>```conda env create -f gatkcondaenv.yml```<BR>Activer l'environnement avec : <BR>```conda activate gatk``` | **v4.5.0.0**
| **Vt**             | https://genome.sph.umich.edu/wiki/Vt#Installation<BR>Bien suivre l'installation wiki via le git clone pour avoir les options requises.| **v0.5**
| **Bcftools**       | https://github.com/samtools/bcftools/releases?page=2 | **v1.10.2**
| **MultiQC**        | https://multiqc.info/docs/getting_started/installation/ | **v1.19**

## **Problèmes rencontrés :** 
### GATK v4.0 
La nouvelle version 4.5.0.0 de GATK installée dans l'environnement conda résoud les erreurs suivantes rencontrées à l'exécution de CNNScoreVariants Model 2D avec les versions v4.0 antérieures :
```
java.lang.RuntimeException: A required Python package ("gatktool") could not be imported into the Python environment. This tool requires that the GATK Python environment is properly established and activated. Please refer to GATK [README.md](http://readme.md/) file for instructions on setting up the GATK Python environment.
```
```
org.broadinstitute.hellbender.utils.python.PythonScriptExecutorException: A nack was received from the Python process (most likely caused by a raised exception caused by): nck received
```

### Vt v0.5
Vt installé sans le git clone renvoie une erreur à l'exécution de Vt normalize mais pas à l'exécution de Vt decompose :
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
Paramètres requis à la bonne exécution de Vt normalize, obtenus en suivant l'installation wiki avec le git clone :
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
## **Lancement du pipeline :**

**Avant d'exécuter le Script, activer l'environnement conda pour le fonctionnement de l'outil gatk :**

```
conda activate gatk
```

**S'assurer d'avoir les fichiers et dossiers suivants pour le fonctionnement des outils :**

- Le fichier bed
- La fichier référence du génome
- La fichier des adaptateurs pour BBDUK
- le fichier vcf avec les SNP connus
- le fichier vcf avec les indels connues
- le dossier package des bases de données Funcotator

**Deux façons de lancer le script du pipeline :**

* Soit on lance le **Script.sh** avec les **options obligatoires suivantes** : 

```
-o : le chemin où vont être rangés les fichiers d'entrée et de sortie du pipeline
-s : le nom du ou des échantillons à analyser
-r : le chemin complet de la référence du génome
-n : le chemin complet du fichier vcf des snp connus
-i : le chemin complet du fichier vcf des indels connues
-f : le chemin complet du package funcotator contenant les bases de données
-a : le chemin complet du fichier contenant les adaptateurs
-v : saisir la version du génome utilisée (hg19 ou hg38)
```
Exemple de lancement :

```
./Script.sh -o /home/elodie/Documents/ \
-s gold_12x_on_data_elodie \
-r /home/elodie/Documents/Genome/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
-n /home/elodie/Documents/Genome/dbsnp_138.hg19.vcf \
-i /home/elodie/Documents/Genome/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-f /home/elodie/Documents/Analysis/funcotator_dataSources.v1.4.20180615 \
-a /home/elodie/Documents/Elodie/adapters.fa \
-v hg19 \
```

* Soit on lance le **Script_lancement.sh** avec un fichier de config (Config.txt) contenant tous les échantillons à analyser. 
Dans le Script_lancement.sh, penser à modifier les chemins dans les variables ```fichier_liste``` et ```script``` pour indiquer les chemins où se trouvent respectivement le fichier de config et le script. Il boucle alors sur le fichier de config tant qu'une ligne n'est pas vide et il exécute le script sur chaque échantillon.

Exemple de fichier de config à faire avec une ligne par échantillon :

```
nom_echantillon_1
nom_echantillon_2
nom_echantillon_3
...
```

Exemple de lancement :

```
./Script_lancement.sh
```

## **Fichiers de sortie** :
- Fichiers BAM finaux, après l'étape BQSR : ```{path}/{sample}_aln_mem_sort_rmduplicates_apply_bqsr.bam```
- Fichiers VCF finaux, après l'annotation fonctionnelle Funcotator : ```{path}/{sample}_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf```
- MultiQC : des fastq, fastq après bbduk et stats samtools

## **Glossaire :**
BAM : Binary Alignment Map

BBDuk : DUK pour Decontamination Using Kmers (de la suite d'outils BBTools).

BQSR : Base Quality Score Recalibration

BWA : Burrows-Wheeler Aligner

duplicats de PCR : paires de lectures qui ont le même début et la même fin d'alignement.

Funcotator : FUNCtional annOTATOR

GATK : Genome Analysis ToolKit

hg19 : human genome GRCh37

MEM : Maximal Exact Matches

NGS : Next Generation Sequencing

PCR : Polymerase Chain Reaction

read : fragment de plusieurs bases générés et lu par un séquenceur appelé aussi une lecture, une séquence

Réseau de neurones convolutionnel : système dont la conception est à l'origine schématiquement inspirée du fonctionnement des neurones biologiques, et qui par la suite s'est rapproché des méthodes statistiques. Lors de l'entraînement, un ensemble de données est utilisé pour permettre au modèle d'apprendre et de reconnaître un objet ou quelque chose en particulier.

SNP : Single Nucleotide Polymorphism, variation d'une seule paire de base entre individus d'une même espèce

VCF : Variant Call Format

Vt : Variant Tool





