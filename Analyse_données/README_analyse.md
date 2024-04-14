# **README : Analyse des données du pipeline**

## **Description**
Ce document détaille l'analyse des fichiers obtenus avec le pipeline et leur comparaison avec le Gold Standard afin de déterminer la performance du pipeline.

## **Les différentes étapes :**
* L'analyse de la couverture est réalisée à partir des BAM avec bamCoverage de deepTools. La couverture est calculée avec le nombre de reads par fenêtre consécutive de 10 nucléotides (bin=10).

* Pour comparer les VCF en sortie du pipeline et le vcf du Gold Standard, **Hap.py** est utilisé (outil d'Illumina) pour déterminer la sensibilité (recall), la spécificité (precision) et le score F1 afin de déterminer si le modèle est performant. Cet outil compare les génotypes au niveau des haplotypes (superlocus de 1 à 1000pb).
Ensuite, une courbe Precision/Recall est réalisé avec R (https://github.com/Illumina/hap.py/blob/master/src/R/rocplot.Rscript).

* R et Rstudio sont utilisés pour l'analyse des vcf dans un Rmarkdown (Analyse.Rmd)

* IGV est utilisé pour visualisé les fichiers bigWig de bamCoverage.

## **Guide d'installation et versions des outils :**
| Outils | Liens d'installation | Versions
|----- |----- |-----
**bamCoverage**     | https://deeptools.readthedocs.io/en/develop/content/installation.html  | **v3.5.4.post1**
**Hap.py**          | https://github.com/Illumina/hap.py?tab=readme-ov-file#installation<BR> | **v0.3.10**
**R**               | https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html#installing-r | **v4.1.2**
**Rstudio**         | https://posit.co/download/rstudio-desktop/                             | **2023.12.0 Build 369**
**IGV**		    | https://data.broadinstitute.org/igv/projects/downloads/2.16/           | **v2.16.2**

## **Problèmes rencontrés :**

### Hap.py (Illumina)
Hap.py utilise une version de Python inférieure à la v3 (v2.7). En installant dans un docker Python 2.7 et Hap.py, des problèmes de mémoires sont rencontrés sur la machine à l'exécution alors la commande a été lancée sur une machine plus puissante.

```
export HGREF=/srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz

hap.py /srv/gold_on_data_elodie_ismael_decomposed_normalize_uniq_copie_chr.vcf.gz \
/srv2/gold_12x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_uniq.vcf.gz \
-r /srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
-o /srv/ref-12x \
-f /srv4/01_data_elodie_V2.bed

hap.py /srv/gold_on_data_elodie_ismael_decomposed_normalize_uniq_copie_chr.vcf.gz \
/srv2/gold_24x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_uniq.vcf.gz \
-r /srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
-o /srv/ref-24x \
-f /srv4/01_data_elodie_V2.bed

hap.py /srv/gold_on_data_elodie_ismael_decomposed_normalize_uniq_copie_chr.vcf.gz \
/srv2/gold_40x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_uniq.vcf.gz \
-r /srv3/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
-o /srv/ref-40x \
-f /srv4/01_data_elodie_V2.bed
```

## **Lancement du Script d'analyse:**

Lancer le Script_analyse.sh sans options.

## **Données de sortie** :
- Fichiers bigWig (BW) de bamCoverage avec piste de couverture.
- Fchiers bedGraph (BD) de bamCoverage avec nombres de reads par intervalle de 10 nucléotides.
- Fichiers JSON, CSV et VCF générés par Hap.py.
- Fichier HTML généré à partir du Rmarkdown avec les résultats d'analyse.

## **Glossaire :**

BAM : Binary Alignment Map

BD : bedGraph

BW : bigWig

CSV : Comma-Separated Values (fichier texte représentant des données tabulaires sous forme de valeurs séparées par des virgules)

IGV : Integrative Genomic Viewer

JSON : JavaScript Object Notation (représentation de données structurées basées sur des paires nom/valeur et des listes ordonnées)

VCF : Variant Call Format
