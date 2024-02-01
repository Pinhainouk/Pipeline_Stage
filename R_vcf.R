vcf_12x = read.table("/home/elodie/Documents/Analysis/gold_12x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf")
vcf_24x = read.table("/home/elodie/Documents/Analysis/gold_24x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf")
vcf_40x = read.table("/home/elodie/Documents/Analysis/gold_40x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf")
vcf_ref = read.table("/home/elodie/Documents/Elodie/Datas_ismael_ref/gold_on_data_elodie_ismael_decomposed_normalize.vcf.gz")

#Sélection des colonnes 1,2,4,5
#library('dplyr)
vcf_12x_select =  vcf_12x %>% select(1,2,4,5)
vcf_24x_select =  vcf_24x %>% select(1,2,4,5)
vcf_40x_select =  vcf_40x %>% select(1,2,4,5)
vcf_ref_select =  vcf_ref %>% select(1,2,4,5)

#Changement de "1" en "chr1" du vcf ref
vcf_ref_select$V1 = sub("1", "chr1", vcf_ref_select$V1)
vcf_ref_select$V1 = sub("4", "chr4", vcf_ref_select$V1)
vcf_ref_select$V1 = sub("8", "chr8", vcf_ref_select$V1)

library('tidyr')
#Concatener les colonnes 1,2,4,5
vcf_12x_concat = unite(vcf_12x_select, Variants, sep = "_", remove = TRUE)
vcf_24x_concat = unite(vcf_24x_select, Variants, sep = "_", remove = TRUE)
vcf_40x_concat = unite(vcf_40x_select, Variants, sep = "_", remove = TRUE)
vcf_ref_concat = unite(vcf_ref_select, Variants, sep = "_", remove = TRUE)
  
#Ajout de l'ID dans les 4 vcf
vcf_12x_concat <- vcf_12x_concat %>% mutate(Id = 1:n())
vcf_24x_concat <- vcf_24x_concat %>% mutate(Id = 1:n())
vcf_40x_concat <- vcf_40x_concat %>% mutate(Id = 1:n())
vcf_ref_concat <- vcf_ref_concat %>% mutate(Id = 1:n())

#Merger les 4 vcf par les variants
vcf_merge <- merge(vcf_12x_concat, vcf_24x_concat, by = "Variants", all.x=TRUE, all.y=TRUE)
names(vcf_merge)[2]= "vcf_12x"
names(vcf_merge)[3]= "vcf_24x"

vcf_merge = merge(vcf_merge, vcf_40x_concat, by = "Variants", all.x=TRUE, all.y=TRUE)
names(vcf_merge)[4]= "vcf_40x"

vcf_merge = merge(vcf_merge, vcf_ref_concat, by = "Variants", all.x=TRUE, all.y=TRUE)
names(vcf_merge)[5]= "vcf_ref"

#Remplacer les NA par 0 et le reste par 1
vcf_merge$vcf_12x <- ifelse(is.na(vcf_merge$vcf_12x), 0, 1)
vcf_merge$vcf_24x <- ifelse(is.na(vcf_merge$vcf_24x), 0, 1)
vcf_merge$vcf_40x <- ifelse(is.na(vcf_merge$vcf_40x), 0, 1)
vcf_merge$vcf_ref <- ifelse(is.na(vcf_merge$vcf_ref), 0, 1)

#Somme des occurrences des variants dans les vcf
vcf_merge$Somme <- rowSums(vcf_merge[, c("vcf_12x", "vcf_24x", "vcf_40x", "vcf_ref")])

#library('VennDiagram')
#library('grid')
#library('futile.logger')

ensemble=list(vcf_12x=vcf_12x_concat$Variants, vcf_24x=vcf_24x_concat$Variants, 
          vcf_40x=vcf_40x_concat$Variants, vcf_ref=vcf_ref_concat$Variants)
          
venn.diagram(ensembl, filename = "Venn Diagramme n&b")

venn.diagram(ensemble, category.names=c("vcf_12x", "vcf_24x", "vcf_40x", "vcf_ref"), filename = "Venn Diagramme couleur", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

display_venn = function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
  
display_venn(ensemble, category.names= c("vcf_12x", "vcf_24x", "vcf_40x", "vcf_ref"), 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))


