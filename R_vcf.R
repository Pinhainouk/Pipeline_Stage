vcf_12x = read.table("/home/elodie/Documents/Analysis/gold_12x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf")
vcf_24x = read.table("/home/elodie/Documents/Analysis/gold_24x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf")
vcf_40x = read.table("/home/elodie/Documents/Analysis/gold_40x_on_data_elodie_filtervarianttranches_2D_decomposed_normalized_funcotated.vcf")
vcf_ref = read.table("/home/elodie/Documents/Elodie/Datas_ismael_ref/gold_on_data_elodie_ismael_decomposed_normalize.vcf.gz")

#Sélection des colonnes 1,2,4,5
library('dplyr')
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

library('VennDiagram')
library('grid')
library('futile.logger')

# VennDiagram n&b
ensemble = list(vcf_12x=vcf_12x_concat$Variants, vcf_24x=vcf_24x_concat$Variants, vcf_40x=vcf_40x_concat$Variants, vcf_ref=vcf_ref_concat$Variants)
venn.diagram(ensemble, category.names=c("vcf_12x" , "vcf_24x" , "vcf_40x", "vcf_ref"), filename = "Venn diagramme vcf_n&b")

# VennDiagram couleur
venn.diagram(ensemble, category.names=c("vcf_12x" , "vcf_24x" , "vcf_40x", "vcf_ref"), filename = "Venn diagramme vcf_couleur", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))


# Fonction d'aide pour afficher le diagramme de Venn
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(ensemble)

display_venn(
  ensemble,
  category.names = c("vcf_12x" , "vcf_24x" , "vcf_40x", "vcf_ref"),
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

 # 01022024
library ("ggVennDiagram")
ggVennDiagram(ensemble, label_alpha = 0, filename = "Venn diagramme vcf_gradient")

library("gplots")
venn.plot = venn.diagram(ensemble, category.names=c("vcf_12x" , "vcf_24x" , "vcf_40x", "vcf_ref"), filename = "Venn diagramme vcf_couleur", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
venndiag = venn(data=ensemble)

cat_ref = as.data.frame(attributes(venndiag)$intersections$vcf_ref)
cat_40x = as.data.frame(attributes(venndiag)$intersections$vcf_40x)
cat_24x = as.data.frame(attributes(venndiag)$intersections$vcf_24x)
cat_12x = as.data.frame(attributes(venndiag)$intersections$vcf_12x)
cat_715 = as.data.frame(attributes(venndiag)$intersections$`vcf_12x:vcf_24x:vcf_40x:vcf_ref`)
cat_24_40_ref = as.data.frame(attributes(venndiag)$intersections$`vcf_24x:vcf_40x:vcf_ref`)
cat_12_40_ref = as.data.frame(attributes(venndiag)$intersections$`vcf_12x:vcf_40x:vcf_ref`)
cat_12_24_40 = as.data.frame(attributes(venndiag)$intersections$`vcf_12x:vcf_24x:vcf_40x`)
cat_40_ref = as.data.frame(attributes(venndiag)$intersections$`vcf_40x:vcf_ref`)
cat_24_ref = as.data.frame(attributes(venndiag)$intersections$`vcf_24x:vcf_ref`)
cat_24_40 = as.data.frame(attributes(venndiag)$intersections$`vcf_24x:vcf_40x`)
cat_12_40 = as.data.frame(attributes(venndiag)$intersections$`vcf_12x:vcf_40x`)
cat_ref$categorie = c("cat_ref")
cat_40x$categorie = c("cat_40x")
cat_24x$categorie = c("cat_24x")
cat_12x$categorie = c("cat_12x")
cat_715$categorie = c("cat_715")
cat_24_40_ref$categorie = c("cat_24_40_ref")
cat_12_40_ref$categorie = c("cat_12_40_ref")
cat_12_24_40$categorie = c("cat_12_24_40")
cat_40_ref$categorie = c("cat_40_ref")
cat_24_ref$categorie = c("cat_24_ref")
cat_24_40$categorie = c("cat_24_40")
cat_12_40$categorie = c("cat_12_40")

list_df = list(cat_ref,cat_40x,cat_24x,cat_12x,cat_715,cat_24_40_ref,cat_12_40_ref,cat_12_24_40,cat_40_ref,cat_24_ref,cat_24_40,cat_12_40)
newcol=c("Variants","Catégorie")
new_list = lapply(list_df, setNames, newcol)
cat=bind_rows(new_list)

cat_split = cat %>% separate(Variants, c("chromosome", "position", "ref", "alt"), "_")

cat_split$type <- ifelse(cat_split$alt == "*", "undetermined",
                        ifelse (nchar(cat_split$ref)==1 & nchar(cat_split$alt)==1, "SNP",
                               ifelse(nchar(cat_split$ref) >= 2 & nchar(cat_split$alt) == 1, "délétion",
                                      ifelse(nchar(cat_split$ref) ==1 & nchar(cat_split$alt)>=2, "insertion", 
                                             "autre"))))
cat_split$gene = ifelse(cat_split$chromosome =="chr4", "EPHA5",
                        ifelse(cat_split$chromosome =="chr1", "UTS2",
                               ifelse(cat_split$chromosome =="chr8", "LPL", "autre")))

cat_split %>% count(cat_split$type,cat_split$gene)
cat_split %>% group_by (cat_split$Catégorie) %>% count(cat_split$type)

fig = ggplot(cat_split, aes(x=cat_split$type)) +
  geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7)
           )+facet_grid(.~cat_split$gene)+theme_classic()

ggplot(cat_split, aes(x=cat_split$gene)) +
  geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7))

ggplot(cat_split, aes(x=cat_split$gene, y=)) +
  geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7))

ggplotly(fig)
