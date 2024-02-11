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
vcf_12x_concat = tidyr::unite(vcf_12x_select, Variants, sep = "_", remove = TRUE)
vcf_24x_concat = tidyr::unite(vcf_24x_select, Variants, sep = "_", remove = TRUE)
vcf_40x_concat = tidyr::unite(vcf_40x_select, Variants, sep = "_", remove = TRUE)
vcf_ref_concat = tidyr::unite(vcf_ref_select, Variants, sep = "_", remove = TRUE)
  
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

cat_split %>% filter(type == "undetermined")

# 08/02/2024

library("ggplot2")
#create datasets
category = cat_split$Catégorie
condition= cat_split$type
nombre = cat_split %>% group_by (cat_split$Catégorie) %>% count(cat_split$type)

data=data.frame(nombre)
plot = ggplot(data, aes(fill=cat_split.type, y=n, x=cat_split.Catégorie))+geom_bar(position = "dodge", stat="identity")
plot = ggplot(data, aes(fill=cat_split.type, y=n, x=cat_split.Catégorie))+geom_bar(position = "stack", stat="identity")

library(plotly)
ggplotly(plot)

vcf_12x_concat$VCF = "vcf12x"
vcf_24x_concat$VCF = "vcf24x"
vcf_40x_concat$VCF = "vcf40x"
vcf_ref_concat$VCF = "vcfref"

vcf_bind = rbind(vcf_12x_concat,vcf_24x_concat,vcf_40x_concat,vcf_ref_concat)
vcf_bind = vcf_bind %>% separate(Variants, c("chromosome", "position", "ref", "alt"), "_")

vcf_bind$Type = ifelse(vcf_bind$alt == "*", "undetermined",
                       ifelse (nchar(vcf_bind$ref)==1 & nchar(vcf_bind$alt)==1, "SNP",
                               ifelse(nchar(vcf_bind$ref) >= 2 & nchar(vcf_bind$alt) == 1, "délétion",
                                      ifelse(nchar(vcf_bind$ref) ==1 & nchar(vcf_bind$alt)>=2, "insertion", 
                                             "autre"))))

vcf_bind$Gene = ifelse(vcf_bind$chromosome == "chr1", "UTS2",
                       ifelse(vcf_bind$chromosome =="chr4", "EPHA5",
                              ifelse(vcf_bind$chromosome =="chr8", "LPL","autre")))

vcf_bind %>% filter(Type == "undetermined")

catvcf = vcf_bind$VCF
condtype = vcf_bind$Type
nb = vcf_bind %>% group_by (vcf_bind$VCF) %>% count(vcf_bind$Type)
data2=data_frame(nb)

# renommer colonnes
names(data2)[1]="Nom_VCF"
names(data2)[2]="Type_Variants"
names(data2)[3]="Nombre_Variants"
plot2 = ggplot(data2, aes(fill=Type_Variants, y=Nombre_Variants, x=Nom_VCF))+geom_bar(position = "dodge", stat="identity")
plot3 = ggplot(data2, aes(fill=Type_Variants, y=Nombre_Variants, x=Nom_VCF))+geom_bar(position = "stack", stat="identity")
library(plotly)
ggplotly(plot3)

#Idem à l'horizontal
plot4 = ggplot(data2, aes(fill=Type_Variants, y=Nom_VCF, x=Nombre_Variants))+geom_bar(position = "stack", stat="identity")
ggplotly(plot4)

#Nombre de variants par VCF et par nom de gènes
nbGene = vcf_bind %>% group_by (vcf_bind$Gene, vcf_bind$VCF) %>% count(vcf_bind$Type)
data3=data_frame(nbGene)
names(data3)[1]="Nom_Genes"
names(data3)[2]="Nom_VCF"
names(data3)[3]="Type_Variants"
names(data3)[4]="Nombre_Variants"
plot5 = ggplot(data3, aes(fill=Nom_Genes, y=Nombre_Variants, x=Nom_VCF))+geom_bar(position = "stack", stat="identity")
ggplotly(plot5)

data3$NomVcf_NomGene = paste(data3$Nom_VCF,data3$Nom_Genes, sep="_")
plot6 = ggplot(data3, aes(fill=Type_Variants, y=Nombre_Variants, x=NomVcf_NomGene))+geom_bar(position = "fill", stat="identity")
plot6 = ggplot(data3, aes(fill=Type_Variants, y=Nombre_Variants, x=NomVcf_NomGene))+
  geom_bar(position = position_fill(reverse=TRUE), stat="identity")+
  labs(title = "Représentation graphique du nombre de variants en fonction des gènes et des VCF par types de variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggplotly(plot6)

cat_split %>% filter(Catégorie == "cat_ref")
cat_split %>% filter(Catégorie == "cat_12x")

#compter la taille des insertions et des délétions
#dataframe avec que délétions, insertions ou snp
vcf_deletion = vcf_bind %>% filter(Type == "délétion")
vcf_insertion = vcf_bind %>% filter(Type=="insertion")
vcf_snp = vcf_bind %>% filter(Type=="SNP")

vcf_deletion$Taille_delins = (nchar(vcf_deletion$ref)-1)
vcf_insertion$Taille_delins = (nchar(vcf_insertion$alt)-1)

cat_split %>% filter(position=="66440385")

#09/02/2024

# Représentation en boxplot des délétions
box_plotVCF <- ggplot(vcf_deletion, aes(x = VCF, y = Taille_delins)) +
  geom_boxplot(outlier.colour = "#756bb1", outlier.shape = 2, outlier.size = 3) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2ca25f") +
  coord_flip() +
  theme_classic() +
  labs(x = "VCF", y = "Taille des délétions", title = "Représentation en boxplot de la taille des délétions par VCF") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

box_plotVCF

box_plotGene <- ggplot(vcf_deletion, aes(x = Gene, y = Taille_delins)) +
  geom_boxplot(outlier.colour = "#756bb1", outlier.shape = 2, outlier.size = 3) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2ca25f") +
  coord_flip() +
  theme_classic() +
  labs(x = "Gènes", y = "Taille des délétions", title = "Représentation en boxplot de la taille des délétions par gènes") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

box_plotGene

#Représentation en boxplot des insertions
box_plotinsVCF <- ggplot(vcf_insertion, aes(x = VCF, y = Taille_delins)) +
  geom_boxplot(outlier.colour = "#d95f0e", outlier.shape = 2, outlier.size = 3) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2b8cbe") +
  coord_flip() +
  theme_classic() +
  labs(x = "VCF", y = "Taille des insertions", title = "Représentation en boxplot de la taille des insertions par VCF") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

box_plotinsVCF

box_plotinsGene <- ggplot(vcf_insertion, aes(x = Gene, y = Taille_delins)) +
  geom_boxplot(outlier.colour = "#d95f0e", outlier.shape = 2, outlier.size = 3) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2b8cbe") +
  coord_flip() +
  theme_classic() +
  labs(x = "Gènes", y = "Taille des délétions", title = "Représentation en boxplot de la taille des insertions par gènes") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

box_plotinsGene

#Renommer les colonnes pour faire un rbind des VCF délétions et insertions en vue de faire un boxplot avec les 2 représentations
colnames(vcf_deletion)[8] = c("Taille_delins")
colnames(vcf_insertion)[8] = c("Taille_delins")
vcfdelins = rbind(vcf_insertion, vcf_deletion)

#Représentation en boxplot des délétions et insertions par VCF
boxplotVCF_delins = ggplot(vcfdelins, aes(x = Taille_delins, y = VCF)) + geom_boxplot(aes(fill=Type)) + theme_classic() + 
  labs(x = "Taille des delins", y = "VCF", title = "Représentation en boxplot de la taille des insertions et délétions par VCF") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

#Représentation en boxplot des délétions et insertions par gènes
boxplotGene_delins = ggplot(vcfdelins, aes(x = Taille_delins, y = Gene)) + geom_boxplot(aes(fill=Type)) + theme_classic() + 
  labs(x = "Gènes", y = "Taille des delins", title = "Représentation en boxplot de la taille des insertions et délétions par Gènes") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))
  
  
#Représentation des délétions par gène et par VCF
  box_plotdelGene_VCF <- ggplot(vcf_deletion, aes(x = VCF, y = Taille_delins)) +
  geom_boxplot(aes(fill=Gene)) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2ca25f")+
  theme_classic() +
  labs(x = "VCF", y = "Taille des délétions", title = "Représentation en boxplot de la taille des délétions par gènes et VCF") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

box_plotdelGene_VCF

box_plotdelVCF_Gene <- ggplot(vcf_deletion, aes(x = Gene, y = Taille_delins)) +
  geom_boxplot(aes(fill=VCF)) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2ca25f")+
  theme_classic() +
  labs(x = "Gènes", y = "Taille des délétions", title = "Représentation en boxplot de la taille des délétions par gènes et VCF") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

box_plotdelVCF_Gene

box_plotinsGene_VCF <- ggplot(vcf_insertion, aes(x = Gene, y = Taille_delins)) +
  geom_boxplot(outlier.colour = "#d95f0e", outlier.shape = 2, outlier.size = 3) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2b8cbe") +
  coord_flip() +
  theme_classic() +
  labs(x = "Gènes", y = "Taille des insertions", title = "Représentation en boxplot de la taille des insertions par gènes") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))

box_plotinsGene_VCF

#Représentation des insertions par gène et par VCF
box_plotinsGene_VCF <- ggplot(vcf_insertion, aes(x = VCF, y = Taille_delins)) +
  geom_boxplot(aes(fill=Gene)) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2b8cbe") +
  coord_flip() +
  theme_classic() +
  labs(x = "VCF", y = "Taille des insertions", title = "Représentation en boxplot de la taille des insertions par gènes") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))
ggplotly(box_plotinsGene_VCF)

box_plotinsVCF_Gene <- ggplot(vcf_insertion, aes(x = Gene, y = Taille_delins)) +
  geom_boxplot(aes(fill=VCF)) +
  geom_point(stat = "summary", fun = "mean", size = 3, color = "#2b8cbe") +
  coord_flip() +
  theme_classic() +
  labs(x = "Gènes", y = "Taille des insertions", title = "Représentation en boxplot de la taille des insertions par gènes") + 
  theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"), plot.subtitle = element_text(hjust = 0.5))
box_plotinsVCF_Gene

#Voir si présence d'inversions
vcfdelins %>% filter(nchar(ref) == nchar(alt))

#Transitions/transversions - ajout colonne dans le df vcf_snp
vcf_snp$Type_substitutions = ifelse (vcf_snp$ref=="A" & vcf_snp$alt=="G", "transition_A>G",
                                     ifelse (vcf_snp$ref=="G" & vcf_snp$alt=="A", "transition_G>A",
                                            ifelse(vcf_snp$ref=="C" & vcf_snp$alt=="T", "transition_C>T",
                                                   ifelse(vcf_snp$ref=="T" & vcf_snp$alt=="C", "transition_T>C",
                                                          ifelse(vcf_snp$ref=="A" & vcf_snp$alt=="C", "transversion_A>C",
                                                                 ifelse(vcf_snp$ref=="C" & vcf_snp$alt=="A", "transversion_C>A",
                                                                        ifelse(vcf_snp$ref=="G" & vcf_snp$alt=="T", "transversion_G>T",
                                                                               ifelse(vcf_snp$ref=="T" & vcf_snp$alt=="G", "transversion_T>G",
                                                                                      ifelse(vcf_snp$ref=="A" & vcf_snp$alt=="T", "transversion_A>T",
                                                                                             ifelse(vcf_snp$ref=="T" & vcf_snp$alt=="A", "transversion_T>A",
                                                                                                    ifelse(vcf_snp$ref=="G" & vcf_snp$alt=="C", "transversion_G>C",
                                                                                                           ifelse(vcf_snp$ref=="C" & vcf_snp$alt=="G", "transversion_C>G","autre"))))))))))))


#Représentation des transitions/transversions
nb_substitutions=vcf_snp %>% group_by (Gene, VCF, Type_substitutions) %>% summarise(Nombre_substitutions=n())

data_subs=data_frame(nb_substitutions)
names(data_subs)[1]="Nom_Gene"
names(data_subs)[2]="VCF"
names(data_subs)[3]="Type_Substitutions"
names(data_subs)[4]="Nombre_Substitutions"

# Nombre de types de substitutions par type et par gène
plot_subs <- ggplot(data_subs, aes(fill = Nom_Gene, y = Nombre_Substitutions, x = Type_Substitutions)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Nombre de substitutions par type et par gène", x = "Type de substitutions", y = "Nombre de substitutions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_subs

# Nombre de types de substitutions_ref>alt par type et par VCF
plot_subs2 <- ggplot(data_subs, aes(fill = VCF, y = Nombre_Substitutions, x = Type_Substitutions)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Nombre de substitutions par type et par VCF", x = "Type de substitutions", y = "Nombre de substitutions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_subs2

#Création d'un nouveau DF avec 2 colonnes supplémentaires : Type de substitutions et ref>alt séparé
data_subs2 = data_subs
data_subs2 = data_subs2 %>% separate(Type_Substitutions, c("Type", "ref>alt"), "_")

#barplot nombre de subtitutions et types par VCF
plot_subs3 <- ggplot(data_subs2, aes(fill = Type, y = Nombre_Substitutions, x = VCF)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Nombre de substitutions par type et par VCF", x = "Type de substitutions", y = "Nombre de substitutions") +
  theme_minimal()
plot_subs3

#barplot nombre de subtitutions et types par gènes
plot_subs4 <- ggplot(data_subs2, aes(fill = Type, y = Nombre_Substitutions, x = Nom_Gene)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Nombre de substitutions par type et par gènes", x = "Type de substitutions", y = "Nombre de substitutions") +
  theme_minimal()
plot_subs4



