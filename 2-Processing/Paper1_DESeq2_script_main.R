###Differential gene expression analysis of ortholgous genes
library(DESeq2)
library(ggplot2)
library(cowplot) #add on to ggplot for better themes and figure customization
library(lemon) #to work with legends and axes in ggplot2
library(dplyr)
library(gdata)
library(RColorBrewer)
library(colorBlindness)
library(colorspace)
library(tidyverse)

theme_set(theme_classic())
setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

###Data prep for Wt x WT
orthologs = read.table("new_1to1_orthologgenelist.txt", sep="\t", head=T, comment.char="#")
head(orthologs)
cre_ortho = orthologs[, 1]
clat_ortho = orthologs[ ,2]
cumulative_ortho = orthologs[,3]

head(cre_ortho)
nrow(orthologs)
head(clat_ortho)

############
##DATA PREP#
############

#C.remanei
cre_counts = read.table("C.remanei_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
cre_counts_ortho = cre_counts[cre_ortho, -c(13:15)] ##not including Whole animal samples
rownames(cre_counts_ortho) = cumulative_ortho #changing row names
nrow(cre_counts_ortho)
colnames(cre_counts_ortho)

#C.latens
clat_counts = read.table("C.latens_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
clat_counts_ortho = clat_counts[clat_ortho, -c(13:15)]
rownames(clat_counts_ortho) = cumulative_ortho
nrow(clat_counts_ortho)
names(clat_counts_ortho)
nrow(clat_counts)
#Combining all data
wt_counts_orthologs = cbind(cre_counts_ortho, clat_counts_ortho)
colnames(wt_counts_orthologs)
colnames(wt_counts_orthologs) = gsub("X", "", colnames(wt_counts_orthologs)) ##searching and replacing the column names 

#Coldata prep
sample_name = colnames(wt_counts_orthologs)
tissue = substr(sample_name, 3,3) 
sex = substr(sample_name, 2,2)
species =  c(rep("Cre", 12), rep("Clat", 12))
batch = c(rep(c(1, 2 , 3), 8))
coldata = data.frame(sample_name, species, sex, tissue, batch)
coldata$species = factor(coldata$species)
coldata$tissue = factor(coldata$tissue)
coldata$sex = factor(coldata$sex)
coldata$batch = factor(coldata$batch)
coldata ##WM not included

#Run deseq2
dds_wt = DESeqDataSetFromMatrix(countData = wt_counts_orthologs, colData = coldata, design = ~ species +
                                  tissue + sex + species:tissue + sex:tissue + species:sex)
dds_wt = DESeq(dds_wt)
resultsNames(dds_wt)

# get the model matrix
mod_mat_wt <- model.matrix(design(dds_wt), colData(dds_wt))

#Effect of Species on gene expression 
#Define coefficient vectors for each condition
cre = colMeans(mod_mat_wt[dds_wt$species == "Cre", ])
clat = colMeans(mod_mat_wt[dds_wt$species == "Clat", ])


crexclat = results(dds_wt, contrast = cre - clat, alpha = 0.05) #baseline clat
crexclat_res <- cbind(as.data.frame(crexclat), as.data.frame(crexclat) %>%
                        mutate(crexclat = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                        dplyr::select(crexclat))
crexclat_res = na.omit(crexclat_res)
crexclat_res1 = crexclat_res %>% filter(crexclat_res$crexclat != 0)
crexclat_res2 = crexclat_res %>% filter(crexclat_res$crexclat == 0)

summary(crexclat)##upregulated in Cre, reference level is Clatens
nrow(crexclat_res)
table(crexclat_res$crexclat) ###Number of species-biased and conserved genes
mod_mat_wt(crexclat)

#############################################################################################################################
#########CONVERTING GENE NAMES FOR GPROFILER#################################################################################
#############################################################################################################################
####Saving species-biased genes names to view in GPROFILER 

###Loading cremanie and clatens gene names for older genomes (from Daniel)
remanei_gene_names = read.table("DESeq2/remanei_gene_names.txt", fill = T)
latens_gene_names = read.table("DESeq2/latens_gene_names.txt", fill = T)
head(latens_gene_names)

###getting species-biased and unbiased genes from DESeq2 (above) analysis
crexclat_up = crexclat_res %>% filter(crexclat_res$crexclat == 1) #Higher in C. remanei
crexclat_down = crexclat_res %>% filter(crexclat_res$crexclat == -1) #Higher in C. latens
crexclat_nobias = crexclat_res %>% filter(crexclat_res$crexclat == 0) #Unbiased
crexclat_speciesbiased = crexclat_res %>% filter(crexclat_res$crexclat != 0) #Unbiased
crexclat_speciesbiased =rownames_to_column(crexclat_speciesbiased)
crexclat_speciesbiased = separate(crexclat_speciesbiased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
crexclat_speciesbiased = separate(crexclat_speciesbiased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
head(crexclat_speciesbiased)

##Saving species-baised and unbiased genes
# write.table(crexclat_up, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/crexclat_up.txt",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)
# 
# write.table(crexclat_down, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/crexclat_down.txt",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)
# 
# write.table(crexclat_speciesbiased, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/species_biased.txt",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)

###Getting gene names for genes with higher expressison in C. remanei
cre_biased = rownames_to_column(crexclat_up)
cre_biased = separate(cre_biased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
cre_biased = separate(cre_biased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2] #run it right after rprevious step
cre_biased_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% cre_biased$Crenames)[,2] #Cre old gene names
cre_biased_clat = subset(latens_gene_names, latens_gene_names$V4 %in% cre_biased$Clatnames)[,2] #corresponding gene name sin C. latens


# write.table(cre_biased_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/cre_biased_cregenes.txt", 
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.table(cre_biased_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/cre_biased_clatgenes.txt", 
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

###Getting gene names for genes with higher expresison in C. latens
clat_biased = rownames_to_column(crexclat_down)
clat_biased = separate(clat_biased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
clat_biased = separate(clat_biased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2] #run it right after rprevious step
clat_biased_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% clat_biased$Crenames)[,2] #Cre old gene names
clat_biased_clat = subset(latens_gene_names, latens_gene_names$V4 %in% clat_biased$Clatnames)[,2] #corresponding gene name sin C. latens

nrow(clat_biased)
length(clat_biased_clat)
# write.table(clat_biased_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/clat_biased_cregenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.table(clat_biased_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/clat_biased_clatgenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)


###Getting gene names for genes with no diffrence in expression across species
sp_unbiased = rownames_to_column(crexclat_nobias)
sp_unbiased = separate(sp_unbiased, rowname, sep = 23, into = c("Crenames", "Clatnames"))
sp_unbiased = separate(sp_unbiased, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2] #run it right after rprevious step
sp_unbiased_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% sp_unbiased$Crenames)[,2] #Cre old gene names
sp_unbiased_clat = subset(latens_gene_names, latens_gene_names$V4 %in% sp_unbiased$Clatnames)[,2] #corresponding gene name sin C. latens
length(sp_unbiased_cre)
nrow(sp_unbiased)
length(sp_unbiased_clat)
head(sp_unbiased)
# write.table(sp_unbiased_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/sp_unbiased_cregenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.table(sp_unbiased_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/sp_unbiased_clatgenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(sp_unbiased, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/sp_unbiased.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

####################################################################################################################################
####PART I##########################################################################################################################
####################################################################################################################################


####FIGURE 1A-VISUALIZING SPECIES-BIASED GENES #####################################################################################
####################################################################################################################################

table(crexclat_res$crexclat)

sp_genes_count = data.frame(gene_category = c("Conserved", "Conserved", "DEG","DEG"),
                            genes_count = c(6865, 0, 3339, 3224), #creating a proxy for conserved deg cetagory
                            DEG_category = c("Clat-b", "Cre-b","Clat-b","Cre-b"))

ggplot(sp_genes_count, aes(x = sp_genes_count$gene_category, y = sp_genes_count$genes_count, fill = sp_genes_count$DEG_category)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,7500,500))+
  ggtitle("Species-biased genes")


####################################################################################################################################
####FIGURE 1B-VISUALIZING SPECIES-BIASED GENES ON CHROMOSOMES#######################################################################
####################################################################################################################################

##adding chromosomal information to the dataframe
chrom_info = read.table("Cremanie_genenames_and_chromosome.txt", fill = T)
head(chrom_info)
chrom_info_X = subset(chrom_info, V1 =="X")
chrom_info_I = subset(chrom_info, V1 =="I")
chrom_info_II = subset(chrom_info, V1 =="II")
chrom_info_III = subset(chrom_info, V1 =="III")
chrom_info_IV = subset(chrom_info, V1 =="IV")
chrom_info_V = subset(chrom_info, V1 =="V")

###get X-linked orthologous genes
ortho_X = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_X$V2)
ortho_X$chromosome = "X"
ortho_I = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_I$V2)
ortho_I$chromosome = "I"
ortho_II = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_II$V2)
ortho_II$chromosome = "II"
ortho_III = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_III$V2)
ortho_III$chromosome = "III"
ortho_IV = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_IV$V2)
ortho_IV$chromosome = "IV"
ortho_V = orthologs %>% filter(orthologs$C..remanei.Gene.name %in% chrom_info_V$V2)
ortho_V$chromosome = "V"
nrow(ortho_X)

orthologs_chr = rbind(ortho_I, ortho_II, ortho_III, ortho_IV, ortho_V, ortho_X)
nrow(orthologs_chr) ##3loss of genes on scaffolds
nrow(orthologs)
nrow(ortho_X)

# write.table(orthologs_chr, file = "orthologs_chr.txt", sep = "\t",quote = FALSE,
#            row.names = FALSE)

##changing rownames to new column
head(crexclat_res)
crexclat_res_a = na.omit(crexclat_res[orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name, ])
crexclat_res_a = rownames_to_column(crexclat_res_a)
nrow(crexclat_res_a)
head(orthologs_chr)
rownames(orthologs_chr) = orthologs_chr$C..remanei.Gene.name_C.latens.Gene.name  
crexclat_res_chr = (na.omit(orthologs_chr[crexclat_res_a$rowname, ]) %>% cbind(crexclat_res_a))
View(crexclat_res_chr)

spgenes_freq = as.data.frame(table(crexclat_res_chr$chromosome, crexclat_res_chr$crexclat))
View(spgenes_freq)
 

####################################################################################################
#####Plotting species-biased genes across chromosomes###############################################
####################################################################################################
ggplot(spgenes_freq, aes(x = spgenes_freq$Var1, y = spgenes_freq$Freq, fill = spgenes_freq$Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  #coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,6500,500))+
  ggtitle("Species-biased genes")

####################################################################################################################################
#####Proportion of species-biased genes on chromosomes###################################################3

##plotting proportions across chromosomes
library(ggstats)

ggplot(spgenes_freq) +
  aes(x = spgenes_freq$Var1, fill = spgenes_freq$Var2, weight = spgenes_freq$Freq, by = spgenes_freq$Var1) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1,0.25))


#########CONVERTING GENE NAMES FOR GPROFILER#################################################################################

####Getting genes in chromsomes 1 and looking at their old assembly gene names for gprofiler
ch1_spgenes = subset(spgenes_chr, spgenes_chr$V1 == "I" & spgenes_chr$crexclat != 0)[,3]
ch1_spgenes = subset(remanei_gene_names, remanei_gene_names$V4 %in% ch1_spgenes)[,2]

chX_spgenes = subset(spgenes_chr, spgenes_chr$V1 == "X" & spgenes_chr$crexclat != 0)[,3]
chX_spgenes = subset(remanei_gene_names, remanei_gene_names$V4 %in% chX_spgenes)[,2]

# write.table(ch1_spgenes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/ch1_spgenes_crenames.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(chX_spgenes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/chX_spgenes_crenames.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################################################################
#######REDRAWING FIGURES FROM GPROFILER FOR SPECIES-BIASED GENES##############################################
##############################################################################################################

#####bidirectional barplot for enrichment from gprofiler
library(ggplot2)

grprofiler_1 = read.table("../New annotation/Results/Plots/Part 1/gprofiler_species_biasedgenes.txt", sep="\t", head=T, comment.char="#")
#View(grprofiler_1)

###This saved txt file was obtained after running grpofiler on species-biased genes with c.remanei gene (old) gene names

grprofiler_1 = grprofiler_1 %>% filter(grprofiler_1$Type != "Underrepresentation")

library(forcats)#to easily reorder varibales 
ggplot(grprofiler_1,aes(x=fct_reorder(grprofiler_1$term_name, grprofiler_1$negative_log10_of_adjusted_p_value), y = intersection_size, fill=negative_log10_of_adjusted_p_value))+ 
  geom_bar(stat="identity", position="identity")+
  scale_fill_gradient(low="#cc9966", high="#663300")+
  scale_y_continuous(limits = c(0, max(grprofiler_1$intersection_size))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+
  coord_flip()



#####################################################################################################
####PART II##########################################################################################
#####################################################################################################
#######VISUALIZING SPECIES BIASED GENES USING WGCNA MODULES
###First run WGCNA_WT.R script or directly use saved network (As below) 
###Data after running WGCNA_WT.R

###loading file which contains genes and their corresponding modules (obtained after running WGCNA_WT.R)
bwnet_df = read.table("Results/Plots/Part 2/WGCNA_modulegenes.txt", fill = T)
#View(bwnet_df)

table(crexclat_res$crexclat) ##table of species biased genes, total 13428 genes including unbiased

sp_genes_modules = subset(bwnet_df, bwnet_df$V1 %in% crexclat_res3$rowname)
sp_genes_modules = subset(crexclat_res3, crexclat_res3$rowname %in% sp_genes_modules$V1) %>% cbind(sp_genes_modules$V2)
sp_genes_modules_freq = as.data.frame(table(sp_genes_modules$crexclat, sp_genes_modules$`sp_genes_modules$V2`))
nrow(sp_genes_modules)
head(sp_genes_modules_freq)
table(sp_genes_modules$crexclat, sp_genes_modules$`sp_genes_modules$V2`)
#View(sp_genes_modules_freq)

###Adding module informaiton 
sp_genes_modules_freq$Module = "M0"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="yellow"] = "M1"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="purple"] = "M2"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="greenyellow"] = "M3"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="blue"] = "M4"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="turquoise"] = "M5"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="red"] = "M6"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="green"] = "M7"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="pink"] = "M8"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="black"] = "M9"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="tan"] = "M10"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="magenta"] = "M11"
sp_genes_modules_freq$Module[sp_genes_modules_freq$Var2=="brown"] = "M12"

######LOOKING AT MODULES UNDERLYING SPECIES-BIASED AND UNBIASED GENES
ggplot(sp_genes_modules_freq, aes(x = sp_genes_modules_freq$Var1, y = sp_genes_modules_freq$Freq, fill = sp_genes_modules_freq$Module)) +
  geom_bar(stat = "identity") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,6500,500))+
  ggtitle("Species-biased genes")

##Visualizing only species-biased genes
ggplot(sp_genes_modules_freq, aes(x = sp_genes_modules_freq$Var1, y = sp_genes_modules_freq$Freq)) +
  geom_bar(stat = "identity") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  coord_cartesian(ylim = c(0, 6500))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 26,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,6500,500))+
  ggtitle("Species-biased genes")


########################################################################################################################################################
#######PART 4
####FINDING INTERACTIONS AND EFFFECTS OF TISSUE, SEX AND SPECIES ON GENE EXPRESSION
###################################################################################################################################################################33

#######VISUALIZING TISSUE-BIASED AND SEX-BIASED GENES UNDERLYING SPECIES-BIASED GENES####

#FINDING GENES THAT ARE TISSUE-BIASED (INCLUDING THOSE THAT ARE DUE TO INTERATION TERM) 

#Define coefficient vectors for each condition
gonad = colMeans(mod_mat_wt[dds_wt$tissue == "G", ])
soma = colMeans(mod_mat_wt[dds_wt$tissue == "S", ])

gxs = results(dds_wt, contrast = soma - gonad, alpha = 0.05) #baseline clat
gxs_res <- cbind(as.data.frame(gxs), as.data.frame(gxs) %>%
                   mutate(gxs = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                   dplyr::select(gxs))
nrow(gxs_res)
summary(gxs)##upregulated in Soma, reference level is gonad
gxs_res = na.omit(gxs_res)
gxs_res1 = gxs_res %>% filter(gxs_res$gxs != 0)
gxs_res2 = gxs_res %>% filter(gxs_res$gxs == 0)

table(gxs_res$gxs)
design(dds_wt)
#write.csv(gxs_res, file = "Results/gxs_res_genes.csv")

################################################################################
#FINDING GENES THAT ARE SEX-BIASED (INCLUDING THOSE THAT ARE DUE TO INTERATION TERM) 
################################################################################

#Define coefficient vectors for each condition
male = colMeans(mod_mat_wt[dds_wt$sex == "M", ])
female = colMeans(mod_mat_wt[dds_wt$sex == "F", ])

mxf = results(dds_wt, contrast = male - female, alpha = 0.05) #baseline clat
mxf_res <- cbind(as.data.frame(mxf), as.data.frame(mxf) %>%
                   mutate(mxf = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                   dplyr::select(mxf))
nrow(mxf_res)##upregulated in Cre, reference level is Clatens
summary(mxf)
mxf_res = na.omit(mxf_res)
mxf_res1 = mxf_res %>% filter(mxf_res$mxf != 0)
mxf_res2 = mxf_res %>% filter(mxf_res$mxf == 0)

table(mxf_res$mxf)
write.csv(mxf_res, file = "Results/mxf_res_genes.csv")

#############################
#####COMBINING DATA
#############################

all_total_biasgene_info = cbind(crexclat_res, gxs_res, mxf_res) ###these DEGs include one swith interactions
table(all_total_biasgene_info$crexclat, all_total_biasgene_info$gxs) ##first data as rownames 


#################################################################################################
#Multiple Effects of Species and Tissue differences on gene expression 
##################################################################################################

#Define coefficient vectors for each condition
# cre_s = colMeans(mod_mat_wt[dds_wt$species == "Cre" & dds_wt$tissue == "S", ])
# cre_g = colMeans(mod_mat_wt[dds_wt$species == "Cre" & dds_wt$tissue == "G", ])
# clat_s = colMeans(mod_mat_wt[dds_wt$species == "Clat" & dds_wt$tissue == "S", ])
# clat_g = colMeans(mod_mat_wt[dds_wt$species == "Clat" & dds_wt$tissue == "G", ])
#spxtsexx = results(dds_wt, contrast = (cre_s - cre_g - clat_s + clat_g), alpha = 0.05)

######################################################################
#Effect of Species and Tissue differences on gene expression 

spxts = results(dds_wt, contrast = list("speciesCre.tissueS"), alpha = 0.05) #baseline clat
spxts_res <- cbind(as.data.frame(spxts), as.data.frame(spxts) %>%
                     mutate(spxts = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                     dplyr::select(spxts))
summary(spxts)
nrow(spxts_res)
spxts_res = na.omit(spxts_res)
spxts_res1 = spxts_res %>% filter(spxts_res$spxts != 0 )
spxts_res2 = spxts_res %>% filter(spxts_res$spxts == 0 )
table(spxts_res1$spxts)

######################################################################
#Effect of Species and Sex differences on gene expression 

#Define coefficient vectors for each condition
spxsex = results(dds_wt, contrast = list("speciesCre.sexM"), alpha = 0.05) #baseline clat
spxsex_res <- cbind(as.data.frame(spxsex), as.data.frame(spxsex) %>%
                      mutate(spxsex = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                      dplyr::select(spxsex))
nrow(spxsex_res)
spxsex_res = na.omit(spxsex_res)
spxsex_res1 = spxsex_res %>% filter(spxsex_res$spxsex != 0 )
spxsex_res2 = spxsex_res %>% filter(spxsex_res$spxsex == 0 )

table(spxsex_res1$spxsex)

######################################################################
#Effect of Tissue and Sex differences on gene expression 

#Define coefficient vectors for each condition
tsxsex = results(dds_wt, contrast = list("tissueS.sexM"), alpha = 0.05) #baseline clat
tsxsex_res <- cbind(as.data.frame(tsxsex), as.data.frame(tsxsex) %>%
                      mutate(tsxsex = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
                      dplyr::select(tsxsex))
summary(tsxsex)##upregulated in Cre, reference level is Clatens
table(tsxsex_res$tsxsex)
tsxsex_res = na.omit(tsxsex_res)
tsxsex_res1 = tsxsex_res %>% filter(tsxsex_res$tsxsex != 0 )
tsxsex_res2 = tsxsex_res %>% filter(tsxsex_res$tsxsex == 0 )

############################################################################################################################################
#####################################################################################################################
#####Finding genes within each category
############################################################################################################################################
############################################################################################################################################

library(UpSetR)
genelist  = list("Sex" = row.names(mxf_res1), "SpxSex" = row.names(spxsex_res1), "TSxSex" = row.names(tsxsex_res1),
                 "Species" = row.names(crexclat_res1),"SpxTS" = row.names(spxts_res1), "Tissue" = row.names(gxs_res1))
View(genelist)
UpSetR::upset(fromList(genelist), nsets = 12, nintersects = NA, order.by = "freq", group.by = "degree",
              mainbar.y.label = "Number of common genes", sets.x.label = "Total number of genes", set_size.show = TRUE, 
              set_size.scale_max = 11000, set_size.numbers_size = 12, line.size = 0.5, point.size = 3)

UpSetR::upset(fromList(genelist), nsets = 12)

nrow(fromList(genelist))

library(ComplexHeatmap)

m3 = make_comb_mat(genelist, min_set_size = 1)
extract_comb(m3, "000001") #get gene names for a particular intersection category

#EXTRA COMMANDS TO GET INFORMATION ABOUT SETS
# m3[1:5] #to get the first 5 least comination/intersection set sizes 
# str(m3)
# set_name(m3)
# comb_name(m3)
# set_size(m3)
# comb_size(m3)
# comb_degree(m3)
# t(m3)
# m3[, 1:3]

##########################################################################################3
#####MAKING AN UPSET PLOT#############################################################3333

UpSet(m3)
UpSet(m3[comb_size(m3) >= 300])
UpSet(m3[comb_degree(m3) == 2])
# UpSet(m3, pt_size = unit(5, "mm"), lwd = 3,
#       comb_col = c("red", "blue", "black", "green", "yellow", "grey", "orange")[comb_degree(m3)])
# UpSet(m3, comb_col = "#0000FF",bg_col = "#F0F0FF", bg_pt_col = "#CCCCFF") #bg_col colours background, bg_pt_col colours unconnected points

# UpSet(m3, top_annotation = upset_top_annotation(m3, 
#                                                 gp = gpar(col = comb_degree(m3))))
# UpSet(m3, right_annotation = upset_right_annotation(m3, 
#                                                     gp = gpar(fill = "green"),
#                                                     annotation_name_side = "bottom",
#                                                     axis_param = list(side = "bottom")))
# UpSet(m3, top_annotation = HeatmapAnnotation(
#   degree = as.character(comb_degree(m3)),
#   "Intersection\nsize" = anno_barplot(comb_size(m3), 
#                                       border = FALSE, 
#                                       gp = gpar(fill = "black"), 
#                                       height = unit(2, "cm")
#   ), 
#   annotation_name_side = "left", 
#   annotation_name_rot = 0))
# UpSet(m3, left_annotation = upset_left_annotation(m3))
# m4 = m3[comb_size(m3) >= 300]
# UpSet(m4, top_annotation = upset_top_annotation(m4, add_numbers = TRUE),
#       right_annotation = upset_right_annotation(m4, add_numbers = TRUE), comb_order = comb_size(m4) )
# # UpSet(m4, comb_order = comb_size(m4))
# UpSet(set_order = order(set_size(m4)))
# cs = comb_size(m3)
# UpSet(m3, comb_order = order(comb_degree(m3), -cs))

###########################################################
##Extracting different gene lists
# set_name(m3)
# str(comb_name(m3))
# set_size(m3)
# str(comb_size(m3))
# str(comb_degree(m3))

combos = as.data.frame(comb_degree(m3), row.names = NULL)
combos = rownames_to_column(combos)

###########################################################
##Extracting genes for a given condition and their modules#
###########################################################

#getting genes and their modules, 
#bwnet_df is the network dataframe with v1 genes names and v2 modules

module.gene.mapping = bwnet_df %>% column_to_rownames(var = "V1") 

########################################
###One main effect######################

tissue_only = extract_comb(m3, "000001")
tissue_only = as.data.frame(tissue_only) 
colnames(tissue_only) = "Genes"
tissue_only$condition = "tissue_only"
tissue_only$Module = module.gene.mapping[tissue_only$Genes, ]
tissue_only = na.omit(tissue_only)
head(tissue_only)

sex_only = extract_comb(m3, "100000")
sex_only = as.data.frame(sex_only)
colnames(sex_only) = "Genes"
sex_only$condition = "sex_only"
sex_only$Module = module.gene.mapping[sex_only$Genes, ]
head(sex_only)
sex_only = na.omit(sex_only)
table(sex_only$Module)
head(sex_only)

species_only = extract_comb(m3, "000100")
species_only = as.data.frame(species_only)
colnames(species_only) = "Genes"
species_only$condition = "species_only"
species_only$Module = module.gene.mapping[species_only$Genes, ]
species_only = na.omit(species_only)
head(species_only)

########################################
###Two main effects#####################

TXS = extract_comb(m3, "100001")
TXS = as.data.frame(TXS)
colnames(TXS) = "Genes"
TXS$condition = "TXS"
TXS$Module = module.gene.mapping[TXS$Genes, ]
head(TXS)
table(TXS$Module)
TXS = na.omit(TXS)

TXSP = extract_comb(m3, "000101")
TXSP = as.data.frame(TXSP)
colnames(TXSP) = "Genes"
TXSP$condition = "TXSP"
TXSP$Module = module.gene.mapping[TXSP$Genes, ]
head(TXSP)
table(TXSP$Module)

SXSP = extract_comb(m3, "100100")
SXSP = as.data.frame(SXSP)
colnames(SXSP) = "Genes"
SXSP$condition = "SXSP"
SXSP$Module = module.gene.mapping[SXSP$Genes, ]
head(SXSP)
table(SXSP$Module)

########################################
####THREE MAIN EFFECTS##################

TXSXSP = extract_comb(m3, "100101")
TXSXSP = as.data.frame(TXSXSP)
colnames(TXSXSP) = "Genes"
TXSXSP$condition = "TXSXSP"
TXSXSP$Module = module.gene.mapping[TXSXSP$Genes, ]
head(TXSXSP)
table(TXSXSP$Module)

########################################
###Main effect and one innteraction

TXSXT_S = extract_comb(m3, "101001")
TXSXT_S = as.data.frame(TXSXT_S)
colnames(TXSXT_S) = "Genes"
TXSXT_S$condition = "TXSXT_S"
TXSXT_S$Module = module.gene.mapping[TXSXT_S$Genes, ]
head(TXSXT_S)
table(TXSXT_S$Module)
nrow(TXSXT_S)

TXSPXSXS_SP = extract_comb(m3, "110101")
TXSPXSXS_SP = as.data.frame(TXSPXSXS_SP)
colnames(TXSPXSXS_SP) = "Genes"
TXSPXSXS_SP$condition = "TXSPXSXS_SP"
TXSPXSXS_SP$Module = module.gene.mapping[TXSPXSXS_SP$Genes, ]
head(TXSPXSXS_SP)
table(TXSPXSXS_SP$Module)
nrow(TXSPXSXS_SP)

TXSPXSXS_T = extract_comb(m3, "101101")
TXSPXSXS_T = as.data.frame(TXSPXSXS_T)
colnames(TXSPXSXS_T) = "Genes"
TXSPXSXS_T$condition = "TXSPXSXS_T"
TXSPXSXS_T$Module = module.gene.mapping[TXSPXSXS_T$Genes, ]
head(TXSPXSXS_T)
table(TXSPXSXS_T$Module)
nrow(TXSPXSXS_T)

########################################
########COMBINING DATA##################
########################################

cumulative_data = bind_rows(tissue_only, sex_only, species_only, TXS,
                            TXSP, SXSP, TXSXSP, TXSXT_S, TXSPXSXS_SP, TXSPXSXS_T)
cumulative_data = na.omit(cumulative_data)

cumulative_datatable = as.data.frame(unclass(table(cumulative_data$condition, cumulative_data$Module)))
cumulative_datatable$sums = rowSums(cumulative_datatable )
cumulative_datatable[nrow(cumulative_datatable) + 1,] = colSums(cumulative_datatable )

cumulative_data$Module = factor(cumulative_data$Module, levels = c("turquoise", "blue", "brown", "yellow", "green", 
                                                                   "pink", "red", "black","magenta", "greenyellow", "purple",
                                                                   "tan", "grey"))
cumulative_data$condition = factor(cumulative_data$condition, levels = c("TXS", "TXSXSP", "tissue_only", "TXSP", "species_only", 
                                                                         "TXSPXSXS_SP", "TXSXT_S", "SXSP", "sex_only", "TXSPXSXS_T"))

fig1_b_data = as.data.frame(unclass(table(cumulative_data$condition, cumulative_data$Module)))
colnames(fig1_b_data) = c("M5", "M4", "M12", "M1", "M7", "M8", "M6", "M9", "M11", "M3", "M2", "M10", "M0")

###########################################################################################
##Upset plot for first 10 categories#######################################################
###########################################################################################

m4 = m3[comb_size(m3) >= 300]
UpSet(m4, top_annotation = upset_top_annotation(m4, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m4, add_numbers = TRUE) )

#######################################################
#Module composition for first 10 categories
#######################################################

ggplot(cumulative_data, aes(x=condition, fill = Module)) +
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("#7beadf", "#87CEFA", "#d06a49","#F0E68C", "#556B2F",
                               "#FFC0CB", "#ff3333", "#000000", "#BA55D3", "#73cf17",
                               "#663399", "#D2B48C","#D3D3D3"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,3500), breaks=seq(0,3500,500))

# ggsave(filename = "DESeq2/modules of intersection.pdf", dpi = 300, width = 9, height = 6 )

#######################################################
##Gene compostion of each module
#######################################################

cumulative_data %>% 
  ggplot(aes(Module, fill = condition)) +
  geom_bar(position = "stack") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,4000), breaks=seq(0,4000,500))


#ggsave(filename = "DESeq2/gene composition of modules.pdf", dpi = 300, width = 12, height = 6 )

#######################################################
###Number of genes in each module
#######################################################
genes_in_modules = as.data.frame((bwnet$colors))

genes_in_modules %>% ggplot(aes(genes_in_modules$`(bwnet$colors)`)) + geom_bar(position = "stack") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,4500), breaks=seq(0,4500,500))


df1=data.frame(Group = c(rep("T", 2), rep("S", 2), rep("Sp", 2), rep("TxS", 2), rep("SpxS", 2), rep("TxSp", 2)),
               No.ofgenes = c(1050, 8909, 1418, 7899, 1266, 5297, 1056, 876, 1241, 544, 450, 0), Category = rep(c("Total", "Top10"), 6))
#data_cord= data.frame(x = c(2,3), y = c(60, 60))

df1$Category = factor(df1$Category, levels = c("Total", "Top10"))
ggplot(df1, aes(x=Group, y=No.ofgenes, fill = Category))+
  geom_bar(stat = "identity", show.legend=T, alpha=1)+
  xlab(" ")+ 
  ylab("Number of genes")+
  labs(title = "Overall distribution of genes")+
  theme(plot.title = element_text(size = 23, face = "bold", hjust = 0.5, vjust = 0.5))+
  coord_cartesian(ylim = c(0, 10000))+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = c("#c08cc0", "#d5d0e2"))+
  scale_x_discrete(limits = c("T","S", "Sp", "TxS", "SpxS", "TxSp")) ##reordering character x-axis
#geom_hline(yintercept = 5000, linetype = "dashed")

# write.csv(cumulative_data, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Cumulative_DGE_genelist.csv", 
#             row.names = FALSE, col.names = TRUE, quote = FALSE)


