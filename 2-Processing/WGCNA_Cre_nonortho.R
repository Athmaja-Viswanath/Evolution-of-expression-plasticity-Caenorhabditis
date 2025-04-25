###########################################################################################
######################PART IIIA############################################################
###########################################################################################

##WGCNA for C. remanei non-orhtologous genes ONLY

library(WGCNA)
library(DESeq2)
#library(GEOquery) for data from NCBI
library(tidyverse)
library(ggplot2)
library(CorLevelPlot)
library(gridExtra)


setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

###Data prep for Wt x WT

#C.remanei
cre_allcounts = read.table("C.remanei_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
cre_allcounts = rownames_to_column(cre_allcounts)
"%notin%" = Negate("%in%")
cre_allcounts = cre_allcounts[cre_allcounts$rowname %notin% orthologs$C..remanei.Gene.name, ]
nrow(cre_allcounts)
cre_allcounts = remove_rownames(cre_allcounts)
cre_allcounts = column_to_rownames(cre_allcounts, var = "rowname")
cre_allcounts = cre_allcounts[ ,-c(13:15)]
names(cre_allcounts)
colnames(cre_allcounts) = gsub("X", "", colnames(cre_allcounts)) ##searching and replacing the column names 
colnames(cre_allcounts)

#COLDATA PREP
sample_name_cre = colnames(cre_allcounts)
tissue_cre = substr(sample_name_cre, 3,3) 
sex_cre = substr(sample_name_cre, 2,2)
batch_cre = c(rep(c(1, 2 , 3), 4))
coldata_cre = data.frame(sample_name_cre, sex_cre, tissue_cre, batch_cre)
View(coldata_cre)


head(cre_allcounts)

###QUALITY CONTROL TO DETECT OUTLIERS

gsg_cre = goodSamplesGenes(t(cre_allcounts))
summary(gsg_cre)

gsg_cre$allOK ##gives outliers, if TRUE = all genes and samples are good and have no outliers

table(gsg_cre$goodGenes) #number of genes that are outliers = FALSE
table(gsg_cre$goodSamples)


cre_allcounts = cre_allcounts[gsg_cre$goodGenes == TRUE, ] ##filtering good genes
nrow(cre_allcounts)


####heirarchical clustering to detect outliers

htree_cre = hclust(dist(t(cre_allcounts)), method = "average")
plot(htree_cre) 


##PCA for detecting outliers

pca_cre = prcomp(t(cre_allcounts))
pca_data_cre = pca_cre$x ####PCA calculations for all the samples

#Calculating the variance explained by each principal component

pca.var_cre = pca_cre$sdev^2
pca.var.percent_cre = round(pca.var_cre/sum(pca.var_cre)*100, digits = 2)

pca_data_cre = as.data.frame(pca_data_cre)

ggplot(pca_data_cre, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data_cre)) +
  labs( x = paste0("PC1:", pca.var.percent_cre[1], "%"), 
        y = paste0("PC2:", pca.var.percent_cre[2], "%"))

##Correcting for batch effects
# library(sva)
# 
# counts_orthologs_adjusted <- ComBat_seq(counts_orthologs, batch=batch, group=NULL)
# counts_orthologs = counts_orthologs_adjusted
# 
# 
# count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
# batch <- c(rep(1, 4), rep(2, 4))
# 
# adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)


##Normalization
###Create a deseq2 dataset

#create coldata if needed

coldata_cre2 = coldata_cre[, -1]     ##changing column to rownames
rownames(coldata_cre2) = coldata_cre[ , 1]

all(rownames(coldata_cre2) %in% colnames(cre_allcounts))
all(rownames(coldata_cre2) == colnames(cre_allcounts))

dds_cre = DESeqDataSetFromMatrix(countData = cre_allcounts, 
                             colData = coldata_cre2, 
                             design = ~1) ##not specifiying a model
View(coldata_cre2)

####removing low count genes

dds75_cre = dds_cre[rowSums(counts(dds_cre) >= 15) >= 6, ]

nrow(dds75_cre)


##variance stabilize transformaiton 

dds_norm_cre = vst(dds75_cre)


#get normalized counts
#transfromed data fro downstream analysis

norma.counts_cre = assay(dds_norm_cre) %>% 
  t()

#View(norma.counts)



#####4. Network construction

##Choose a set of soft threshold powers

power = c(c(1:10), seq(from = 12, to = 50, by = 2))

sft_cre = pickSoftThreshold(norma.counts_cre,
                        powerVector = power,
                        networkType = "signed",
                        verbose = 5)

sft.data_cre = sft_cre$fitIndices ##we will use the R suqared value and mean connectivity , max R^2 and min mean connectivity

###Visualization to pick the right power
a1 = ggplot(sft.data_cre, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()


a2 = ggplot(sft.data_cre, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  #geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)  ##we need higher R^2 and low mean connectivity


###converting matric to numeric

norma.counts_cre[] = sapply(norma.counts_cre, as.numeric)  ##different step than official tutorial

soft_power_cre = 26  #3old was 24

temp_cor = cor #to prevent WGCNa from using other cor function

cor = WGCNA::cor


###Memory estimate wrt blocksize

##can take 20-30 mins to run

bwnet_cre = blockwiseModules(norma.counts_cre, 
                         maxBlockSize = 14000, ##depends on the ram of the system 4gb = 8-10k, 16gb = 20,000, 232gb = 30,000 
                         TOMType = "signed",
                         networkType = "signed",
                         power = soft_power_cre,
                         mergeCutHeight = 0.25,#threshold that we want to merge similar modules at
                         numericLabels = FALSE, #want the module names to be colours if not, then say TRUE
                         randomSeed = 1234,
                         verbose = 3) 

cor = temp_cor

#5. Module eigenvalue

module_eigengenes_cre = bwnet_cre$MEs
View(module_eigengenes_cre)



#write.csv(module_eigengenes_cre, file = "WGCNA//For WT_without WM/Non-orthologous genes_Cre/C.remanie module eigengenes values.csv")

###CALCULATE AVERAGE ACROSS REPLICATES FOR EACH MODULE

##get number of genes in each module
table(bwnet_cre$colors)

names(bwnet_cre$colors)


##############################################################################################################
####Visualize the numebr of genes in each module##################################################################
bwnet_df_cre = as.data.frame(bwnet_cre$colors)
modules_counts_cre = as.data.frame(table(bwnet_df_cre$`bwnet_cre$colors`))

ggplot(modules_counts_cre, aes(x=Var1, y = Freq, fill = Var1))+
  geom_bar(stat = "identity") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  coord_cartesian(ylim = c(0, 1300))+
  #scale_fill_manual(values = c( "#afced0", "#4b7d81", "#696967", '#DCDDDF'))+
  theme(axis.title.y = element_text(size = 14, hjust = 0.5))+
  theme(axis.text.x = element_text(size = 10,face = "bold", colour = "black"))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 20, hjust = 0.5))+
  theme(axis.text.y = element_text(color="black", size=15))+
  guides(fill=guide_legend(title=" "))+
  # scale_colour_brewer(palette = 2)+
  #scale_x_discrete(limits = c("Conserved","Differentially expressed", "C. nigoni dominant", "Ambiguous"))+ ##reordering character x-axis
  expand_limits(y=0)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,1300,100))+
  ggtitle("Cre-specific WGCNA Modules")

# pick out a few modules of interest here

# module_df_cre <- data.frame(
#   gene_id = names(bwnet_cre$colors),
#   colors = labels2colors(bwnet_cre$colors)) #it doesnt match the right module with the right gene

#write_delim(module_df, file = "WGCNA/gene_modules.txt", delim = "\t")

# Pull out list of genes in that module

# modules_of_interest_cre = c("brown")
# 
# genes_of_interest_cre = module_df_cre %>%
#   subset(colors %in% modules_of_interest_cre)
# 
# submod_cre = module_df_cre %>%
#   subset(colors %in% modules_of_interest_cre)
# 
# row.names(module_df_cre) = module_df_cre$gene_id
# 
# # Get normalized expression for those genes
# norma.counts.df = as.data.frame(t(norma.counts))
# brown_list = norma.counts.df[submod$gene_id, ]
# brown_list[1:5, 1:5]
# 
# TOM = TOMsimilarityFromExpr(t(brown_list),
#                             power = soft_power)
# 

# Add gene names to row and columns
# row.names(TOM) = row.names(brown_list)
# colnames(TOM) = row.names(brown_list)
# 
# 
# brown_list = data.frame(TOM) %>%
#   mutate(
#     gene1 = row.names(.)
#   ) %>%
#   pivot_longer(-gene1) %>%
#   dplyr::rename(gene2 = name, correlation = value) %>%
#   unique() %>%
#   subset(!(gene1==gene2)) %>%
#   mutate(
#     module1 = module_df[gene1,]$colors,
#     module2 = module_df[gene2,]$colors
#   )
# 
# head(brown_list)

# Export Network file to be read into Cytoscape, VisANT, etc
# write_delim(brown_list,
#             file = "WGCNA/brown_list.tsv",
#             delim = "\t")

##plotting dedongram
plotDendroAndColors(bwnet_cre$dendrograms[[1]], cbind(bwnet_cre$unmergedColors, bwnet_cre$colors), 
                    c("unmerged", "merged"), 
                    dendroLabels = FALSE,
                    oddGuide = TRUE,
                    hang = 0.03, 
                    guideHang = 0.05)


dendo1 = plotDendroAndColors(bwnet_cre$dendrograms[[1]], cbind(bwnet_cre$colors), 
                             c("merged"), 
                             dendroLabels = FALSE,
                             oddGuide = TRUE,
                             hang = 0.03, 
                             guideHang = 0.05)

##If using module numbers nsted of colours
# plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)


###Finding modules significantly associated with certain conditions
##1.Sexes - has three categories - M, F, W
#coldata2$sex = factor(coldata2$sex, levels = c("M", "F", "W")) #M will nto showup after binarised below

sex.out_cre = binarizeCategoricalColumns(coldata_cre2$sex,
                                     includePairwise = TRUE,
                                     includeLevelVsAll = FALSE)
#minCount = 1


row.names(sex.out_cre) = row.names(coldata_cre2) #need to change rownames when binarizing categorical variables


##3. Tissues
tissues_cre = binarizeCategoricalVariable(coldata_cre2$tissue,
                                      includePairwise = TRUE,
                                      includeLevelVsAll = FALSE) ##o is the reference category
row.names(tissues_cre) = row.names(coldata_cre2) #need to change rownames when binarizing categorical variables


traits_cre = cbind(sex.out_cre, tissues_cre)

##define number of genes and samples

nSample_cre = nrow(norma.counts_cre)
nGenes_cre = ncol(norma.counts_cre)

##calculating corrrealtions between eingenes and traits

##Extra steps if module number is used instead of module colour
# MEs0 = moduleEigengenes(norma.counts, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)

modules.trait.correlation.sex_cre = cor(module_eigengenes_cre, sex.out_cre, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.tissue_cre = cor(module_eigengenes_cre, tissues_cre, use = "p") #correlating eingengenes adn traits

modules.trait.correlation_cre = cor(module_eigengenes_cre, traits_cre, use = "p") 

#Calculating p value for the correlations
modules.trait.corr.pvals.sex_cre = corPvalueStudent(modules.trait.correlation.sex_cre, nSample_cre)
modules.trait.corr.pvals.tissue_cre = corPvalueStudent(modules.trait.correlation.tissue_cre, nSample_cre)

modules.trait.corr.pvals_cre = corPvalueStudent(modules.trait.correlation_cre, nSample_cre)

#visualize module trait association as a heatmap

heatmap_data_cre = merge(module_eigengenes_cre, traits_cre, by = "row.names")
heatmap_data_cre = heatmap_data_cre %>% 
column_to_rownames(var = "Row.names")

colnames(heatmap_data_cre)
CorLevelPlot(heatmap_data_cre,
             x = names(heatmap_data_cre)[11], 
             y = names(heatmap_data_cre)[1:10],
             col = c("red", "pink","white", "pink", "red"),
             #col = c("#da1c5c","#e05a8d", "white","#77c6c6", "#16b2b2"),
             main = "C. remanei SEX")

CorLevelPlot(heatmap_data_cre,
             x = names(heatmap_data_cre)[12], 
             y = names(heatmap_data_cre)[1:10],
             col = c("#004b8d","skyblue", "white", "skyblue", "#004b8d" ),
             #col = c("#557A46","#A6CF98", "white", "#DFA878", "#6C3428" ),
             main = "C. remanei TISSUES")


CorLevelPlot(heatmap_data_cre,
             x = names(heatmap_data_cre)[11:12], 
             y = names(heatmap_data_cre)[1:10],
             col = c("blue1", "skyblue", "white", "pink", "red" ),
             main = "C. remanei")

###Modules in red/blue are significantly associated with one of the sex/species over the other
colnames(heatmap_data_cre)
##Checking genes in a module

module.gene.mapping_cre = as.data.frame(bwnet_cre$colors) #getting genes and theri modules

module.gene.mapping_cre %>%
  filter(`bwnet_cre$colors` == "green") %>%
  rownames() #list of genes in ME green module

###Identifying intramodular hub genes using module membership

module.mem.measure_cre = cor(module_eigengenes_cre, norma.counts_cre, use = "p")

modules.mem.measure.pval_cre = corPvalueStudent(module.mem.measure_cre, nSample_cre)


modules.mem.measure.pval_cre[1:5, 1:5] # depending on the p valu , we can identify whihc genes have sig values for membership meaning can be hub genes

View(t(modules.mem.measure.pval_cre))

##Calculating gene significance associated with difference in sex

gene.sig.cor_cre = cor(norma.counts_cre, sex.out_cre, use = "p")
gene.sig.cor.pval_cre = corPvalueStudent(gene.sig.cor_cre, nSample_cre)


#top gene significantly associated with species trait
gs.sig.sex_cre = gene.sig.cor.pval_cre %>%
  as.data.frame() %>%
  arrange(data.M.vs.F) %>% ##from lowest to highest p-value
  filter(data.M.vs.F < 0.05)

nrow(gs.sig.sex_cre)
head(gs.sig.sex_cre)



eingen_dendo_cre = plotEigengeneNetworks(module_eigengenes_cre, "Cre Eigengene dendrogram", marDendro = c(0,4,2,0),
                                     plotHeatmaps = FALSE)
eingen_heatmap_cre = plotEigengeneNetworks(module_eigengenes_cre, "Cre Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                                       plotDendrograms = FALSE, xLabelsAngle = 90, colorLabels = FALSE)

###############################################################################################




#################################################################################################################
##LOOKING AT MODULES MAKING UP DIFFERENT GENE CATEGORIES

####Make a text file with average module eigengene values across replicates for each module

library(lme4)
library(dplyr)

average_counts_cre = read.table("WGCNA/For WT_without WM/Non-orthologous genes_Cre/Cre_Average module eigenegen value.txt", sep="\t", head=T, comment.char="#")
View(average_counts_cre)

average_counts_cre$Sex = c(rep(c("F", "F", "M", "M"), 1))
average_counts_cre$Tissue = c(rep(c("G", "S"), 2))

str(average_counts_cre)

View(module_eigengenes_cre)
table(bwnet_cre$colors)

####################
##Turquoise
table(bwnet_cre$colors)
turq_cre = average_counts_cre %>% 
  ggplot(aes(Tissue, MEturquoise, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEturquoise N = 1544")

brown_cre = average_counts_cre %>% 
  ggplot(aes(Tissue, MEbrown, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEbrown N = 422")

blue_cre = average_counts_cre %>% 
  ggplot(aes(Tissue, MEblue, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEblue N = 789")

green_cre = average_counts_cre %>% 
  ggplot(aes(Tissue, MEgreen, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEgreen N = 78")

red_cre = average_counts_cre %>% 
  ggplot(aes(Tissue, MEred, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEred N = 45")

yellow_cre = average_counts_cre %>% 
  ggplot(aes(Tissue, MEyellow, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEyellow N = 161")

grey_cre = average_counts_cre %>% 
  ggplot(aes(Tissue, MEgrey, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEgrey N = 397")


################################################################################################
###LOOKING AT GPROFILER ENRICHMENT FOR NON-ORTHOLOGOUS GENES
######################################################################################################333

####Saving genes names to view in GPROFILER 

###Loading cremanie and clatens gene names for older genomes (from Daniel)
remanei_gene_names = read.table("DESeq2/remanei_gene_names.txt", fill = T)
head(remanei_gene_names)

###getting species-biased and unbiased genes from DESeq2 (above) analysis
cre_allcounts ##is the variable with orthologs
cre_allcounts_a = rownames_to_column(cre_allcounts)

cre_1 = subset(remanei_gene_names, remanei_gene_names$V4 %in% cre_allcounts_a$rowname)[,2] #corresponding gene name sin C. latens


##Saving species-baised and unbiased genes
write.table(cre_1, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/cre_nonortho_genenames.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
