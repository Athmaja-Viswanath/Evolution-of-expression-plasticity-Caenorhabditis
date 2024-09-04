###########################################################################################
######################PART IIIB############################################################
###########################################################################################

##WGCNA for C. latens non-orhtologous genes ONLY

library(WGCNA)
library(DESeq2)
#library(GEOquery) for data from NCBI
library(tidyverse)
library(ggplot2)
library(CorLevelPlot)
library(gridExtra)

setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

###Data prep for Wt x WT

#C.latens
clat_allcounts = read.table("C.latens_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
clat_allcounts = rownames_to_column(clat_allcounts)
"%notin%" = Negate("%in%")
clat_allcounts = clat_allcounts[clat_allcounts$rowname %notin% orthologs$C.latens.Gene.name, ]
nrow(clat_allcounts)
clat_allcounts = remove_rownames(clat_allcounts)
clat_allcounts = column_to_rownames(clat_allcounts, var = "rowname")
clat_allcounts = clat_allcounts[ ,-c(13:15)]
names(clat_allcounts)
colnames(clat_allcounts) = gsub("X", "", colnames(clat_allcounts)) ##searching and replacing the column names 
colnames(clat_allcounts)

#COLDATA PREP
sample_name_clat = colnames(clat_allcounts)
tissue_clat = substr(sample_name_clat, 3,3) 
sex_clat = substr(sample_name_clat, 2,2)
batch_clat = c(rep(c(1, 2 , 3), 4))
coldata_clat = data.frame(sample_name_clat, sex_clat, tissue_clat, batch_clat)
View(coldata_clat)


head(clat_allcounts)

###QUALITY CONTROL TO DETECT OUTLIERS

gsg_clat = goodSamplesGenes(t(clat_allcounts))
summary(gsg_clat)

gsg_clat$allOK ##gives outliers, if TRUE = all genes and samples are good and have no outliers

table(gsg_clat$goodGenes) #number of genes that are outliers = FALSE
table(gsg_clat$goodSamples)


clat_allcounts = clat_allcounts[gsg_clat$goodGenes == TRUE, ] ##filtering good genes
nrow(clat_allcounts)


####heirarchical clustering to detect outliers

htree_clat = hclust(dist(t(clat_allcounts)), method = "average")
plot(htree_clat) 


##PCA for detecting outliers

pca_clat = prcomp(t(clat_allcounts))
pca_data_clat = pca_clat$x ####PCA calculations for all the samples

#Calculating the variance explained by each principal component

pca.var_clat = pca_clat$sdev^2
pca.var.percent_clat = round(pca.var_clat/sum(pca.var_clat)*100, digits = 2)

pca_data_clat = as.data.frame(pca_data_clat)

ggplot(pca_data_clat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data_clat)) +
  labs( x = paste0("PC1:", pca.var.percent_clat[1], "%"), 
        y = paste0("PC2:", pca.var.percent_clat[2], "%"))

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

coldata_clat2 = coldata_clat[, -1]     ##changing column to rownames
rownames(coldata_clat2) = coldata_clat[ , 1]

all(rownames(coldata_clat2) %in% colnames(clat_allcounts))
all(rownames(coldata_clat2) == colnames(clat_allcounts))

dds_clat = DESeqDataSetFromMatrix(countData = clat_allcounts, 
                                 colData = coldata_clat2, 
                                 design = ~1) ##not specifiying a model
View(coldata_clat2)

####removing low count genes

dds75_clat = dds_clat[rowSums(counts(dds_clat) >= 15) >= 6, ]

nrow(dds75_clat) 


##variance stabilize transformaiton 

dds_norm_clat = vst(dds75_clat)


#get normalized counts
#transfromed data fro downstream analysis

norma.counts_clat = assay(dds_norm_clat) %>% 
  t()

#View(norma.counts)


#####4. Network construction

##Choose a set of soft threshold powers

power = c(c(1:10), seq(from = 12, to = 50, by = 2))

sft_clat = pickSoftThreshold(norma.counts_clat,
                            powerVector = power,
                            networkType = "signed",
                            verbose = 5)

sft.data_clat = sft_clat$fitIndices ##we will use the R suqared value and mean connectivity , max R^2 and min mean connectivity

###Visualization to pick the right power
a1 = ggplot(sft.data_clat, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()


a2 = ggplot(sft.data_clat, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  #geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)  ##we need higher R^2 and low mean connectivity


###converting matric to numeric

norma.counts_clat[] = sapply(norma.counts_clat, as.numeric)  ##different step than official tutorial

soft_power_clat = 26

temp_cor = cor #to prevent WGCNa from using other cor function

cor = WGCNA::cor


###Memory estimate wrt blocksize

##can take 20-30 mins to run

bwnet_clat = blockwiseModules(norma.counts_clat, 
                             maxBlockSize = 14000, ##depends on the ram of the system 4gb = 8-10k, 16gb = 20,000, 232gb = 30,000 
                             TOMType = "signed", 
                             power = soft_power_clat,
                             networkType = "signed",
                             mergeCutHeight = 0.25,#threshold that we want to merge similar modules at
                             numericLabels = FALSE, #want the module names to be colours if not, then say TRUE
                             randomSeed = 1234,
                             verbose = 3) 

cor = temp_cor

#5. Module eigenvalue

module_eigengenes_clat = bwnet_clat$MEs
View(module_eigengenes_clat)


#write.csv(module_eigengenes_clat, file = "WGCNA//For WT_without WM/Non-orthologous genes_Clat//C.latens module eigengenes values.csv")


###CALCULATE AVERAGE ACROSS REPLICATES FOR EACH MODULE

##get number of genes in each module
table(bwnet_clat$colors)

names(bwnet_clat$colors)

##############################################################################################################
####Visualize the numebr of genes in each module##################################################################
bwnet_df_clat = as.data.frame(bwnet_clat$colors)
modules_counts_clat = as.data.frame(table(bwnet_df_clat$`bwnet_clat$colors`))

ggplot(modules_counts_clat, aes(x=Var1, y = Freq, fill = Var1))+
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
  ggtitle("Clat-specific WGCNA Modules")
# pick out a few modules of interest here

# module_df_clat <- data.frame(
#   gene_id = names(bwnet_clat$colors),
#   colors = labels2colors(bwnet_clat$colors)) #it doesnt match the right module with the right gene

#write_delim(module_df, file = "WGCNA/gene_modules.txt", delim = "\t")

# Pull out list of genes in that module

# modules_of_interest_clat = c("brown")
# 
# genes_of_interest_clat = module_df_clat %>%
#   subset(colors %in% modules_of_interest_clat)
# 
# submod_clat = module_df_clat %>%
#   subset(colors %in% modules_of_interest_clat)
# 
# row.names(module_df_clat) = module_df_clat$gene_id
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
plotDendroAndColors(bwnet_clat$dendrograms[[1]], cbind(bwnet_clat$unmergedColors, bwnet_clat$colors), 
                    c("unmerged", "merged"), 
                    dendroLabels = FALSE,
                    oddGuide = TRUE,
                    hang = 0.03, 
                    guideHang = 0.05)


dendo1_clat = plotDendroAndColors(bwnet_clat$dendrograms[[1]], cbind(bwnet_clat$colors), 
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

sex.out_clat = binarizeCategoricalColumns(coldata_clat2$sex,
                                         includePairwise = TRUE,
                                         includeLevelVsAll = FALSE)
#minCount = 1


row.names(sex.out_clat) = row.names(coldata_clat2) #need to change rownames when binarizing categorical variables


##3. Tissues
tissues_clat = binarizeCategoricalVariable(coldata_clat2$tissue,
                                          includePairwise = TRUE,
                                          includeLevelVsAll = FALSE) ##o is the reference category
row.names(tissues_clat) = row.names(coldata_clat2) #need to change rownames when binarizing categorical variables


traits_clat = cbind(sex.out_clat, tissues_clat)

##define number of genes and samples

nSample_clat = nrow(norma.counts_clat)
nGenes_clat = ncol(norma.counts_clat)

##calculating corrrealtions between eingenes and traits

##Extra steps if module number is used instead of module colour
# MEs0 = moduleEigengenes(norma.counts, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)

modules.trait.correlation.sex_clat = cor(module_eigengenes_clat, sex.out_clat, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.tissue_clat = cor(module_eigengenes_clat, tissues_clat, use = "p") #correlating eingengenes adn traits

modules.trait.correlation_clat = cor(module_eigengenes_clat, traits_clat, use = "p") 

#Calculating p value for the correlations
modules.trait.corr.pvals.sex_clat = corPvalueStudent(modules.trait.correlation.sex_clat, nSample_clat)
modules.trait.corr.pvals.tissue_clat = corPvalueStudent(modules.trait.correlation.tissue_clat, nSample_clat)

modules.trait.corr.pvals_clat = corPvalueStudent(modules.trait.correlation_clat, nSample_clat)

#visualize module trait association as a heatmap

heatmap_data_clat = merge(module_eigengenes_clat, traits_clat, by = "row.names")
heatmap_data_clat = heatmap_data_clat %>% 
  column_to_rownames(var = "Row.names")

colnames(heatmap_data_clat)

CorLevelPlot(heatmap_data_clat,
             x = names(heatmap_data_clat)[16], 
             y = names(heatmap_data_clat)[1:15],
             #col = c("#da1c5c","#e05a8d", "white","#77c6c6", "#16b2b2"),
             main = "C. latens SEX", 
             col = c("red", "pink","white", "pink", "red"))

CorLevelPlot(heatmap_data_clat,
             x = names(heatmap_data_clat)[17], 
             y = names(heatmap_data_clat)[1:15],
             col = c("#004b8d","skyblue", "white", "skyblue", "#004b8d" ),
             #col = c("#557A46","#A6CF98", "white", "#DFA878", "#6C3428" ),
             main = "C. latens TISSUES")


CorLevelPlot(heatmap_data_clat,
             x = names(heatmap_data_clat)[16:17], 
             y = names(heatmap_data_clat)[1:15],
             col = c("blue1", "skyblue", "white", "pink", "red" ), 
             main = "C. latens")

###Modules in red/blue are significantly associated with one of the sex/species over the other
colnames(heatmap_data_clat)
##Checking genes in a module

module.gene.mapping_clat = as.data.frame(bwnet_clat$colors) #getting genes and theri modules

module.gene.mapping_clat %>%
  filter(`bwnet_clat$colors` == "green") %>%
  rownames() #list of genes in ME green module

###Identifying intramodular hub genes using module membership

module.mem.measure_clat = cor(module_eigengenes_clat, norma.counts_clat, use = "p")

modules.mem.measure.pval_clat = corPvalueStudent(module.mem.measure_clat, nSample_clat)


modules.mem.measure.pval_clat[1:5, 1:5] # depending on the p valu , we can identify whihc genes have sig values for membership meaning can be hub genes

View(t(modules.mem.measure.pval_clat))

##Calculating gene significance associated with difference in species

gene.sig.cor_clat = cor(norma.counts_clat, sex.out_clat, use = "p")
gene.sig.cor.pval_clat = corPvalueStudent(gene.sig.cor_clat, nSample_clat)


#top gene significantly associated with species trait
gs.sig.sex_clat = gene.sig.cor.pval_clat %>%
  as.data.frame() %>%
  arrange(data.M.vs.F) %>% ##from lowest to highest p-value
  filter(data.M.vs.F < 0.05)

nrow(gs.sig.sex_clat)
head(gs.sig.sex_clat)



eingen_dendo_clat = plotEigengeneNetworks(module_eigengenes_clat, "Clat Eigengene dendrogram", marDendro = c(0,4,2,0),
                                         plotHeatmaps = FALSE)
eingen_heatmap_clat = plotEigengeneNetworks(module_eigengenes_clat, "Clat Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                                           plotDendrograms = FALSE, xLabelsAngle = 90, colorLabels = FALSE)


#################################################################################################################
##LOOKING AT MODULES MAKING UP DIFFERENT GENE CATEGORIES

####Make a text file with average module eigengene values across replicates for each module

library(lme4)
library(dplyr)

average_counts_clat = read.table("WGCNA/For WT_without WM/Non-orthologous genes_Clat/Clat_average module eigengene values.txt", sep="\t", head=T, comment.char="#")
View(average_counts_clat)

average_counts_clat$Sex = c(rep(c("F", "F", "M", "M"), 1))
average_counts_clat$Tissue = c(rep(c("G", "S"), 1))

str(average_counts_clat)

View(module_eigengenes_clat)
View(table(bwnet_clat$colors))

####################
##Turquoise

turq_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEturquoise, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEturquoise N = 1694")

brown_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEbrown, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEbrown N = 565")

blue_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEblue, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEblue N = 825")

green_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEgreen, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEgreen N = 205")

red_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEred, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEred N = 185")

yellow_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEyellow, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEyellow N = 264")

grey_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEgrey, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEgrey N = 610")

gyellow_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEgreenyellow, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEgreenyellow N = 28")

magenta_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEmagenta, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEmagenta N = 97")


yellow_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEyellow, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEyellow N = 264")


pink_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEpink, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEpink N = 112")

purple_clat = average_counts_clat %>% 
  ggplot(aes(Tissue, MEpurple, col = Sex, group = Sex)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  stat_smooth(method = "lm") +
  ggtitle("MEpurple N = 88")


################################################################################################
###LOOKING AT GPROFILER ENRICHMENT FOR NON-ORTHOLOGOUS GENES
######################################################################################################333

####Saving genes names to view in GPROFILER 

###Loading cremanie and clatens gene names for older genomes (from Daniel)
latens_gene_names = read.table("DESeq2/latens_gene_names.txt", fill = T)
head(latens_gene_names)

###getting species-biased and unbiased genes from DESeq2 (above) analysis
clat_allcounts ##is the variable with orthologs
clat_allcounts_a = rownames_to_column(clat_allcounts)

clat_1 = subset(latens_gene_names, latens_gene_names$V4 %in% clat_allcounts_a$rowname)[,2] #corresponding gene name sin C. latens


##Saving species-baised and unbiased genes
# write.table(clat_1, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/DESeq2/clat_nonortho_genenames.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)




