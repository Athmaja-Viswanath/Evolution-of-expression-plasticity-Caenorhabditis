##WGCNA for WT DATA ONLY


#####################################################################################################
#####PART II-PAPER 1#################################################################################
#####################################################################################################
library(WGCNA)
library(DESeq2)
#library(GEOquery) for data from NCBI
library(tidyverse)
library(ggplot2)
library(CorLevelPlot)
library(gridExtra)


setwd("C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/")

###Data prep for Wt x WT
orthologs = read.table("new_1to1_orthologgenelist.txt", sep="\t", head=T, comment.char="#")
head(orthologs)
cre_ortho = orthologs[, 1]
clat_ortho = orthologs[ ,2]
cumulative_ortho = orthologs[,3]

head(cre_ortho)
nrow(orthologs)
str(cumulative_ortho)

#C.remanei
cre_counts = read.table("C.remanei_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
cre_counts_ortho = cre_counts[cre_ortho, -c(13:15)]
rownames(cre_counts_ortho) = cumulative_ortho #changing row names
nrow(cre_counts_ortho)

#C.latens
clat_counts = read.table("C.latens_cumulative.txt", sep="\t", head=T, row.name=1, comment.char = "#")
clat_counts_ortho = clat_counts[clat_ortho, -c(13:15)]
rownames(clat_counts_ortho) = cumulative_ortho
nrow(clat_counts)

#Combining all data

counts_orthologs = cbind(cre_counts_ortho, clat_counts_ortho)

colnames(counts_orthologs)
colnames(counts_orthologs) = gsub("X", "", colnames(counts_orthologs)) ##searching and replacing the column names 
colnames(counts_orthologs)

#COLDATA PREP
sample_name = colnames(counts_orthologs)
tissue = substr(sample_name, 3,3) 
#tissue = gsub("M", "W", tissue)
sex = substr(sample_name, 2,2)
species =  c(rep("Cre", 12), rep("Clat", 12))
batch = c(rep(c(1, 2 , 3), 4))
coldata = data.frame(sample_name, species, sex, tissue, batch)
#View(coldata)

###QUALITY CONTROL TO DETECT OUTLIERS

gsg = goodSamplesGenes(t(counts_orthologs))
summary(gsg)

gsg$allOK ##gives outliers, if TRUE = all genes and samples are good and have no outliers

table(gsg$goodGenes) #number of genes that are outliers = FALSE
table(gsg$goodSamples)


counts_orthologs = counts_orthologs[gsg$goodGenes == TRUE, ] ##filtering good genes
dim(counts_orthologs)


####heirarchical clustering to detect outliers

htree = hclust(dist(t(counts_orthologs)), method = "average")
plot(htree) 


##PCA for detecting outliers

pca = prcomp(t(counts_orthologs))
pca_data = pca$x ####PCA calculations for all the samples

#Calculating the variance explained by each principal component

pca.var = pca$sdev^2
pca.var.percent = round(pca.var/sum(pca.var)*100, digits = 2)

pca_data = as.data.frame(pca_data)

ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs( x = paste0("PC1:", pca.var.percent[1], "%"), 
        y = paste0("PC2:", pca.var.percent[2], "%"))

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

coldata2 = coldata[, -1]     ##changing column to rownames
rownames(coldata2) = coldata[ , 1]

all(rownames(coldata2) %in% colnames(counts_orthologs))
all(rownames(coldata2) == colnames(counts_orthologs))

dds = DESeqDataSetFromMatrix(countData = counts_orthologs, 
                             colData = coldata2, 
                             design = ~1) ##not specifiying a model
#View(coldata2)

####removing low count genes

dds75 = dds[rowSums(counts(dds) >= 15) >= 12, ]

nrow(dds75)#11935 genes 


##variance stabilize transformaiton 

dds_norm = vst(dds75)


#get normalized counts
#transfromed data fro downstream analysis

norma.counts = assay(dds_norm) %>% 
  t()

#View(norma.counts)

##### Network construction

##Choose a set of soft threshold powers

power = c(c(1:10), seq(from = 12, to = 50, by = 2))

sft = pickSoftThreshold(norma.counts,
                        powerVector = power,
                        networkType = "signed",
                        verbose = 5)

sft.data = sft$fitIndices ##we will use the R suqared value and mean connectivity , max R^2 and min mean connectivity

###Visualization to pick the right power
a1 = ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()


a2 = ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  #geom_hline(yintercept = 0.8, colour = "red") +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)  ##we need higher R^2 and low mean connectivity


###converting matric to numeric

norma.counts[] = sapply(norma.counts, as.numeric)  ##different step than official tutorial

soft_power = 18

temp_cor = cor #to prevent WGCNa from using other cor function

cor = WGCNA::cor


###Memory estimate wrt blocksize

##can take 20-30 mins to run

bwnet = blockwiseModules(norma.counts, 
                         maxBlockSize = 14000, ##depends on the ram of the system 4gb = 8-10k, 16gb = 20,000, 232gb = 30,000 
                         TOMType = "signed", 
                         power = soft_power,
                         networkType = "signed",
                         mergeCutHeight = 0.25,#threshold that we want to merge similar modules at
                         numericLabels = FALSE, #want the module names to be colours if not, then say TRUE
                         randomSeed = 1234,
                         verbose = 3) 

cor = temp_cor

#5. Module eigenvalue

module_eigengenes = bwnet$MEs
#View(module_eigengenes)

#write.csv(module_eigengenes, file = "Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/module eigengenes values_WT.csv")

###CALCULATE AVERAGE ACROSS REPLICATES FOR EACH MODULE

##get number of genes in each module
table(bwnet$colors)

names(bwnet$colors)

bwnet_df = as.data.frame(bwnet$colors)
bwnet_df = rownames_to_column(bwnet_df)
#View(bwnet_df)

# 
# write.table(bwnet_df, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2//WGCNA_modulegenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)



##############################################################################################################
####Visualize the numebr of genes in each module##################################################################
##############################################################################################################

###loading file which contains genes and their corresponding modules (obtained from above)
bwnet_df = read.table("Results/Plots/Part 2/WGCNA_modulegenes.txt", fill = T)

modules_counts = as.data.frame(table(bwnet_df$V2))

ggplot(modules_counts, aes(x=Var1, y = Freq, fill = Var1))+
  geom_bar(stat = "identity") +
  xlab(" ")+ 
  ylab("Number of genes")+
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5, vjust = 0.5))+
  coord_cartesian(ylim = c(0, 3500))+
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,3500,500))+
  ggtitle("WGCNA Modules")
  

# pick out a few modules of interest here

# module_df <- data.frame(
#   gene_id = names(bwnet$colors),
#   colors = labels2colors(bwnet$colors)) #it doesnt match the right module with the right gene
# 
# #write_delim(module_df, file = "WGCNA/gene_modules.txt", delim = "\t")
# 
# # Pull out list of genes in that module
# 
# modules_of_interest = c("brown")
# 
# genes_of_interest = module_df %>%
#   subset(colors %in% modules_of_interest)
# 
# submod = module_df %>%
#   subset(colors %in% modules_of_interest)
# 
# row.names(module_df) = module_df$gene_id
# 
# # Get normalized expression for those genes
# norma.counts.df = as.data.frame(t(norma.counts))
# brown_list = norma.counts.df[submod$gene_id, ]
# brown_list[1:5, 1:5]
# 
# TOM = TOMsimilarityFromExpr(t(brown_list),
#                             power = soft_power)
# Add gene names to row and columns
# row.names(TOM) = row.names(brown_list)
# colnames(TOM) = row.names(brown_list)
# 
# # brown_list = data.frame(TOM) %>%
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

# Export Network file to be read into Cytoscape, VisANT, etc
# write_delim(brown_list,
#             file = "WGCNA/brown_list.tsv",
#             delim = "\t")

##plotting dedrongram
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors), 
                    c("unmerged", "merged"), 
                    dendroLabels = FALSE,
                    oddGuide = TRUE,
                    hang = 0.03, 
                    guideHang = 0.05)


dendo1 = plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$colors), 
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

##############################################################################################################
###Finding modules significantly associated with certain conditions
##############################################################################################################

##1.Sexes - has three categories - M, F, W
#coldata2$sex = factor(coldata2$sex, levels = c("M", "F", "W")) #M will nto showup after binarised below

sex.out = binarizeCategoricalColumns(coldata2$sex,
                                     includePairwise = TRUE,
                                     includeLevelVsAll = FALSE)
#minCount = 1


row.names(sex.out) = row.names(coldata2) #need to change rownames when binarizing categorical variables

##2. Species
# species = coldata2 %>%
#   mutate(species_binary = ifelse(grepl("Cre", species), 1, 0)) %>% #merges M & W into 0s
#   select(5)

species.factor <- factor(coldata2$species,
                         levels = c('Clat','Cre','H1','H2')) #changing the reference level
species = binarizeCategoricalVariable(coldata2$species,
                                      includePairwise = TRUE,
                                      includeLevelVsAll = FALSE) ##o is the reference category
row.names(species) = row.names(coldata2) #need to change rownames when binarizing categorical variables

##3. Tissues
tissues = binarizeCategoricalVariable(coldata2$tissue,
                                      includePairwise = TRUE,
                                      includeLevelVsAll = FALSE) ##o is the reference category
row.names(tissues) = row.names(coldata2) #need to change rownames when binarizing categorical variables


traits = cbind(sex.out, species, tissues)

##define number of genes and samples

nSample = nrow(norma.counts)
nGenes = ncol(norma.counts)

##############################################################################################################
##calculating corrrealtions between eingenes and traits
##############################################################################################################

##Extra steps if module number is used instead of module colour
# MEs0 = moduleEigengenes(norma.counts, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)
modules.trait.correlation.sp = cor(module_eigengenes, species, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.sex = cor(module_eigengenes, sex.out, use = "p") #correlating eingengenes adn traits
modules.trait.correlation.tissue = cor(module_eigengenes, tissues, use = "p") #correlating eingengenes adn traits

modules.trait.correlation = cor(module_eigengenes, traits, use = "p") 

#Calculating p value for the correlations
modules.trait.corr.pvals.sp = corPvalueStudent(modules.trait.correlation.sp, nSample)
modules.trait.corr.pvals.sex = corPvalueStudent(modules.trait.correlation.sex, nSample)
modules.trait.corr.pvals.tissue = corPvalueStudent(modules.trait.correlation.tissue, nSample)

modules.trait.corr.pvals = corPvalueStudent(modules.trait.correlation, nSample)

##############################################################################################################
#visualize module trait association as a heatmap
##############################################################################################################

heatmap_data = merge(module_eigengenes, traits, by = "row.names")
heatmap_data = heatmap_data %>% 
  column_to_rownames(var = "Row.names")

colnames(heatmap_data)
CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[13], 
             y = names(heatmap_data)[1:12],
             col = c("red", "pink","white", "pink", "red"),
             main = "SEX")

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[14], 
             y = names(heatmap_data)[1:12],
             col = c("#ffeb3b","#fff280", "white", "#fff280", "#ffeb3b"),
             main = "SPECIES")

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[15], 
             y = names(heatmap_data)[1:12],
             col = c("#004b8d","skyblue", "white", "skyblue", "#004b8d" ),
             main = "TISSUES")

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[13:15], 
             y = names(heatmap_data)[1:12],
             col = c("blue1", "skyblue", "white", "pink", "red" ))

###Modules in red/blue are significantly associated with one of the sex/species over the other
colnames(heatmap_data)

##############################################################################################################
##Checking genes in a module

#module.gene.mapping = as.data.frame(bwnet$colors) #getting genes and their modules
module.gene.mapping = bwnet_df %>% column_to_rownames(var = "V1")
module.gene.mapping %>%
  filter(V2 == "green") %>%
  rownames() #list of genes in ME green module

##############################################################################################################
###Identifying intramodular hub genes using module membership
##############################################################################################################

module.mem.measure = cor(module_eigengenes, norma.counts, use = "p")

modules.mem.measure.pval = corPvalueStudent(module.mem.measure, nSample)

module.mem.measure = as.data.frame(module.mem.measure) #need to be transposed

modules.mem.measure.pval[1:5, 1:5] # depending on the p value , we can identify whihc genes have sig values for membership meaning can be hub genes
# 
# View(t(modules.mem.measure.pval))
# View(t(module.mem.measure))

##############################################################################################################
##Calculating gene significance associated with difference in species
##############################################################################################################

gene.sig.cor = cor(norma.counts, species, use = "p")
gene.sig.cor.pval = corPvalueStudent(gene.sig.cor, nSample)
gene.sig.cor = as.data.frame(gene.sig.cor)
#view(gene.sig.cor)
#top gene significantly associated with species trait
gs.sig.sp = gene.sig.cor.pval %>%
  as.data.frame() %>%
  arrange(Cre.vs.Clat) %>% ##from lowest to highest p-value
  filter(Cre.vs.Clat < 0.05)

nrow(gs.sig.sp)
nrow(gene.sig.cor)

##########################################################################################################################
##################IDENTIFYING HUB GENES IN YELLOW AND BLACK MODULES CORRELATED WITH SPECIES###############################

#getting genes in yellow module
yellow_module = module.gene.mapping %>% filter(V2 == "yellow") %>% rownames()

##Getting module membership values for genes in yellow module (M10)
#View(t(module.mem.measure[2, yellow_module])) #since module.mem.measure is transposed, we have rows as module names and columns as gene names

mem_yellow = as.data.frame(t(module.mem.measure[2,yellow_module]))
gs_yellow = as.data.frame(abs(gene.sig.cor[yellow_module, 1]))

yelow_gs_meme = cbind(mem_yellow, gs_yellow)

#creating a new column based on other column values
yelow_gs_meme$category = with(yelow_gs_meme, ifelse(MEyellow>0.9 & `abs(gene.sig.cor[yellow_module, 1])`>0.9, "hub genes", "non-hub"))
#View(yelow_gs_meme)

#####Gene significance vs module member ship plot
ggplot(yelow_gs_meme, aes(x=yelow_gs_meme$MEyellow, y=yelow_gs_meme$`abs(gene.sig.cor[yellow_module, 1])`, colour = yelow_gs_meme$category)) +
  geom_point(alpha = 0.6, show.legend = FALSE)+
  coord_cartesian(ylim = c(0, 1))+
  coord_cartesian(xlim = c(0, 1))+
  geom_hline(yintercept = 0.9, colour = "red")+
  geom_vline(xintercept = 0.9, colour = "red")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))+
  ggtitle("Gene Significance vs Module membership for M10")

# verboseScatterplot((t(module.mem.measure[2,yellow_module])),
#                    abs(gene.sig.cor[yellow_module, 1]),
#                    xlab = paste("Module Membership in", "Yellow", "module"),
#                    ylab = "Gene significance for Yellow Module",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")


####Number of hubgenes in yellwo/M10 module

table(yelow_gs_meme$category)

# write.table(yelow_gs_meme, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/yelow_gs_meme.txt",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)

##############################################################################################################
########################GETTING GENE LISTS FOR HUB GENES TO SEE IN GPROFILER
##############################################################################################################

remanei_gene_names = read.table("DESeq2/remanei_gene_names.txt", fill = T)
latens_gene_names = read.table("DESeq2/latens_gene_names.txt", fill = T)

str(yelow_gs_meme)

yellow_genes = subset(yelow_gs_meme, yelow_gs_meme$category=="hub genes") %>%  rownames_to_column() %>%
       separate(rowname, sep = 23, into = c("Crenames", "Clatnames"))
yellow_genes = separate(yellow_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]

yellow_genes_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% yellow_genes$Crenames)[,2]
yellow_genes_clat = subset(latens_gene_names, latens_gene_names$V4 %in% yellow_genes$Clatnames)[,2]

length(yellow_genes_cre)

# write.table(yellow_genes_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/M10_yellow_crehubgenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.table(yellow_genes_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/M10_yellow_clathubgenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################################################################
########OVERLAP BETWEEN HUB GENES in M10/YELLOW AND SPECIES BIASED GENES#####################################
##############################################################################################################

cre_biased = read.table("DESeq2/crexclat_up.txt") %>% row.names()
clat_biased = read.table("DESeq2/crexclat_down.txt") %>% row.names()
M10_hubgenes = subset(yelow_gs_meme, yelow_gs_meme$category=="hub genes") %>% row.names()
M10_allgenes  = yelow_gs_meme %>% row.names()
#View(M10_hubgenes)


################################################################################################3
######FIGURE 2D################################################################################
################################################################################################3

library(VennDiagram)

temp = venn.diagram(list("Clat_biased" =  clat_biased, "M10 hub genes" = M10_hubgenes, "M10 all genes" = M10_allgenes),filename = NULL)
temp = draw.triple.venn(979, 3339, 271, 952, 271, 271, 271, fileName = NULL)

overrideTriple = TRUE
draw.triple.venn(
  area1=200, area2=300, area3=400, n12=20, n23=30, n13=40, n123=10, 
  category=c("#1","#2","#3"),
  col="Black",fill=c("White","Red","Yellow"),
  cex=3, cat.cex=2, euler.d=TRUE, scaled=TRUE)

draw.triple.venn(
  area1=979, area2=3339, area3=271, n12=952, n23=271, n13=271, n123=271, 
  category=c("#1","#2","#3"),
  col="Black",fill=c("White","Red","Yellow"),
  cex=3, cat.cex=2, euler.d=TRUE, scaled=TRUE)
Sys.sleep(20)
grid.draw(temp)
pdf(file="Results/Plots/Part 2/DEGvsHubgenes_M10.pdf")
  grid.draw(temp)
dev.off()

pdf(file = "3-Output/DEGvsHubgenes_M10.pdf")
grid.draw(temp)
dev.off()
##############################################################################################################
####getting files for cytoscape##############################

head(subset(bwnet_df, bwnet_df$`bwnet$colors` == "yellow"))
head(bwnet_df)

row.names(bwnet_df )= bwnet_df$rowname

##############################################################################################################
# Get normalized expression for those genes
norma.counts.df = as.data.frame(t(norma.counts))
y_m10_list = norma.counts.df[M10_hubgenes, ]
y_m10_list[1:5, 1:5]
View(norma.counts.df)
# 
TOM_yellow_M10 = TOMsimilarityFromExpr(t(y_m10_list),
                             power = soft_power,
                             networkType = "signed", 
                             TOMType = "signed")

# # Add gene names to row and columns
row.names(TOM_yellow_M10) = row.names(y_m10_list)
colnames(TOM_yellow_M10) = row.names(y_m10_list)
head(TOM_yellow_M10)

y_m10_list = data.frame(TOM_yellow_M10) %>%
   mutate(
     gene1 = row.names(.)
   ) %>%
   pivot_longer(-gene1) %>%
   dplyr::rename(gene2 = name, correlation = value) %>%
   unique() %>%
   subset(!(gene1==gene2)) %>%
   mutate(
     module1 = bwnet_df[gene1,]$'bwnet$colors',
     module2 = bwnet_df[gene2,]$'bwnet$colors'
   )

View(y_m10_list)

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(y_m10_list,
            file = "Results/Plots/Part 2/Cytoscape/Yellow_hub_M10.tsv",
            delim = "\t")

########################################################################################################
####VISUALIZING MEM VS GS FOR BLACK MODULE
##############################################################################################################

module = "MEblack"
#getting genes in black module
black_module = module.gene.mapping %>% filter(`bwnet$colors` == "black") %>% rownames()

length(black_module)

##############################################################################################################
##Getting module membership values for genes in black module (M10)

View(t(module.mem.measure[module, black_module])) #since module.mem.measure is transposed, we have rows as module names and columns as gene names

mem_black = as.data.frame(t(module.mem.measure[module,black_module]))
gs_black = as.data.frame(abs(gene.sig.cor[black_module, 1]))

head(black_gs_meme)
black_gs_meme = cbind(mem_black, gs_black)

#creating a new column based on other column values
black_gs_meme$category = with(black_gs_meme, ifelse(MEblack>0.9 & `abs(gene.sig.cor[black_module, 1])`>0.9, "hub genes", "non-hub"))
View(black_gs_meme)

##############################################################################################################
#####Gene significance vs module member ship plot
ggplot(black_gs_meme, aes(x=black_gs_meme$MEblack, y=black_gs_meme$`abs(gene.sig.cor[black_module, 1])`, colour = black_gs_meme$category)) +
  geom_point(alpha = 0.6, show.legend = FALSE)+
  coord_cartesian(ylim = c(0, 1),xlim = c(0, 1))+
  geom_hline(yintercept = 0.9, colour = "red")+
  geom_vline(xintercept = 0.9, colour = "red")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))+
  ggtitle("Gene Significance vs Module membership for M4")

# verboseScatterplot((t(module.mem.measure[module,black_module])),
#                    abs(gene.sig.cor[black_module, 1]),
#                    xlab = paste("Module Membership in", "black", "module"),
#                    ylab = "Gene significance for black Module",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")

##############################################################################################################
####Number of hubgenes in black/M4 module

table(black_gs_meme$category)

# write.table(black_gs_meme, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/black_gs_meme.txt",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)

##############################################################################################################
########################GETTING GENE LISTS FOR HUB GENES TO SEE IN GPROFILER
##############################################################################################################

remanei_gene_names = read.table("DESeq2/remanei_gene_names.txt", fill = T)
latens_gene_names = read.table("DESeq2/latens_gene_names.txt", fill = T)


table(black_gs_meme$category)

black_genes = subset(black_gs_meme, black_gs_meme$category=="hub genes") %>%  rownames_to_column() %>%
  separate(rowname, sep = 23, into = c("Crenames", "Clatnames"))
black_genes = separate(black_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]

black_genes_cre = subset(remanei_gene_names, remanei_gene_names$V4 %in% black_genes$Crenames)[,2]
black_genes_clat = subset(latens_gene_names, latens_gene_names$V4 %in% black_genes$Clatnames)[,2]

length(black_genes_cre)
length(black_genes_clat)

# write.table(black_genes_cre, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/M4_black_crehubgenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.table(black_genes_clat, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/M4_black_clathubgenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################################################################
########OVERLAP BETWEEN HUB GENES in M4/BLACK AND SPECIES BIASED GENES#####################################
##############################################################################################################

cre_biased = read.table("DESeq2/crexclat_up.txt") %>% row.names()
clat_biased = read.table("DESeq2/crexclat_down.txt") %>% row.names()
M4_hubgenes = subset(black_gs_meme, black_gs_meme$category=="hub genes") %>% row.names()
M4_allgenes  = black_gs_meme %>% row.names()
View(M10_hubgenes)

################################################################################################3
######FIGURE 2F################################################################################
################################################################################################3

library(VennDiagram)

temp = venn.diagram(list("Cre_biased" =  cre_biased, "M4 hub genes" = M4_hubgenes, "M4 all genes" = M4_allgenes), 
                    filename = NULL )

grid.draw(temp)
pdf(file="Results/Plots/Part 2/DEGvsHubgenes_M4.pdf")
grid.draw(temp)
dev.off()

##########################################################################################################################

#######################################################################################################################

#Calculating gene significance associated with difference in sex

gene.sig.cor.sex = cor(norma.counts, sex.out, use = "p")
gene.sig.cor.pval.sex = corPvalueStudent(gene.sig.cor.sex, nSample)

head(gene.sig.cor.pval.sex)
#top gene significantly associated with species trait
gs.sig.sex = gene.sig.cor.pval.sex %>%
  as.data.frame() %>%
  arrange(data.M.vs.F) %>% ##from lowest to highest p-value
  filter(data.M.vs.F < 0.05)

nrow(gs.sig.sex)

#######DENDROGRAMS AND ADJACENCY MATRIX#######################################################################################
eingen_dendo = plotEigengeneNetworks(module_eigengenes, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE, excludeGrey = FALSE)
eingen_heatmap = plotEigengeneNetworks(module_eigengenes, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90, colorLabels = FALSE)

# write.table(bwnet_df, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/WGCNA_modulegenes.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

##################################################################################################################
##VISUALIZING MODULES##############################################################################################################
##############################################################################################################
####Make a text file with average module eigen gene values across replicates for each module

library(lme4)
library(dplyr)

average_counts = read.table("Results/Plots/Part 2/Average_MEvalues.txt", sep="\t", head=T, comment.char="#")
View(average_counts)
#colnames(average_counts) = c("Samples", "M12", "M11", "M10", "M9", "M8", "M7", 
#                             "M6", "M5", "M4", "M3", "M2", "M1", "M0")

average_counts$Species = c(rep("Cre", 4), rep("Clat", 4))
average_counts$Sex = c(rep(c("F", "F", "M", "M"), 2))
average_counts$Tissue = c(rep(c("G", "S"), 4))
average_counts$S_Sp = c("Cre_F", "Cre_F", "Cre_M", "Cre_M", "Clat_F", "Clat_F", "Clat_M", "Clat_M")

str(average_counts)

#View(module_eigengenes)

##Turquoise
table(bwnet$colors)
average_counts %>% 
  ggplot(aes(Tissue, MEturquoise, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M5 N = 3156")


##BLUE

average_counts %>% 
  ggplot(aes(Tissue, MEblue, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M9 N = 2761")

###GREENYELLOW
average_counts %>% 
  ggplot(aes(Tissue, MEgreenyellow, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M3 N = 133")

###PURPLE

average_counts %>% 
  ggplot(aes(Tissue, MEpurple, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M11 N = 164")

######RED
average_counts %>% 
  ggplot(aes(Tissue, MEred, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M8 N = 652")

##BLACK

average_counts %>% 
  ggplot(aes(Tissue, MEblack, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M4 N = 583")

###GREEN

average_counts %>%
  ggplot(aes(Tissue, MEgreen, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M6 N = 822")


###YELLOW

average_counts %>% 
  ggplot(aes(Tissue, MEyellow, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M10 N = 979")

###PINK

average_counts %>% 
  ggplot(aes(Tissue, MEpink, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M7 N = 398")

####BROWN

average_counts %>% 
  ggplot(aes(Tissue, MEbrown, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M2 N = 1770")

###MAGENTA

average_counts %>% 
  ggplot(aes(Tissue, MEmagenta, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M1 N = 274")


###$###GREY
average_counts %>% 
  ggplot(aes(Tissue, MEgrey, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M0 N = 112")

######################################################################################################################
###GETTIN GENES IN EACH MODULE##################################################################################
##############################################################################################################

remanei_gene_names = read.table("DESeq2/remanei_gene_names.txt", fill = T)
latens_gene_names = read.table("DESeq2/latens_gene_names.txt", fill = T)

table(bwnet$colors)
names(bwnet$colors)
bwnet_df = as.data.frame(bwnet$colors)
bwnet_df = rownames_to_column(bwnet_df)

#Magenta MODULE 1
M1_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "magenta", ])
M1_genes = separate(M1_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M1_genes = separate(M1_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M1_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M1_genes$Crenames)[,2]
M1_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M1_genes$Clatnames)[,2]

write.table(M1_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M1_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M1_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M1_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#BROWN MODULE 2

M2_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "brown", ])
M2_genes = separate(M2_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M2_genes = separate(M2_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M2_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M2_genes$Crenames)[,2]
M2_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M2_genes$Clatnames)[,2]

length(M2_genes$rowname)
write.table(M2_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M2_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M2_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M2_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


#greenyellow MODULE 3
M3_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "greenyellow", ])
M3_genes = separate(M3_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M3_genes = separate(M3_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M3_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M3_genes$Crenames)[,2]
M3_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M3_genes$Clatnames)[,2]

length(M3_Cre_genes)
write.table(M3_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M3_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M3_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M3_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#BLACK MODULE 4
M4_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "black", ])
M4_genes = separate(M4_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M4_genes = separate(M4_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M4_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M4_genes$Crenames)[,2]
M4_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M4_genes$Clatnames)[,2]

length(M4_Cre_genes)
write.table(M4_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M4_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M4_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M4_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##Turquoise M5
M5_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "turquoise", ])
M5_genes = separate(M5_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M5_genes = separate(M5_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M5_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M5_genes$Crenames)[,2]
M5_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M5_genes$Clatnames)[,2]

length(M5_Cre_genes)
write.table(M5_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M5_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M5_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M5_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#green MODULE 6
M6_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "green", ])
M6_genes = separate(M6_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M6_genes = separate(M6_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M6_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M6_genes$Crenames)[,2]
M6_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M6_genes$Clatnames)[,2]

length(M6_Cre_genes)
write.table(M6_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M6_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M6_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M6_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#pink MODULE 7
M7_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "pink", ])
M7_genes = separate(M7_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M7_genes = separate(M7_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M7_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M7_genes$Crenames)[,2]
M7_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M7_genes$Clatnames)[,2]

length(M7_Cre_genes)
write.table(M7_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M7_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M7_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M7_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#red MODULE 8
M8_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "red", ])
M8_genes = separate(M8_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M8_genes = separate(M8_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M8_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M8_genes$Crenames)[,2]
M8_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M8_genes$Clatnames)[,2]

length(M8_Cre_genes)
write.table(M8_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M8_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M8_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M8_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#blue MODULE 9
M9_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "blue", ])
M9_genes = separate(M9_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M9_genes = separate(M9_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M9_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M9_genes$Crenames)[,2]
M9_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M9_genes$Clatnames)[,2]

length(M9_Cre_genes)
write.table(M9_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M9_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M9_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M9_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#yellow MODULE 10
M10_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "yellow", ])
M10_genes = separate(M10_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M10_genes = separate(M10_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M10_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M10_genes$Crenames)[,2]
M10_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M10_genes$Clatnames)[,2]

length(M10_Cre_genes)
write.table(M10_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M10_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M10_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M10_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#purple MODULE 11
M11_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "purple", ])
M11_genes = separate(M11_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M11_genes = separate(M11_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M11_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M11_genes$Crenames)[,2]
M11_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M11_genes$Clatnames)[,2]

length(M11_Cre_genes)
write.table(M11_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M11_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M11_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M11_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#grey MODULE 0
M0_genes = as.data.frame(bwnet_df[bwnet_df$`bwnet$colors` == "grey", ])
M0_genes = separate(M0_genes, rowname, sep = 23, into = c("Crenames", "Clatnames"))
M0_genes = separate(M0_genes, Clatnames, sep = 1, into = c("random", "Clatnames"))[,-2]
M0_Cre_genes  = subset(remanei_gene_names, remanei_gene_names$V4 %in% M0_genes$Crenames)[,2]
M0_Clat_genes  = subset(latens_gene_names, latens_gene_names$V4 %in% M0_genes$Clatnames)[,2]

length(M0_Cre_genes)
write.table(M0_Cre_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M0_cre_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M0_Clat_genes, file = "C:/Users/athma/Desktop/RNASeq_results/DGE analysis and Data/New annotation/WGCNA/For WT_without WM/Orthologous genes/Module gene list/M0_clat_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)



