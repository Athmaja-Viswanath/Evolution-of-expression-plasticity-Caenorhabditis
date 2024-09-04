##################################################################################################################
##VISUALIZING MODULES##############################################################################################################
##############################################################################################################
####Make a text file with average module eigen gene values across replicates for each module

#Packages to be loaded
library(lme4)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(gridExtra)

###Save the module eingen gene values in csv/txt once you get them (usign the command below)
#module_eigengenes = bwnet$MEs
#View(module_eigengenes)

#write.csv(module_eigengenes, file = "Desktop/RNASeq_results/DGE analysis and Data/New annotation/Results/Plots/Part 2/module eigengenes values_WT.csv")

##In your case, you should have one value for each sample for each module


###I had to calculate average (ie., the average_counts file below) as I had 3 replicates/sample 
#but you dont need to do that, you shoudl be directly able to use the excel sheet that you saved (above)

#read in the file, you dont have to run wgcna again if you saved the excel file
average_counts = read.table("Results/Plots/Part 2/Average_MEvalues.txt", sep="\t", head=T, comment.char="#")
View(average_counts)
#colnames(average_counts) = c("Samples", "M12", "M11", "M10", "M9", "M8", "M7", 
#                             "M6", "M5", "M4", "M3", "M2", "M1", "M0")

###Adding columns to represent sample information (you will have the different treatment names)
average_counts$Species = c(rep("Cre", 4), rep("Clat", 4))
average_counts$Sex = c(rep(c("F", "F", "M", "M"), 2))
average_counts$Tissue = c(rep(c("G", "S"), 4))
average_counts$S_Sp = c("Cre_F", "Cre_F", "Cre_M", "Cre_M", "Clat_F", "Clat_F", "Clat_M", "Clat_M")

str(average_counts)

#View(module_eigengenes)


###Making the plot for each module

##Turquoise

average_counts %>% 
  ggplot(aes(x=Tissue, y=MEturquoise, col = S_Sp, group = S_Sp)) +
  geom_hline(yintercept=0, color = "black", linetype="dashed")+
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  #facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("M5 N = 3156")

##NOTE: YOU CAN TRY FLIPPING AXES AND COLOUR COLUMNS TO SEE WHICH COMBO IS BEST TO VISUALIZE

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