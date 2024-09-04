##############################################################################################################
#######REDRAWING FIGURES FROM GPROFILER FOR SPECIES-BIASED GENES##############################################
##############################################################################################################

#####bidirectional barplot for enrichment from gprofiler
library(ggplot2)

grprofiler_1 = read.table("../New annotation/Results/Plots/Part 1/gprofiler_species_biasedgenes.txt", sep="\t", head=T, comment.char="#")
View(grprofiler_1)

###This saved txt file was obtained after running grpofiler on species-biased genes with c.remanei gene (old) gene names

grprofiler_1 = grprofiler_1 %>% filter(grprofiler_1$Type != "Underrepresentation")

library(forcats)#to easily reorder varibales 
ggplot(grprofiler_1,aes(x=fct_reorder(grprofiler_1$term_name, grprofiler_1$negative_log10_of_adjusted_p_value), y = intersection_size, fill=negative_log10_of_adjusted_p_value))+ 
  geom_bar(stat="identity", position="identity")+
  scale_fill_gradient(low="#cc9966", high="#663300")+
  scale_y_continuous(limits = c(0, max(grprofiler_1$intersection_size))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+
  coord_flip()

######################ONLY MF##############3
grprofiler_mf = read.table("../New annotation/Results/Plots/Part 1/MF_Gprofiler_species_biasedenrichment.txt", sep="\t", head=T, comment.char="#")
View(grprofiler_mf)

###This saved txt file was obtained after running grpofiler on species-biased genes with c.remanei gene (old) gene names


ggplot(grprofiler_mf,aes(x=fct_reorder(grprofiler_mf$term_name, grprofiler_mf$negative_log10_of_adjusted_p_value), y = intersection_size, fill=negative_log10_of_adjusted_p_value))+ 
  geom_bar(stat="identity", position="identity")+
  scale_fill_gradient(low="#cc9966", high="#663300")+
  scale_y_continuous(limits = c(0, max(grprofiler_mf$intersection_size))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+
  coord_flip()

######################ONLY BP##############3
grprofiler_bp = read.table("../New annotation/Results/Plots/Part 1/BP_Gprofiler_species_biasedenrichment.txt", sep="\t", head=T, comment.char="#")
colnames(grprofiler_bp) = colnames(grprofiler_mf)

View(grprofiler_bp)

###This saved txt file was obtained after running grpofiler on species-biased genes with c.remanei gene (old) gene names


ggplot(grprofiler_bp,aes(x=fct_reorder(grprofiler_bp$term_name, grprofiler_bp$negative_log10_of_adjusted_p_value), y = intersection_size, fill=negative_log10_of_adjusted_p_value))+ 
  geom_bar(stat="identity", position="identity")+
  scale_fill_gradient(low="#cc9966", high="#663300")+
  scale_y_continuous(limits = c(0, max(grprofiler_bp$intersection_size))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+
  coord_flip()

#LOOKING AT EACH MODULE SPECIFICALLY

##Linear model fitting on modules to know how they are effected 

library(lme4)
library(tidyverse)
library(dplyr)

model.matrix(~coldata2$sex) # looking at model matrix 
contrasts(species.factor) #lookign at matrix for factorized metadata

View(coldata2)
View(total_data)

total_data = bind_cols(module_eigengenes, coldata2)

#####average_counts is the file 
average_counts = read.table("WGCNA/Average_ME_across_reps.txt", sep="\t", head=T, row.name=1, comment.char = "#")

total_data %>% filter(species == "Cre") %>% 
  ggplot(aes(tissue, MEmagenta, col = sex)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)


total_data %>% 
  ggplot(aes(sex, MEmagenta, col = tissue)) +
  geom_point() +
  facet_wrap(~species)


#write.csv(total_data, file = "WGCNA/total_data_with ME.csv", sep = "\t",
#          row.names = TRUE, col.names = NA)

# When building linear model, there are different ways to encode categorical variables, 
# known as contrast coding systems. The default option in R is to use the first level of 
# the factor as a reference and interpret the remaining levels relative to this level

plot(module_eigengenes)


average_counts = read.table("WGCNA/Average_ME_across_reps.txt", sep="\t", head=T, comment.char="#")


View(average_counts)



#Calculating the effects of different factors for different modules
#MEMagenta
anova(lm(module_eigengenes$MEmagenta ~ coldata2$species + coldata2$sex + coldata2$tissue + coldata2$batch + coldata2$sex*coldata2$tissue +
           coldata2$species*coldata2$tissue + coldata2$sex*coldata2$species + coldata2$sex*coldata2$tissue*coldata2$species)) # -1 is to remove intercept and include the reference level

summary(lm(module_eigengenes$MEmagenta ~ coldata2$species + coldata2$sex + coldata2$tissue + coldata2$batch + coldata2$sex*coldata2$tissue +
             coldata2$species*coldata2$tissue + coldata2$sex*coldata2$species + coldata2$sex*coldata2$tissue*coldata2$species  )) # -1 is to remove intercept and include the reference level

magenta = average_counts %>% filter(average_counts$Species != "H1", 
                                    average_counts$Sex != "W", average_counts$Tissue != "W") %>% 
  ggplot(aes(Sex, MEmagenta, col = Tissue, group = Tissue)) +
  geom_point(show.legend=T, alpha=1, size = 2.5) +
  theme_bw() +
  facet_wrap(~Species) +
  stat_smooth(method = "lm") +
  ggtitle("MEMagenta N = 273")


file$module_name = "cyan"
