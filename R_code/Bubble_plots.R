#Joy Buongiorno
#jbuongior21@gmail.com
#make bubble plot figure for paper
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)

#norm_AB<-read.csv("AB_norm.csv")
#head(norm_AB)
#norm_AB_tidy<-gather(norm_AB, depth, abundance, X0.5:X18.5)
#head(norm_AB_tidy)
#write.csv(norm_AB_tidy, "norm_AB_tidy.csv")

bub_cols<-c("-9" = "#6346f2", "1"="#01b64e","-10" = "#01b64e", "2"= "#6346f2", "-8" = "#cbc600", "3"="#cbc600", "-7" ="#cb0049", "4"="#cb0049", "-6" = "#3aedc7", "5"="#3aedc7", "6"="#ffb89c")

############### AB bubble plot #############
norm_AB_tidy_fixed<-read.csv("norm_AB_tidy_fixed.csv", stringsAsFactors = FALSE)
head(norm_AB_tidy_fixed)
norm_AB_tidy_fixed$Genus <- reorder(norm_AB_tidy_fixed$Genus,norm_AB_tidy_fixed$guild)
norm_AB_tidy_fixed$Percent_Abundance<-as.numeric(as.character(norm_AB_tidy_fixed$Percent_Abundance))
norm_AB_tidy_fixed$Genus <- reorder(norm_AB_tidy_fixed$Genus,norm_AB_tidy_fixed$Percent_Abundance) #don't do this if you want the sums separated
norm_AB_tidy_fixed$guild <- as.factor(norm_AB_tidy_fixed$guild)


ggplot(norm_AB_tidy_fixed, aes(x=depth, y=Genus, color=guild, size=Percent_Abundance)) +
  geom_point(aes(fill=guild)) +
  scale_color_manual(values=bub_cols) +
  scale_size(range = c(0,10), breaks = c(0.001, 0.01, 0.1, 1, 1.5, 2, 4)) +
  coord_flip() +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Depth (cmbsf)", y= "Genus") +
  scale_x_reverse(limits= c(20,0))


##############################################
#norm_AC<-read.csv("AC_norm.csv")
#head(norm_AC)
#norm_AC_tidy<-gather(norm_AC, depth, abundance, X0.5:X19.5)
#head(norm_AC_tidy)
#write.csv(norm_AC_tidy, "norm_AC_tidy.csv")

###################AC Bubble plot ############
norm_AC_tidy_fixed<-read.csv("norm_AC_tidy_fixed.csv", stringsAsFactors = FALSE)
head(norm_AC_tidy_fixed)
norm_AC_tidy_fixed$Genus <- reorder(norm_AC_tidy_fixed$Genus,norm_AC_tidy_fixed$guild)
norm_AC_tidy_fixed$Percent_Abundance<-as.numeric(as.character(norm_AC_tidy_fixed$Percent_Abundance))
norm_AC_tidy_fixed$Genus <- reorder(norm_AC_tidy_fixed$Genus,norm_AC_tidy_fixed$Percent_Abundance)
norm_AC_tidy_fixed$guild <- as.factor(norm_AC_tidy_fixed$guild)


ggplot(norm_AC_tidy_fixed, aes(x=depth, y=Genus, color=guild, size=Percent_Abundance)) +
  geom_point(aes(fill=guild)) +
  scale_color_manual(values=bub_cols) +
  scale_size(range = c(0,10), breaks = c(0.001, 0.01, 0.1, 1, 1.5, 2, 4)) +
  coord_flip() +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Depth (cmbsf)", y= "Genus") +
  scale_x_reverse(limits= c(20,0)) 
