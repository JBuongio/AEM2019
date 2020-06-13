#Joy Buongiorno
#jbuongior21@gmail.com
#stacked bar figure

####Family level, collected with get.lineage command in mothur####
library(tidyr)
library(dplyr)
library(ggplot2)
AB_FeS<-read.csv("AB_Fe_S_fams_new_seqs.csv")
AB_FeS_tidy<-gather(AB_FeS, depth, abundance, AB_1.5:AB_14.5)
##### don't rewrite this ####
write.csv(AB_FeS_tidy, "AB_FeS_tidy.csv")
AB_FeS_sum <- AB_FeS_tidy %>%
  group_by(depth) %>%
  summarise(sum(abundance))
write.csv(AB_FeS_sum, file="AB_FeS_tidy_sum.csv")
AB_FeS_tidy_clean<-read.csv("AB_FeS_tidy.csv")

library(RColorBrewer)
AB_FeS_tidy_clean$Family_name <- reorder(AB_FeS_tidy_clean$Family_name,AB_FeS_tidy_clean$abundance)
AB_FeS_tidy_clean$Family_name <- factor(AB_FeS_tidy_clean$Family_name, levels=rev(levels(AB_FeS_tidy_clean$Family_name)))
length(unique(AB_FeS_tidy_clean$Family_name))

AB_FeS_plot<-ggplot(AB_FeS_tidy_clean, aes(x = depth, y = abundance, fill = Family_name)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1))+ 
  xlab(AB_FeS_tidy_clean$depth) +
  labs(x="Sample", y="Abundance (%)") +
  scale_x_discrete(limits=rev(levels(AB_FeS_tidy_clean$depth))) +
  coord_flip() +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=colorRampPalette(c("gray", "maroon4", "wheat3", "peachpuff4", "lightblue1", "lemonchiffon3", "midnightblue","hotpink","mintcream", "purple", "green", "royalblue", "orange", "turquoise2" , "red", "purple", "salmon2", "blue"))( 44))

AC_FeS<-read.csv("AC_Fe_S_fams_new_seqs.csv")
AC_FeS_tidy<-gather(AC_FeS, depth, abundance, AC_0.5:AC_16.5)
write.csv(AC_FeS_tidy, "AC_FeS_tidy.csv")
AC_FeS_sum<-AC_FeS_tidy %>%
  group_by(depth) %>%
  summarise(sum(abundance))
write.csv(AC_FeS_sum, file="AC_FeS_tidy_sum.csv")
AC_FeS_tidy_clean<-read.csv("AC_FeS_tidy.csv")
AC_FeS_tidy_clean$Family_name <- reorder(AC_FeS_tidy_clean$Family_name,AC_FeS_tidy_clean$abundance)
AC_FeS_tidy_clean$Family_name <- factor(AC_FeS_tidy_clean$Family_name, levels=rev(levels(AC_FeS_tidy_clean$Family_name)))
length(unique(AC_FeS_tidy_clean$Family_name))

AC_FeS_plot<-ggplot(AC_FeS_tidy_clean, aes(x = depth, y = abundance, fill = Family_name)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1))+ 
  xlab(AC_FeS_tidy_clean$depth) +
  labs(x="Sample", y="Abundance (%)") +
  scale_x_discrete(limits=rev(levels(AC_FeS_tidy_clean$depth))) +
  coord_flip() +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=colorRampPalette(c("gray", "maroon4", "wheat3", "peachpuff4", "lightblue1", "lemonchiffon3", "midnightblue","hotpink","mintcream", "purple", "green", "royalblue", "orange", "turquoise2" , "red", "purple", "salmon2", "blue"))( 42))

