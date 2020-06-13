#Joy Buongiorno
#jbuongior21@gmail.com
#stacked bar plot in supplementals

####Family level
library(tidyr)
library(dplyr)
library(ggplot2)
VK<-read.csv("AB_AC_lvl5.csv")
VK_tidy<-gather(VK, depth, abundance, AB_1.5:AC_16.5)
write.csv(VK_tidy, "VK_tidy.csv")
VK_tidy_top<-VK_tidy %>%
  group_by(depth) %>%
  top_n(n = 30, abundance)
VK_tidy_top_sum <- VK_tidy_top %>%
  group_by(depth) %>%
  summarise(sum(abundance))
write.csv(VK_tidy_top, file="VK_tidy_top.csv")
write.csv(VK_tidy_top_sum, file="VK_tidy_top_sum.csv")

VK_tidy_top<-read.csv("VK_tidy_top.csv")
VK_tidy_top_sum <- VK_tidy_top %>%
  group_by(depth) %>%
  summarise(sum(abundance))
library(RColorBrewer)
VK_tidy_top$Family_name <- reorder(VK_tidy_top$Family_name, VK_tidy_top$abundance)
VK_tidy_top$Family_name <- factor(VK_tidy_top$Family_name, levels=rev(levels(VK_tidy_top$Family_name)))
length(unique(VK_tidy_top$Family_name))
VK_plot<-ggplot(VK_tidy_top, aes(x = sample, y = abundance, fill = Family_name)) + theme(axis.text.x=element_text(angle=90,hjust=1))+ xlab(VK_tidy_top$sample) + geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values=colorRampPalette(c("gray", "maroon4", "wheat3", "peachpuff4", "lightblue1", "lemonchiffon3", "midnightblue","hotpink","mintcream", "purple", "green", "royalblue", "orange", "turquoise2" , "red", "purple", "salmon2", "blue"))( 66))
