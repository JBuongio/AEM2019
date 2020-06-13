#### Fe/S after second MiSeq run ###
# 1. Make file of Fe_S only famliies: awk {'print $3'} Sva_All.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.pick.taxonomy | cut -d ";" -f 5-6 | cut -d "(" -f 1 > Fe_S_fams.txt
# 2. Remove "Unknown_Families": grep -v "Unknown" Fe_S_fams.txt > Fe_S_fams_no_unknown.txt
# 3. Make file of only families from each station in Excel, and calculate relative abundance for each library (or make mothur calc it for you in a summary file or normalized.shared file instead of doing it in Excel).
# 4. Grep Fe_S families from rel abund file:  grep -f Fe_S_fams_no_unknown.txt AB1_all_families.txt > AB1_Fe_S_fams.txt
# 5. Add headers that will allow proper sorting (ex. AB1_00.5).
# 6. Follow directions below, Thinking to combine all data together so taht colours remain uniform .
# 7a. Remove families in low abundance after getting sum files so that No_FeS math doesn't get messed up with removal of low abundance families. Cut off is <1 % summed across all libraries. that do not turn up in more than one library (save in sep file "..._too_low.txt"). do this by going back to the original file you read in, sum across libraries, sort, and grab the names of the families <1%.
# 7b. Put grabbed family names from each replicate into seperate files, but then get shared names between each file with: comm -12 <(sort AB2_Fe_S_fams_too_low.txt) <(sort AB1_Fe_S_fams_too_low.txt) > AB12_Fe_S_fams_too_low.txt
# 7c. Grep out these families: grep -vf AB12_Fe_S_fams_too_low.txt AB12_Fe_S_tidy_fixed.csv > AB12_Fe_S_tidy_above_1percent.txt. Save that as csv after opening in excel. 

# This code provides relabund of Fe/S taxa only, all adding to 100%.

library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)

fam_cols<-c("Other"= "white", "Woeseiaceae/JTB255" = "#7B4F4B", "Desulfobacteraceae"	= "#0CBD66" ,"Nitrosomonadales_unclassified" = "#636375", "Desulfobulbaceae" = "purple4","Sva1033"	= "#FF34FF", "Bacteriovoracaceae" =	"#FF4A46","BVA18" =	"#008941", "Acidimicrobiaceae" =	"#006FA6","Desulfovibrionaceae"	= "#A30059", "Piscirickettsiaceae" =	"#FFDBE5","Rhodobacteraceae" = "#7A4900","Helicobacteraceae" =	"#0000A6",   "Alteromonadaceae" =	"#63FFAC","Desulfarculaceae" =	"#B79762","GR-WP33-58" =	"#004D43","Desulfuromonadales_unclassified" =	"#8FB0FF","M113"	= "#997D87","M20-Pitesti" = "#5A0007", "Geobacteraceae"  =	"#809693", "Desulfuromonadaceae"  =	"#FEFFE6", "C8S-102" =	"#1B4400","Nitrosomonadaceae" =	"#4FC601","Thermoplasmatales_Incertae_Sedis" =	"#3B5DFF", "Campylobacteraceae"	= "#4A3B53", "Gallionellaceae" =	"#FF2F80", "Leptospiraceae"	= "#61615A", "AMOS1A-4113-D04"	= "#BA0900", "CCA47"	= "#6B7900", "VC2.1_Arc6"	= "#00C2A0", "MKCST-A3" =	"#FFAA92","Nitrosomondales_unclassified"	= "#FF90C9", "Hydrogenophilaceae" =	"#B903AA","Mariprofundaceae"	= "#D16100", "Thiotrichaceae" =	"#DDEFFF","Ectothiorhodospiraceae" =	"#c3e840", "Thermoplasmatales_unclassified"	= "#A1C299","Syntrophaceae" =	"#300018", "Comamonadaceae" =	"#0AA6D8","Xanthomonadaceae" =	"#013349", "ANT06-05" =	"#00846F","Marine_Group_II" =	"#f45042", "Acidithiobacillaceae" =	"#FFB500", "Marine_Group_III"	= "#C2FFED", "Thiotrichales_Incertae_Sedis" =	"#A079BF", "Marine_Benthic_Group_D_and_DHVEG-1"	= "#CC0744", "Campylobacterales_unclassified" =	"#3A2465"
)
           
              
 #AB1
AB1_FeS<-read.csv("AB1_Fe_S_fams.csv")
head(AB1_FeS)
AB1_FeS_tidy<-gather(AB1_FeS, depth, abundance, AB1_00.5:AB1_14.5)
head(AB1_FeS_tidy)
write.csv(AB1_FeS_tidy, "AB1_Fe_S_tidy.csv")
AB1_Fe_S_sum<- aggregate(AB1_FeS_tidy$abundance, by=list(Category=AB1_FeS_tidy$depth), FUN=sum)
write.csv(AB1_Fe_S_sum, "AB1_Fe_S_sum.csv") #use this to figureo out "other" and add that number to the tidy file (unless it is really small)
AB1_FeS_tidy_fixed<-read.csv("AB1_Fe_S_tidy_fixed.csv")

library(RColorBrewer)
AB1_FeS_tidy_fixed$Family_name <- reorder(AB1_FeS_tidy_fixed$Family_name,AB1_FeS_tidy_fixed$abundance)
AB1_FeS_tidy_fixed$Family_name <- factor(AB1_FeS_tidy_fixed$Family_name, levels=rev(levels(AB1_FeS_tidy_fixed$Family_name)))

length(unique(AB1_FeS_tidy_fixed$Family_name))

AB1_FeS_plot<-ggplot(AB1_FeS_tidy_fixed, aes(x = depth, y = abundance, fill = Family_name)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1))+ 
  xlab(AB1_FeS_tidy_fixed$depth) +
  labs(x="Sample", y="Abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(limits=rev(levels(AB1_FeS_tidy_fixed$depth))) +
  geom_bar(stat = "identity") + 
  coord_flip(ylim=c(0,20)) +
  scale_fill_manual(values=fam_cols)

#AB2
AB2_FeS<-read.csv("AB2_Fe_S_fams.csv")
head(AB2_FeS)
AB2_FeS_tidy<-gather(AB2_FeS, depth, abundance, AB2_00.5:AB2_18.5)
head(AB2_FeS_tidy)
write.csv(AB2_FeS_tidy, "AB2_Fe_S_tidy.csv")
AB2_Fe_S_sum<-aggregate(AB2_FeS_tidy$abundance, by=list(Category=AB2_FeS_tidy$depth), FUN=sum)
write.csv(AB2_Fe_S_sum, "AB2_Fe_S_sum.csv")
head(AB2_Fe_S_sum)
AB2_FeS_tidy_fixed<-read.csv("AB2_Fe_S_tidy_fixed.csv")

AB2_FeS_tidy_fixed$Family_name <- reorder(AB2_FeS_tidy_fixed$Family_name,AB2_FeS_tidy_fixed$abundance)
AB2_FeS_tidy_fixed$Family_name <- factor(AB2_FeS_tidy_fixed$Family_name, levels=rev(levels(AB2_FeS_tidy_fixed$Family_name)))

length(unique(AB2_FeS_tidy_fixed$Family_name))

AB2_FeS_plot<-ggplot(AB2_FeS_tidy_fixed, aes(x = depth, y = abundance, fill = Family_name)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1))+ 
  xlab(AB2_FeS_tidy_fixed$depth) +
  labs(x="Sample", y="Abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(limits=rev(levels(AB2_FeS_tidy_fixed$depth))) +
  geom_bar(stat = "identity") + 
  coord_flip(ylim=c(0,20)) +
  scale_fill_manual(values=fam_cols)

#AC1

AC1_FeS<-read.csv("AC1_Fe_S_fams.csv")
head(AC1_FeS)
AC1_FeS_tidy<-gather(AC1_FeS, depth, abundance, AC1_00.5:AC1_18.5)
head(AC1_FeS_tidy)
write.csv(AC1_FeS_tidy, "AC1_Fe_S_tidy.csv")
AC1_Fe_S_sum<-aggregate(AC1_FeS_tidy$abundance, by=list(Category=AC1_FeS_tidy$depth), FUN=sum)
write.csv(AC1_Fe_S_sum, "AC1_Fe_S_sum.csv")
head(AC1_Fe_S_sum)
AC1_FeS_tidy_fixed<-read.csv("AC1_Fe_S_tidy_fixed.csv")


AC1_FeS_tidy_fixed$Family_name <- reorder(AC1_FeS_tidy_fixed$Family_name,AC1_FeS_tidy_fixed$abundance)
AC1_FeS_tidy_fixed$Family_name <- factor(AC1_FeS_tidy_fixed$Family_name, levels=rev(levels(AC1_FeS_tidy_fixed$Family_name)))

length(unique(AC1_FeS_tidy_fixed$Family_name))

AC1_FeS_plot<-ggplot(AC1_FeS_tidy_fixed, aes(x = depth, y = abundance, fill = Family_name)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1))+ 
  xlab(AC1_FeS_tidy_fixed$depth) +
  labs(x="Sample", y="Abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(limits=rev(levels(AC1_FeS_tidy_fixed$depth))) +
  geom_bar(stat = "identity") + 
  coord_flip(ylim=c(0,20)) +
  scale_fill_manual(values=fam_cols)

#AC2
AC2_FeS<-read.csv("AC2_Fe_S_fams.csv")
head(AC2_FeS)
AC2_FeS_tidy<-gather(AC2_FeS, depth, abundance, AC2_01.5:AC2_19.5)
head(AC2_FeS_tidy)
write.csv(AC2_FeS_tidy, "AC2_Fe_S_tidy.csv")
AC2_Fe_S_sum<-aggregate(AC2_FeS_tidy$abundance, by=list(Category=AC2_FeS_tidy$depth), FUN=sum)
write.csv(AC2_Fe_S_sum, "AC2_Fe_S_sum.csv")
head(AC2_Fe_S_sum)
AC2_FeS_tidy_fixed<-read.csv("AC2_Fe_S_tidy_fixed.csv")

AC2_FeS_tidy_fixed$Family_name <- reorder(AC2_FeS_tidy_fixed$Family_name,AC2_FeS_tidy_fixed$abundance)
AC2_FeS_tidy_fixed$Family_name <- factor(AC2_FeS_tidy_fixed$Family_name, levels=rev(levels(AC2_FeS_tidy_fixed$Family_name)))

length(unique(AC2_FeS_tidy_fixed$Family_name))

AC2_FeS_plot<-ggplot(AC2_FeS_tidy_fixed, aes(x = depth, y = abundance, fill = Family_name)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1))+ 
  xlab(AC2_FeS_tidy_fixed$depth) +
  labs(x="Sample", y="Abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(limits=rev(levels(AC2_FeS_tidy_fixed$depth))) +
  geom_bar(stat = "identity") + 
  coord_flip(ylim=c(0,20)) +
  scale_fill_manual(values=fam_cols)

