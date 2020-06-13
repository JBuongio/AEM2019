getwd()
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
theme_set(theme_bw())
sharedfile = "Shared_all.shared" # What is in dissertation is the pick.shared file that is only Fe/S groups. Changing for publication. 
taxfile = "Sva.trimmed.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy"
mapfile = "Metadata.csv"
mothur_data <- import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile) # Import mothur data
map <-read.csv(mapfile) # Import sample metadata
head(map)
map <- sample_data(map)
rownames(map) <- map$Sample.Id # Assign rownames to be Sample ID's
moth_merge <- merge_phyloseq(mothur_data, map) # Merge mothurdata object with sample metadata
moth_merge
colnames(tax_table(moth_merge))
colnames(tax_table(moth_merge)) <-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
colnames(tax_table(moth_merge))
sample_sum_df <- data.frame(sum = sample_sums(moth_merge)) # Make a data frame with a column for the read counts of each sample
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
# mean, max and min of sample read counts
smin<-min(sample_sums(moth_merge))
smin #24
smax<-max(sample_sums(moth_merge))
smax #348261
smean<-mean(sample_sums(moth_merge))
smean #194122.6
sva_phylum<- moth_merge %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Phylum)
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)
ggplot(sva_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Svalbard sediments \n Bacterial Communities by Sampling Site") +
  theme(axis.text.x = element_text(angle =90, hjust = 1))


###### ordination with discrete scale, 60,000 min lib size (how data was normalized within Mothur)

# Source code files downloaded from ~/git_repos/MicrobeMiseq/R/miseqR.R
source("C:/Users/JoySpin/Documents/miseqR.R")

minlib = 60000
sva_scale<-scale_reads(moth_merge, minlib) #scale reads to even depth
sample_data(sva_scale)$Depth <- factor(
  sample_data(sva_scale)$Depth,
  levels = c(0.5,
             1.5,
             2.5,
             3.5,
             4.5,
             5.5,
             6.5,
             7.5,
             8.5,
             9.5,
             10.5,
             11.5,
             12.5,
             13.5,
             14.5,
             15.5,
             16.5,
             17.5,
             18.5,
             19.5)
)
require(devtools)

install_version("vegan", version ="2.4-5", repos = "http://cran.us.r-project.org") #Phyloseq qas written with dependency on an older Vegan package
#Restart R
library(vegan)

sva_pcoa<-ordinate(
  physeq = sva_scale,
  method = "PCoA",
  distance = "bray"
)
palette<-colfunc <- colorRampPalette(c("lightpink", "brown"))

plot_ordination(
  physeq = sva_scale,
  ordination = sva_pcoa,
  color = "Depth",
  shape = "Station",
  title = "PCoA Bray Curtis"
) +
  theme(text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(size=3) + scale_color_manual(values = palette(c(20)))

otu_table(sva_scale)

############### PCoA, reads > 60,000. continuous scale
sva_60000<-prune_samples(sample_sums(moth_merge)>60000, moth_merge)
sva_60000
otu_table(sva_60000)
sva_60000_pcoa<-ordinate(
  physeq = sva_scale,
  method = "PCoA",
  distance = "bray"
)

plot_ordination(
  physeq = sva_60000,
  ordination = sva_60000_pcoa,
  color = "Depth",
  shape = "Station",
  title = "PCoA Bray Curtis Stations AB and AC"
) +
  theme(text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(size=5)  + scale_color_gradient(low="lightpink", high = "brown") +
  geom_text(mapping = aes(label = Depth), size = 5, vjust = 1.5) 


write.csv(sva_60000_pcoa$vectors, file = "sva_60000_pcoa.csv")

########################### PCoA of individual sites
sva_AB<-moth_merge %>%
  subset_samples(Station=="AB")
sva_AC<-moth_merge %>%
  subset_samples(Station=="AC")

sva_AB_pcoa<-ordinate(
  physeq = sva_AB,
  method = "PCoA",
  distance = "bray"
)

plot_ordination(
  physeq = sva_AB,
  ordination = sva_AB_pcoa,
  color = "Depth",
  title = "PCoA of VK stn AB communities Bray Curtis"
) +
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(size = 4) + scale_colour_gradient(high = "brown", low = "lightpink")


sva_AC_pcoa<-ordinate(
  physeq = sva_AC,
  method = "PCoA",
  distance = "bray"
)

plot_ordination(
  physeq = sva_AC,
  ordination = sva_AC_pcoa,
  color = "Depth",
  title = "PCoA of VK stn AC communities Bray Curtis"
) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(size = 4) + scale_colour_gradient(high = "brown", low = "lightpink")


set.seed(1)
#################### NMDS plot updated with normalized libraries, 11/27/18 ###################################
library(grid)
sva_nmds<-ordinate(
  physeq = sva_60000, 
  method = "NMDS",
  distance = "bray"
)
#plot without vectors
NMDS_plot<-plot_ordination(
  physeq = sva_60000,
  ordination = sva_nmds,
  color = "Depth",
  shape = "Station",
  title = "NMDS of Svalbard bacterial Communities"
) +
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(size = 4) + scale_colour_gradient(high = "brown", low = "lightpink") +
  geom_text(mapping = aes(label = Depth), size = 5, vjust = 1.5) 

stressplot(sva_nmds)

geochem<-read.csv("Metadata_orgchem.csv") #read in geochemistry again to use in envfit
head(geochem)

samples<-get_sample(sva_60000) #get samples that made the 60000 read cut off
write.csv(samples, "samples_sva_60000.csv") #export to transpose into new file
good_samples<-read.csv("Good_samples.csv") #import to complete interjoin
good_samples_meta<-left_join(good_samples, geochem, by = "Sample.Id")
good_samples_meta
#using Vegan envfit function
(fit <-envfit(sva_nmds, good_samples_meta, perm =999))
fit_scores<-as.data.frame(scores(fit, display = "vectors"))
fit_scores<-cbind(fit_scores, env.variables = rownames(fit_scores))

#make distance matrix outside of the ones made embedded in ordiation calls
distance_matrix<-distance(sva_60000, method = "bray")
good_samples_meta_num<-good_samples_meta[, sapply(good_samples_meta, class) == "numeric"] #only want numbers for bioenv
sva_nmds_meta_bioenv<-bioenv(distance_matrix,good_samples_meta_num, use='p', metric = "manhattan")
summary(sva_nmds_meta_bioenv)
dbrda<-dbrda(distance_matrix ~ Depth + C.N +d13Corg, data= good_samples_meta)
anova(dbrda, by ='margin')
adonis(distance_matrix ~ Depth + C.N +d13Corg, data= good_samples_meta) #How significantly different are geochemistry measurements by station?
anisom_group<-get_variable(sva_60000, "Station")
anosim(distance_matrix,anisom_group)

#plot with vectors we just identified with envfit
NMDS_plot_vectors<-plot_ordination(
  physeq = sva_60000,
  ordination = sva_nmds,
  color = "Depth",
  shape = "Station",
  title = "NMDS of Svalbard bacterial Communities"
) +
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(size = 4) + scale_colour_gradient(high = "brown", low = "lightpink") +
  geom_text(mapping = aes(label = Depth), size = 5, vjust = 1.5) +
  geom_segment(data=fit_scores,aes(x=0,xend=NMDS1,y=0,yend=NMDS2, shape = NULL, color = NULL),
               color = "black", size = 1, linetype = 1, arrow = arrow(length = unit(0.1,"cm"))) + 
  geom_text(data = fit_scores, 
            aes(x = NMDS1, y = NMDS2, label=env.variables, shape = NULL, color = NULL),
            size = 6,
            hjust = .75)


########### Doing new distance matrix only on libraries with more than 60,000 reads for CAP ordination (10/19/2018) ###################
VK_60000_not_na<-sva_60000 %>%
  subset_samples(
    !is.na(Hydrogen) &
      !is.na(d13Corg) & 
      !is.na("%C") & 
      !is.na("C/N")
  )

head(sample_data(sva_60000))
colnames(sample_data(VK_60000_not_na)) <-c("Sample Id",	"Fjord",	"Station", "Replicate", "Depth",	"Hydrogen",	"d13Corg", "Percent_C",	"CtoN", "Fe", "Mn")
VK_bray_60000_not_na<-phyloseq::distance(VK_60000_not_na, method = "bray")
sampledf_VK_60000<-data.frame(sample_data(sva_60000))  
adonis(VK_bray_60000_not_na ~ Station, data= sampledf_VK_60000)  
beta_VK_60000<-betadisper(VK_bray_60000_not_na, sampledf_VK_60000$Station)  
permutest(beta_VK_60000)
sva_60000_bray<-phyloseq::distance(physeq=sva_60000, method="bray")

############ CAP ord plot, remade with "proximity to glacier" removed from metadata 11/26/18 ###################
cap_ord_VK <- ordinate(
  physeq = VK_60000_not_na, 
  method = "CAP",
  distance = VK_bray_60000_not_na,
  formula = ~ Hydrogen + d13Corg + Percent_C + CtoN + Depth + Fe + Mn)

cap_plot <- plot_ordination(
  physeq = VK_60000_not_na, 
  ordination = cap_ord_VK, 
  color = "Depth",
  axes = c(1,2)
) + 
  aes(shape = Station) + 
  theme(text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(size = 4) + scale_colour_gradient(high = "brown", low = "lightpink") +
  geom_text(mapping = aes(label = Depth), size = 5, vjust = 1.5) 
arrowmat <- vegan::scores(cap_ord_VK, display = "bp")
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)
label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)
arrowhead = arrow(length = unit(0.02, "npc")) 

cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .75, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 5,
    data = arrowdf,
    color = "dodgerblue3",
    show.legend = FALSE
  )

anova(cap_ord_VK, by="terms", perm.max=500)
anova(cap_ord_VK, by = "margin")
anova(cap_ord_VK)


sample_dist<-vegdist(tax_table(sva_60000), method = "bray")
mantel(sva_bray, bray_not_na)


########## diversity

pal="Set1"
moth_merge
plot_richness(sva_60000)
plot_richness(sva_60000, measures = c("Chao1", "Shannon"))
plot_richness(sva_60000, x="Depth", measures = c("Chao1", "Shannon"))
sample_data(sva_60000)$fjord<-get_variable(sva_60000, "Station") %in% c("AB", "AC")
plot_richness(sva_60000, x="Depth", color="Station", measures = c("Chao1", "Shannon"))
sample_data(sva_60000)$fjord<-get_variable(sva_AB_AC, "Station") %in% c("AB", "AC")
plot_richness(sva_60000, x="Depth", color="Station", measures = c("Chao1","Shannon"))
number_ticks<-function(n) {function(limits) pretty (limits, n)}
p<-plot_richness(sva_60000, x="Depth", color="Station", measures = c("Chao1", "Shannon", "Simpson")) + geom_line() +
  scale_x_reverse() +
  coord_flip()
plots<-layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
       widths=c(3,1), heights=c(1,2))


#### plots used in Appendix
d<-plot_richness(sva_60000, x="Depth", color="Station", measures = c("Shannon")) + 
  geom_line() +
  geom_point(aes(size=7)) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

K<-plot_richness(sva_60000, x="Depth", color="Station", measures = c("Chao1")) + 
  geom_line() +
  geom_point(aes(size=7)) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

L<-plot_richness(sva_60000, x="Depth", color="Station", measures = c("Simpson")) + 
  geom_line() +
  geom_point(aes(size=7)) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
    coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

richness_AB_AC<-estimate_richness(sva_60000)
write.csv(richness_AB_AC, "Richness_AB_AC_60000.csv")


