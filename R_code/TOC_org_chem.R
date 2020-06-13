#Joy Buongiorno
#jbuongior21@gmail.com
#Organic geochem plots

library(ggplot2)
install.packages("ggThemeAssist")
library(ggThemeAssist)
library(plyr)
org_chem<-read.csv("TOC.csv")
#####################################################################
TOC_P<-org_chem[1:11,4]
depth_P<-org_chem[1:11,3]
TOC_Pbind<-cbind(TOC_P, depth_P)
TOC_F<-org_chem[12:26, 4]
depth_F<-org_chem[12:26,3]
TOC_Fbind<-cbind(TOC_F,depth_F)

################## Plots for all TOC data across KF and VK together ################
TOC_plot<-ggplot(org_chem, aes(x=Depth, y=TOC, shape=Station, color=Fjord)) + 
  geom_point(aes(size=7)) +
  theme(text = element_text(size=20)) +
  geom_line() +
  scale_x_reverse() +
  coord_flip()
##########################CtoN vs d13C plots for all fjords ###################
CtoN_plot<-ggplot(org_chem, aes(x=d13C, y=CtoN, shape=Station, color=Fjord)) +
  geom_point(aes(size=7)) +  
  theme(text = element_text(size=20)) +
  coord_flip()
  
##########################organic isotopes, KF ####################
KF<-org_chem[1:26,]
KF_plot<-ggplot(KF, aes(x=Depth, y=d13C, shape=Station, color=Station)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(23, 22)) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(20, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  coord_flip() +
labs(x="Depth", y=expression(paste(delta^{13}, "C (\u2030 vs. PDB)")))
##################### TOC, KF ###################3
KF_TOC_plot<-ggplot(KF, aes(x=Depth, y=TOC, shape=Station, color=Station)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(23, 22)) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(20, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Total Organic Carbon (TOC)") +
  coord_flip()


################C/N, KF ###################

KF_cton_plot<-ggplot(KF, aes(x=Depth, y=CtoN, shape=Station, color=Station)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(23, 22)) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(20, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "C/N") +
  coord_flip()

########## Isolate geochem from each station ####################
VK<-org_chem[27:84,]
VKAB<-VK[which(VK$Station == "AB"), names(VK) %in% c("Fjord", "Station", "Depth", "TOC", "CtoN", "d13C", "d15N")]
VKAC<-org_chem[44:63,]
VKHA<-VK[which(VK$Station == "HA"), names(VK) %in% c("Fjord", "Station", "Depth", "TOC", "CtoN", "d13C", "d15N")]
################ Cook's distance to determine outliers, change data set each time to create plots ################
mod<-lm(CtoN ~ Depth , data=VKHA)
cooksd<-cooks.distance(mod)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance, CtoN in HA")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels)
influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])  # influential row numbers from original data frame

##########t tests for data without outliers ################
org_chem_no_out<-read.csv("TOC_outliers_removed.csv")
VK_no_out<-org_chem_no_out[27:84,]
VKAB_no_out<-VK_no_out[which(VK_no_out$Station == "AB"), names(VK_no_out) %in% c("Fjord", "Station", "Depth", "TOC", "CtoN", "d13C", "d15N")]
VKAC_no_out<-org_chem_no_out[44:63,]
VKHA_no_out<-VK_no_out[which(VK_no_out$Station == "HA"), names(VK_no_out) %in% c("Fjord", "Station", "Depth", "TOC", "CtoN", "d13C", "d15N")]
shapiro.test(VKAB_no_out$TOC)
shapiro.test(VKAC_no_out$TOC) #All arrays are normally distributed.
shapiro.test(VKHA_no_out$TOC)

shapiro.test(VKAB_no_out$CtoN)
shapiro.test(VKAC_no_out$CtoN)
shapiro.test(VKHA_no_out$CtoN)

shapiro.test(VKAB_no_out$d13C)
shapiro.test(VKAC_no_out$d13C) #sig diff from normal
shapiro.test(VKHA_no_out$d13C)


t.test(VKAB_no_out$TOC,VKAC_no_out$TOC) #Welsh two sample t test -  0.004435
t.test(VKAB_no_out$TOC, VKHA_no_out$TOC) #Welsh two sample t test - 0.0002843
t.test(VKHA_no_out$TOC, VKAC_no_out$TOC) #Welsh two sample t test - 0.228

t.test(VKAB_no_out$d13C,VKAC_no_out$d13C) #Welsh two sample t test - 0.005297
t.test(VKAB_no_out$d13C, VKHA_no_out$d13C) #Welsh two sample t test - 0.0006274
t.test(VKHA_no_out$d13C, VKAC_no_out$d13C) #Welsh two sample t test - 0.06565

t.test(VKAB_no_out$CtoN,VKAC_no_out$CtoN) #Welsh two sample t test - 0.05487
t.test(VKAB_no_out$CtoN, VKHA_no_out$CtoN) #Welsh two sample t test - 0.005338
t.test(VKHA_no_out$CtoN, VKAC_no_out$CtoN) #Welsh two sample t test - 0.558

##############crossplots with TOC for supp####

VK_toc_d13C<-ggplot(VK_no_out, aes(x=TOC, y=d13C, shape=Station, color=Station)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_colour_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  geom_smooth(method=lm, se=FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="TOC (wt%)", y=expression(paste(delta^{13}, "C (\u2030 vs. PDB)")))

VK_CtoN2<-ggplot(VK_no_out, aes(x=TOC, y=CtoN, shape=Station, color=Station)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_colour_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  geom_smooth(method=lm, se=FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="TOC wt(%)", y="C/N")

### for defense talk ###
VK_CtoN<-ggplot(VK_no_out, aes(x=CtoN, y=d13C, shape=Station, color=Station)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_colour_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 30) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="C/N", y=expression(paste(delta^{13}, "C (\u2030 vs. PDB)")))

#######################################################################


VK_AC_plot<-ggplot(VK_AC, aes(x=Depth, y=d13C, shape=Station, color=Station)) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24)) +
  theme(text = element_text(size=20)) +
  scale_y_continuous(limits=c(-27.5, -25.5)) +
  scale_x_reverse() +
  geom_line(aes(color=Station)) +
  geom_line(size=1.5) +
  labs(x="Depth (cmbsf)", y=expression(paste(delta^{13}, "C (\u2030 vs. PDB)"))) +
  coord_flip()
############### downcore plots for fig 2, outliers removed ##############
VK_cton_plot<-ggplot(VK_no_out, aes(x=Depth, y=CtoN, shape=Station, color=Station)) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Station), size=0.75) +
  labs(x="Depth (cmbsf)", y="CtoN") +
  coord_flip()

VK_d13C_plot<-ggplot(VK_no_out, aes(x=Depth, y=d13C, shape=Station, color=Station)) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits=c(-28, -24)) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Station), size=0.75) +
  labs(x="Depth (cmbsf)", y=expression(paste(delta^{13}, "C (\u2030 vs. PDB)"))) +
  coord_flip()

VK_TOC<-ggplot(VK_no_out, aes(x=Depth, y=TOC, shape=Station, color=Station)) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_x_reverse() +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_line(aes(color=Station, linetype=Station), size=0.75) +
  labs(x="Depth (cmbsf)", y= "TOC (wt%)") +
  coord_flip()

############################################################################

VK_CtoN<-ggplot(VK, aes(x=CtoN, y=d13C, shape=Station, color=Station)) +
  scale_color_manual(values=c("turquoise3", "turquoise4")) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24)) +
  theme(text = element_text(size=20)) +
  scale_x_continuous(limits = c(0,80)) +
  scale_y_continuous(limits = c(-35,-5)) +
  labs(y=expression(paste(delta^{13}, "C (\u2030 vs. PDB)")))



VK_CtoN_zoomedin<-ggplot(VK, aes(x=CtoN, y=d13C, shape=Station, color=Station)) +
  scale_color_manual(values=c("turquoise3", "turquoise4")) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24)) +
  theme(text = element_text(size=20)) +
  labs(y=expression(paste(delta^{13}, "C (\u2030 vs. PDB)")))

library(scales)
show_col(hue_pal()(9))
cols<-c("AB.2" = "#F8766D","AB.1" = "#F8766D", "AC.1" = "#C77CFF", "AC.2" = "#C77CFF", "AB" = "#F8766D", "AC" = "#C77CFF", "HA" = "#00B9E3")


VKAC_TOC<-ggplot(VK_AC, aes(x=Depth, y=TOC, shape=Station, color=Station)) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24)) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  geom_line(aes(color=Station)) +
  geom_line(size=1.5) +
  labs(x="Depth (cmbsf)", y= "TOC (wt%)") +
  coord_flip()


VK_H<-read.csv("Hydrogen.csv")
VK_H_plot<-ggplot(VK_H, aes(x=Depth, y=Hydrogen, shape=Rep, color=Station)) +
  geom_point(size=4) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  geom_line(aes(color=Station)) +
  geom_line(size=1.5) +
  labs(x="Depth (cmbsf)", y= "Hydrgen (nM)") +
  coord_flip()



KF_plot_N<-ggplot(KF, aes(x=Depth, y=d15N, shape=Station, color=Station)) +

  geom_point(aes(size=7)) +
  scale_shape_manual(values=c(15,3)) +
  scale_color_manual(values=c("salmon", "salmon4")) +
  theme(text = element_text(size=20)) +
  scale_y_continuous(limits=c(0,14), breaks = c(0,2,4,6,8,10,12,14)) +
  scale_x_reverse() +
  geom_line() +
  coord_flip()

VK_plot_N<-ggplot(VK, aes(x=Depth, y=d15N, shape=Station, color=Station)) +
  geom_point(aes(size=7)) +
  scale_color_manual(values=c("turquoise3", "turquoise4")) +
  theme(text = element_text(size=20)) +
  scale_y_continuous(limits=c(0,14), breaks = c(0,2,4,6,8,10,12,14)) +
  scale_x_reverse() +
  geom_line() +
  coord_flip()



