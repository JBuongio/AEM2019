#Joy Buongiorno
#jbuongior21@gmail.com
#Plot for qPCR data

library(ggplot2)
library(plyr)

qPCR_melted<-read.csv("qPCR_geochem.csv", header=TRUE, stringsAsFactors = FALSE)
AC_vals<-subset(qPCR_melted, Station %in% c("AC"))
AB_vals<-subset(qPCR_melted, Station %in% c("AB"))
AC_Bac<-subset(AC_vals, Domain %in% c("Bacteria"))
AC_Arc<-subset(AC_vals, Domain %in% c("Archaea"))
AB_Bac<-subset(AB_vals, Domain %in% c("Bacteria"))
AB_Arc<-subset(AB_vals, Domain %in% c("Archaea"))
VK_Bac<-subset(qPCR_melted, Domain %in% c("Bacteria"))
VK_Arc<-subset(qPCR_melted, Domain %in% c("Archaea"))

head (VK_Bac)


##########both bac and arch together is too cluttered#######
AC_qPCR<-ggplot(AC_vals[!is.na(AC_vals$Average),], aes(x=Depth, y=Average, color=Domain, shape=Domain)) +
  geom_point(aes(fill=Subcore), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=cols) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  geom_line(aes(color=Subcore, linetype=Domain)) +
  scale_y_log10(limits = c(1e2,1e12)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "qPCR (average copies/g sediment)") +
  coord_flip()
  
##########################All Archaeal data #########################################
Arc_qPCR<-ggplot(VK_Arc[!is.na(VK_Arc$Average),], aes(x=Depth, y=Average, color=Station, shape=Subcore)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 21, 24)) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Domain)) +
  scale_y_log10(limits = c(1e2,1e12)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Archaeal 16S (average copies/g sediment)") +
  coord_flip()

##########################All Bacterial data #########################################
Bac_qPCR<-ggplot(VK_Bac[!is.na(VK_Bac$Average),], aes(x=Depth, y=Average, color=Station, shape=Subcore)) +
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 21, 24)) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Domain)) +
  scale_y_log10(limits = c(1e4,1e12)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Bacterial 16S (average copies/g sediment)") +
  coord_flip()

######## downcore plots for figure ################


AB_Bacteria<-ggplot(AB_Bac[!is.na(AB_Bac$Average),], aes(x=Depth, y=Average, color=Subcore, Shape=Subcore)) +
  geom_point(aes(fill=Subcore, shape=Subcore),colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21,22)) +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_reverse(limits=c(20,0)) +
  geom_line(aes(color=Subcore)) +
  scale_y_log10(limits = c(1e4,1.5e11)) +
  geom_line(size=1.5) +
  labs(x="Depth (cmbsf)", y= "Bacteria (average copies/g sediment)") +
  coord_flip() +
  geom_text(aes(x = 2, y = 10, label = rout[[1]]), hjust = 0) +
  geom_text(aes(x = 2, y = 9.5, label = rout[[2]]), hjust = 0, parse = TRUE )


AC_Bacteria<-ggplot(AC_Bac[!is.na(AC_Bac$Average),], aes(x=Depth, y=Average, color=Subcore)) +
  geom_point(aes(fill=Subcore, shape=Subcore), colour="black", size=4, stroke=2) +
  scale_colour_manual(values=cols) +
  scale_shape_manual(values = c(23,24)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_reverse(limits=c(20,0)) +
  geom_line(aes(color=Subcore)) +
  scale_y_log10(limits = c(1e4,1.5e11)) +
  geom_line(size=1.5) +
  labs(x="Depth (cmbsf)", y= "Bacteria (average copies/g sediment)") +
  coord_flip() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed")


AB_Archaea<-ggplot(AB_Arc[!is.na(AB_Arc$Average),], aes(x=Depth, y=Average, color=Subcore)) +
  geom_point(aes(shape=Subcore), colour="black", size=4, stroke=2) +
  scale_colour_manual(values=cols) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_reverse(limits=c(20,0)) +
  geom_line(aes(color=Subcore)) +
  scale_y_log10(limits = c(1e4,1.5e11)) +
  geom_line(size=1.5) +
  labs(x="Depth (cmbsf)", y= "Archaea (average copies/g sediment)") +
  coord_flip()

AC_Archaea<-ggplot(AC_Arc[!is.na(AC_Arc$Average),], aes(x=Depth, y=Average, color=Subcore)) +
  geom_point(aes(shape=Subcore), colour="black", size=4, stroke=2) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values = c(23,24)) +
  scale_colour_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_reverse(limits=c(20,0)) +
  geom_line(aes(color=Subcore)) +
  scale_y_log10(limits = c(1e4,1.5e11)) +
  geom_line(size=1.5) +
  labs(x="Depth (cmbsf)", y= "Archaea (average copies/g sediment)") +
  coord_flip()

############ Regression and spearman 
#subset to get regression stats
AB1<-subset(AB_Bac, Subcore %in% c("AB.1")) 
AB2<-subset(AB_Bac, Subcore %in% c("AB.2"))
AC1<-subset(AC_Bac, Subcore %in% c("AC.1"))
AC2<-subset(AC_Bac, Subcore %in% c("AC.2"))

AB1A<-subset(AB_Arc, Subcore %in% c("AB.1")) 
AB2A<-subset(AB_Arc, Subcore %in% c("AB.2"))
AC1A<-subset(AC_Arc, Subcore %in% c("AC.1"))
AC2A<-subset(AC_Arc, Subcore %in% c("AC.2"))

# for R2 values
mod<-lm(Average ~ Depth, data = AC2) #Have to change each time (ie AB1 then AB2, AC1 then AC2.)

rout <- list(paste('Fitted model: ', round(coef(mod)[1], 3), ' + ',
                   round(coef(mod)[2], 3), ' x', sep = ''),
             paste('R^2 == ', round(summary(mod)[['r.squared']], 3),
                   sep = '')  )

#spearman corrleation of average value with depth
m <- cor.test(x=AC2A$Average, y=AC2A$Depth, method = "spearman") #change each time for each subset
n<-lm(AC1A$Average ~AC1A$Depth) #same
summary(n)



##########################################
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

########  Bac qPCR vs d13C, trendline thru all data points #############
VK_bac_d13C<-ggplot(VK_Bac, aes(x=qPCR, y=d13C)) +
  geom_point(aes(shape=Station, color=Subcore, size =4)) +
  theme(text = element_text(size=20)) +
  ylim(-28,-24) +
  scale_x_log10(limits=c(1e5,1e12)) +
  geom_smooth(method="lm") +
  labs(x="Bacterial 16S copies per gram sediment", y=expression(paste(delta^{13}, "C" [org],"(\u2030 vs. PDB)")))
VK_Bac_no_na<-na.omit(VK_Bac)
write.csv(VK_Bac_no_na, "VK_Bac.csv")
ggplotRegression(lm(d13C~qPCR, data=VK_Bac))
########  Bac qPCR vs d13C, no trendline #############
VK_bac_d13C<-ggplot(VK_Bac, aes(x=qPCR, y=d13C)) +
  geom_point(aes(shape=Station, color=Subcore, size=4)) +
  theme(text = element_text(size=20)) +
  ylim(-28,-24) +
  scale_x_log10(limits=c(1e5,1e12)) +
  labs(x="Bacterial 16S copies per gram sediment", y=expression(paste(delta^{13}, "C" [org],"(\u2030 vs. PDB)")))

########  Bac qPCR vs C/N, no trendline #############
VK_bac_CtoN<-ggplot(VK_Bac, aes(x=CtoN, y=d13C)) +
  geom_point(aes(shape=Station, color=Subcore, size =4)) +
  theme(text = element_text(size=20)) +
  ylim(-28,-24) +
  labs(x="Bacterial 16S copies per gram sediment", y="C/N")

VK_bac_CtoN<-ggplot(VK_Bac, aes(x=qPCR, y=CtoN, shape=Station, color=Subcore)) +
  geom_point(aes(shape=Station, color=Subcore)) +
  geom_point(size =4) +
  theme(text = element_text(size=20)) +
  scale_x_log10(limits=c(1e6,1e12)) +
  labs(y="C/N", x="Bacterial 16S copies per gram sediment")

######## Arch qPCR vs d13C ###################
VK_arc_d13C<-ggplot(VK_Arc, aes(x=qPCR, y=d13C)) +
  geom_point(aes(shape=Station, color=Subcore, size =4)) +
  theme(text = element_text(size=20)) +
  ylim(-28,-24) +
  scale_x_log10() +
  labs(x="Archaeal 16S copies per gram sediment", y=expression(paste(delta^{13}, "C" [org],"(\u2030 vs. PDB)")))

######## Arch qPCR vs C/N ###################
VK_arc_CtoN<-ggplot(VK_Arc, aes(x=qPCR, y=CtoN)) +
  geom_point(aes(shape=Station, color=Subcore, size =4)) +
  theme(text = element_text(size=20)) +
  scale_x_log10(limits=c(1e5,1e12))+
  labs(x="Archaeal 16S copies per gram sediment", y="C/N")
