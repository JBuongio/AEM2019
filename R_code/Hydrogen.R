#Joy Buongiorno
#jbuongior21@gmail.com
#Plots for hydrogen data

library(ggplot2)

#Kongsfjorden
hyd<-read.csv("Hydrogen.csv")
KF_H<-subset(hyd, Fjord %in% c("KF"))
head(KF_H)
#### finding outliers -- P has none, F has one right at the top
KFP_H<-subset(hyd, Station %in% c("P")) 
mod<-lm(Hydrogen ~ Depth , data=KFP_H)
cooksd<-cooks.distance(mod)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance, H in P")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels, these la
#### plotting, with outliers included
KF_H_plot<-ggplot(KF_H, aes(x=Depth, y=Hydrogen, shape=Station, color=Station)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(23, 22)) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(25, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Hydrogen (nM)") +
  coord_flip()

VK_H_plot<-ggplot(VK_H, aes(x=Depth, y=Hydrogen, shape=Station, color=Station)) + 
  geom_point(aes(size=7)) +
  theme(text = element_text(size=20)) +
  scale_color_manual(values=c("turquoise3", "turquoise4")) +
  geom_line() +
  scale_y_continuous(limits=c(0, 2)) +
  scale_x_reverse() +
  coord_flip()
  
####################################################
#Van Keulenfjorden
cols<-c("AB.2" = "#F8766D","AB.1" = "#F8766D", "AC.1" = "#C77CFF", "AC.2" = "#C77CFF", "AB" = "#F8766D", "AC" = "#C77CFF", "HA" = "#00B9E3")

H_AB<-subset(hyd, Station %in% c("AB"))

ggplot(H_AB, aes(x = Depth, y = Hydrogen, color=Station, shape= Rep)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(23, 22)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(30, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Hydrogen (nM)") +
  coord_flip()


H_AC<-subset(hyd, Station %in% c("AC"))

ggplot(H_AC, aes(x = Depth, y = Hydrogen, color=Station, shape= Rep)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(23, 22)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(30, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Hydrogen (nM)") +
  coord_flip()


H_HA<-subset(hyd, Station %in% c("HA"))

ggplot(H_HA, aes(x = Depth, y = Hydrogen, color=Station, shape= Rep)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(23, 22)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(30, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Hydrogen (nM)") +
  coord_flip()

######## All Van Keulenfjorden together
VK_hyd<-subset(hyd, Fjord %in% c("VK"))
VK_hyd_plot<-ggplot(VK_hyd, aes(x=Depth, y=Hydrogen, shape=Station, color=Station)) +
  geom_point(aes(shape=Station, fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Station), size=0.75) +
  labs(x="Depth (cmbsf)", y= "Hydrogen (nM)") +
  coord_flip()
