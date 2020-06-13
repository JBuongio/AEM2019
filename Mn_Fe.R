#Joy Buongiorno
#jbuongior21@gmail.com
#Plot for Mn and Fe concetrations

library(ggplot2)
library(scales)
show_col(hue_pal()(9))
cols<-c("AB.2" = "#F8766D","AB.1" = "#F8766D", "AC.1" = "#C77CFF", "AC.2" = "#C77CFF", "AB" = "#F8766D", "AC" = "#C77CFF", "HA" = "#00B9E3")
getwd()
Fe<-read.csv("Fe.csv")
VK_Fe_plot<-ggplot(Fe, aes(x=Depth, y=Fe, shape=Station, color=Station)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Station), size=0.75) +
  coord_flip() +
  labs(x="Depth (cmbsf)", y= "Fe (nM)")

VK_Mn<-ggplot(Fe, aes(x=Depth, y=Mn, shape=Station, color=Station)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Station), size=0.75) +
  coord_flip() +
  labs(x="Depth (cmbsf)", y= "Mn (nM)")

