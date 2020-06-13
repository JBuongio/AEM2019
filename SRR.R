# figure of SRR with depth #
library(ggplot2)
SRR<-read.csv("SRR.csv", header=TRUE)

######### Both together ##############
ggplot(SRR, aes(x = Depth, y = SRR, color=Station, shape= Station)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  theme(text = element_text(size=20)) +
  scale_x_reverse() +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Sulfate Reduction Rate (nmol cm-3 d-1)") +
  coord_flip()

########### separated ###############
SRR_AB<-subset(SRR, Station %in% c("AB"))
ggplot(SRR_AB, aes(x = Depth, y = SRR, color=Station, shape= Station)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(21)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(25, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  labs(x="Depth (cmbsf)", y= "Sulfate Reduction Rate (nmol cm-3 d-1)") +
  coord_flip()

SRR_AC<-subset(SRR, Station %in% c("AC"))

ggplot(SRR_AC, aes(x = Depth, y = SRR, color=Station, shape= Station)) + 
  geom_point(aes(fill=Station), colour="black", size=4, stroke=2) +
  scale_shape_manual(values = c(24)) +
  scale_fill_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=cols) +
  theme(text = element_text(size=20)) +
  scale_x_reverse(limits=c(25, 0)) +
  geom_line(aes(color=Station, linetype=Station)) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=SRR_AC$SRR-SRR_AC$sd, ymax=SRR_AC$SRR+SRR_AC$sd)) +
  labs(x="Depth (cmbsf)", y= "Sulfate Reduction Rate (nmol cm-3 d-1)") +
  coord_flip()
