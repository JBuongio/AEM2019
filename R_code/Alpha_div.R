#Joy Buongiorno
#jbuongior21@gmail.com
#calculate alpha diversity

library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)

AB1_alpha<-read.csv("Alpha_AB1.csv") #updated with normalized data at 60,000 sequences
head(AB1_alpha)
AB1_aplha_tidy<-gather(AB1_alpha, sample, Sobs, AB.1.0.1: AB.1.9.10)
head (AB1_aplha_tidy)
AB1_alpha_plot<-ggplot(AB1_aplha_tidy, aes(x=Number_of_sequences, y=Sobs, color=sample)) +
  geom_point(aes(fill=sample))


AB2_alpha<-read.csv("Alpha_AB2.csv")
head(AB2_alpha)
AB2_aplha_tidy<-gather(AB2_alpha, sample, Sobs, AB.2.1.2: AB.2.9.10)
head (AB2_aplha_tidy)
AB2_alpha_plot<-ggplot(AB2_aplha_tidy, aes(x=Number_of_sequences, y=Sobs, color=sample)) +
  geom_point(aes(fill=sample))

AC1_alpha<-read.csv("Alpha_AC1.csv")
head(AC1_alpha)
AC1_aplha_tidy<-gather(AC1_alpha, sample, Sobs, AC.1.0.1: AC.1.8.9)
head (AC1_aplha_tidy)
AC1_alpha_plot<-ggplot(AC1_aplha_tidy, aes(x=Number_of_sequences, y=Sobs, color=sample)) +
  geom_point(aes(fill=sample))

AC2_alpha<-read.csv("Alpha_AC2.csv")
head(AC2_alpha)
AC2_aplha_tidy<-gather(AC2_alpha, sample, Sobs, AC.2.17.18: AC.2.8.9)
head (AC2_aplha_tidy)
AC2_alpha_plot<-ggplot(AC2_aplha_tidy, aes(x=Number_of_sequences, y=Sobs, color=sample)) +
  geom_point(aes(fill=sample))
