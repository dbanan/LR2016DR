#data_prep_plant

#gather and prepare plant level harvest, PE, and mid season phenotype data 

###phenotype list: 
#panicle emergence 
#harvest weights 
#harvest architecture 
#midseason architecture? (height, tiller n)
#midseason PAI 


library(tools)
library(plyr)
library(ggplot2)
library(reshape2)
library(gridExtra)



#still need to fix path etc 

#####LOAD#####
#bring in the equivalent of step2 data (formatted, duplicates removed, outliers flagged)
setwd("/rsync/box/Setaria/2016 Setaria/data analysis")
load("step2_16DR_harvest_traits.RData")
load("step2_16DR_harvest_weights.RData")
load("step2_16DR_panicle_emergence.RData")
load("step2_16DR_midseason_architecture.RData")



#####FORMAT#####
#harvest traits 
together_l3$abb<-together_l3$trait
together_l3$trait[together_l3$trait=="BN"]<-"branch_number"
together_l3$trait[together_l3$trait=="BW"]<-"basal_width"
together_l3$trait[together_l3$trait=="CH"]<-"culm_height"
together_l3$trait[together_l3$trait=="TH"]<-"tiller_height"
together_l3$trait[together_l3$trait=="TN"]<-"tiller_number"

harvest_traits<-together_l3[,c(1,5,6,7,8,10,11)]

#get to subplot_id average 
harvest_traits1<-ddply(harvest_traits, c("subplot_id","rep","treatment","genotype","trait"), summarise, average=mean(data))
colnames(harvest_traits1)[6]<-"data"

#harvest weights 
harvest_weights<-stack_l2[,c(1,2,4,5,6,7)]

#panicle emergence
panicle_emergence<-pe16.d1[,c(1,2,3,5,8,19)]
colnames(panicle_emergence)[6]<-"data"

#midseason architecture 
clean$data<-as.numeric(clean$data)
#get subplot_id averages 
clean1<-ddply(clean, c("subplot_id", "rep", "treatment", "genotype", "trait"), summarise, average=mean(data))
colnames(clean1)[6]<-"data"

