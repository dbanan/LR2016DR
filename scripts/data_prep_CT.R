#data_prep_CT.R 

#get Parthi's CT data for subset15 to compare to leaf rolling and productivity traits 

library(reshape2)
library(ggplot2)

#infile
CT16dr<-read.csv("./data/raw_data/16DR_canopy_temp.csv", header=T, stringsAsFactors=F, na.strings=".")

#planting is 163 DOY 
#July 21 is 40 DAS
#July 11 is 30 DAS

#trim to relevant columns 
CT16dr<-CT16dr[,c(1,4,5,7,8)]
#re-do colnames 
colnames(CT16dr)<-c("subplot_id","treatment","genotype","CT40DAS","CT30DAS")

#long by trait 
CT16dr1<-melt(CT16dr, id.vars=c("subplot_id","treatment","genotype"), measure.vars=c("CT40DAS","CT30DAS"), variable.name="trait", value.name="data")

#quick visualization 
ggplot(data=CT16dr1, aes(treatment, data))+
  geom_boxplot(aes(fill=treatment))+
  facet_wrap(~trait)


###OUTPUT###
save(CT16dr1, file="./data/clean_data/data_subset_CT.Rdata")



