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





###Full pop### 
#do this again but on whole population data and relate to score
CTfull<-read.csv("./data/raw_data/16DR_canopyTemp_fullPop.csv", header=T, stringsAsFactors=F, na.strings=".")

#trim and rename columns 
CTfull<-CTfull[,c(1,4,5,7,8)]
colnames(CTfull)<-c("subplot_id","treatment","genotype","CT40DAS","CT30DAS")

#long by trait 
CTfulll<-melt(CTfull, id.vars=c("subplot_id","treatment","genotype"), measure.vars=c("CT40DAS","CT30DAS"), variable.name="trait", value.name="data")

#remove extreme temperature points 
CTfulll$data[CTfulll$data>41]<-NA

#calculate genotype average canopy temperature 
CTfullg<-ddply(CTfulll, c("genotype","treatment","trait"), summarise, average=mean(data,na.rm=TRUE))



######################
###STOMATAL DENSITY###
######################

SDfull<-read.csv("./data/raw_data/Setaria_stomata_2016_GWAS_phenotype.csv", header=T, stringsAsFactors=F, na.strings=".")

#P1 is July 14-20
#P2 is July 26 and 27 

SDfull<-SDfull[,c(1,4,5,7,8,9)]
colnames(SDfull)<-c("subplot_id","treatment","genotype","SDraw","collection","SDpred")

#hist() doesn't show major outliers. raw has a more normal distribution than pred 

#make raw and pred long as trait 
SDfulll<-melt(SDfull, id.vars=c("subplot_id","treatment","genotype","collection"), measure.vars=c("SDraw","SDpred"), variable.name="trait", value.name="data")

#take genotype*treatment*collection averages 
SDg<-ddply(SDfulll, c("genotype","treatment","collection","trait"), summarise, average=mean(data, na.rm=TRUE))









