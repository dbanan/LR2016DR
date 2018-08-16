#data_prep_LR

#gather and prepare leaf rolling data from 2016 field experiment 
#15 genotype in depth test
#whole population score survey 


library(tools)
library(plyr)
library(ggplot2)
library(reshape2)
library(gridExtra)


#this script should output (bottom of script) cleaned, workable datasets for:
# (1) population level leaf roll score
# (2) 15 genotype subset leaf roll evaluation traits (hemi, stamp, roll)




###DESIGN###----
#infile experimental design
design16dr<-read.csv("./data/raw_data/16DR_subplot_genotype_list_FINAL.csv",header=T, stringsAsFactors=F,na.strings=".")
#format column names
colnames(design16dr)<-c("subplot_id","rep","awning","treatment","genotype")


###HEMI###----
#infile image id key
hemikey<-read.csv("./data/raw_data/16DR_LR_hemi_key.csv",header=T,stringsAsFactors=F)
#infile hemiview output values 
hemiraw<-read.csv("./data/raw_data/16DR_LR_hemi_data.csv",header=T,stringsAsFactors=F)

#remove first column from each (dummy count)
hemikey<-hemikey[-c(1)]
hemiraw<-hemiraw[-c(1)]

#match hemi data to key
hemiraw1<-merge(hemikey,hemiraw,by=c("Label"))
hemiraw1$subsample[hemiraw1$camera=="LGP_001"]<-1
hemiraw1$subsample[hemiraw1$camera=="LGP_002"]<-2
colnames(hemiraw1)[6]<-"time"
hemiraw1$date[hemiraw1$time=="midday"]<-"7/22/2016"
hemiraw1$date[hemiraw1$time=="dawn"]<-"7/23/2016"

#**gonna keep this in for now** (4/2/18)
#remove outlier ID'd in boxplots below?
#hemiraw1<-hemiraw1[!(hemiraw1$LAI>3),]

#plot camera performance (evaluated by LAI and GSF)
#split data camera wide
camera.LAI<-dcast(hemiraw1, subplot_id+time~camera, value.var="LAI")
camera.GSF<-dcast(hemiraw1, subplot_id+time~camera, value.var="GSF")
camLAI<-ggplot(camera.LAI,aes(x=LGP_001,y=LGP_002))+geom_point()+labs(title="LAI camera comparison")+geom_smooth(method="lm")
camGSF<-ggplot(camera.GSF,aes(x=LGP_001,y=LGP_002))+geom_point()+labs(title="GSF camera comparison")+geom_smooth(method="lm")

#boxplots at individual observation level 
hemiraw2<-merge(hemiraw1, design16dr, by=c("subplot_id"))
ggplot()+geom_boxplot(data=hemiraw2, aes(factor(time), LAI, fill=factor(treatment)))
ggplot()+geom_boxplot(data=hemiraw2, aes(factor(time), GSF, fill=factor(treatment)))
ggplot()+geom_boxplot(data=hemiraw2, aes(factor(time), DSF, fill=factor(treatment)))
ggplot()+geom_boxplot(data=hemiraw2, aes(factor(time), ISF, fill=factor(treatment)))

#trim
hemi<-hemiraw1[,c(5,6,33,10,12,14,28)]

#idenfity the subplots used in the 15 genotypes subset 
subset15<-unique(hemi$subplot_id)

#export plots
camLAI
camGSF

#hemi END



###POPULATION SCORE###----
#infile roll score data
scoreraw<-read.csv("./data/raw_data/Behavior_962016_095456.csv",header=F, stringsAsFactors=F,na.strings=".")
#column names
colnames(scoreraw)[2]<-"trait"
colnames(scoreraw)[3]<-"subplot_id"
colnames(scoreraw)[4]<-"subsample"
colnames(scoreraw)[5]<-"score"
#make score numeric
scoreraw$score<-as.numeric(scoreraw$score)
#subset to just entries with expected data (leaf rolling trait, first and only Janam subsample, dry plots)
scoreraw<-subset(scoreraw, trait=="leaf_rolling")
scoreraw<-subset(scoreraw, subsample=="1")
#merge with experimental design
scoreraw<-merge(scoreraw, design16dr, by=c("subplot_id"))
#just dry plots 
scoreraw<-subset(scoreraw, treatment=="dry")
#remove B100, BLANK, and TB_0430 
scoreraw<-scoreraw[!scoreraw$genotype=="B100",]
scoreraw<-scoreraw[!scoreraw$genotype=="BLANK",]
scoreraw<-scoreraw[!scoreraw$genotype=="TB_0430",]
#remove NAs 
scoreraw<-scoreraw[complete.cases(scoreraw$score),]

#distribution at subplot_id level
hist(scoreraw$score)

#genotype averages 
scoreave<-ddply(scoreraw,c("genotype"),summarise,N=length(score),mean=mean(score),min=min(score),max=max(score))

#population distribution of genotype averages 
ggplot(scoreave, aes(factor(mean)))+geom_bar()

#plot rank order
scoreave<-scoreave[order(scoreave$mean),]
scoreave$genotype1<-factor(scoreave$genotype, levels=unique(as.character(scoreave$genotype)))

score_rank<-ggplot(scoreave, aes(x=genotype1, y=mean))+geom_point()+geom_errorbar(aes(ymin=min, ymax=max))+
  labs(title="Plot roll score rank ordered",x="genotype",y="score")+
  theme(panel.background=element_rect(fill=NA),
        panel.border=element_rect(fill=NA,color="black"),
        axis.text=element_text(size=8),
        axis.title=element_text(size=14), 
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, hjust = 1)
  )
score_rank


#####COMPRESSED SCORE#####
#compress score scale to binary (0,1) or ternary (0,1,2) or (0,1,2,3) at subplot level 

#visual score 
visual_score<-scoreraw[,c(1,15,14,12,5)]

#trial various compressed score scales 
visual_score$score4<-visual_score$score
visual_score$score4[visual_score$score>3]<-3

visual_score$score3<-visual_score$score
visual_score$score3[visual_score$score>2]<-2

visual_score$score2<-visual_score$score
visual_score$score2[visual_score$score>1]<-1

visual_score1<-melt(visual_score, id.vars=c("subplot_id", "rep", "treatment", "genotype"), measure.vars=c("score", "score2", "score3", "score4"), variable.name="trait",  value.name="data")



###SUBSET SCORE###----
#SCORE FOR 15 GENOTYPE SUBSET 
#subset to observations that match with stamp and hemi 
scoresub<-subset(scoreraw, subplot_id %in% subset15)
scoresub$time<-"midday"
scoresub$date<-"7/21/2016"
scoresub$trait<-"score"

#trim
scorelong<-scoresub[,c(1,16,4,3,5)]
colnames(scorelong)[5]<-"data"
rownames(scorelong)<-c()

#subset score distribution
hist(scoresub$score)

#score END



###STAMP###----
#infile roll and inclination data data 
stampraw<-read.csv("./data/raw_data/16DR leaf rolling heat dome diurnal.csv",header=T, stringsAsFactors=F,na.strings=".")
#make data columns numeric
stampraw$inclination<-as.numeric(stampraw$inclination)
stampraw$obtuse<-as.numeric(stampraw$obtuse)
#correct obtuse roll angles, make all values positive
stampraw$roll1<-abs(stampraw$roll-stampraw$obtuse)

#trim the fat, rename columns 
stamp<-stampraw[,c(2,4,6,8,10)]
colnames(stamp)[5]<-"roll"

#stamp END


###SUBSET COMBO###----
#make them all long
hemilong<-melt(hemi,id.vars=c("subplot_id","time","subsample"),measure.vars=c("DSF","ISF","GSF","LAI"),variable.name="trait",value.name="data")
stamplong<-melt(stamp,id.vars=c("subplot_id","time","subsample"),measure.vars=c("inclination","roll"),variable.name="trait",value.name="data")
#stack 
combo<-rbind(hemilong, stamplong, scorelong)
combo$subplot_id<-as.numeric(combo$subplot_id)
#all the boxplots
combo1<-merge(combo, design16dr, by=c("subplot_id"))
combobox<-ggplot()+geom_boxplot(data=combo1, aes(factor(time), data, fill=factor(treatment)))+facet_wrap(~trait, scale="free")+theme(legend.position="none")
combobox

#trim boxplots to just traits of interest, add zeros to score
combo2<-subset(combo1, trait=="GSF"|trait=="inclination"|trait=="roll"|trait=="score")

#dummy score data 
plots<-subset(design16dr, subplot_id %in% subset15)
plots_d<-plots
plots_d$time<-"dawn"
plots_d$trait<-"score"
plots_d$data<-0
plots_d$subsample<-1
plots_d$rep<-1
plots_m<-subset(plots, treatment=="wet")
plots_m$time<-"midday"
plots_m$trait<-"score"
plots_m$data<-0
plots_m$subsample<-1
plots_m$rep<-1


#join dummy score data with rest of subset traits 
combo3<-rbind(combo2, plots_d, plots_m)

#plots
combobox_trim<-ggplot()+geom_boxplot(data=combo3, aes(factor(time), data, fill=factor(treatment)))+facet_wrap(~trait, scale="free")+theme(legend.position="none")
combobox_trim

box_GSF<-ggplot(subset(combo3,trait %in% c("GSF")))+geom_boxplot(aes(factor(time), data, fill=factor(treatment)))+labs(list(title="(C) Global site factor",x=NULL,y="proportion PPFD"))+theme(legend.position="none")+
  theme(panel.background=element_rect(fill=NA),
        panel.border=element_rect(fill=NA,color="black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        title=element_text(size=14)
  )
box_inclination<-ggplot(subset(combo3,trait %in% c("inclination")))+geom_boxplot(aes(factor(time), data, fill=factor(treatment)))+labs(list(title="(A) Leaf inclination angle",x=NULL,y="degrees"))+theme(legend.position="none")+
  theme(panel.background=element_rect(fill=NA),
        panel.border=element_rect(fill=NA,color="black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        title=element_text(size=14)
  )
box_roll<-ggplot(subset(combo3,trait %in% c("roll")))+geom_boxplot(aes(factor(time), data, fill=factor(treatment)))+labs(list(title="(B) Leaf roll angle",x=NULL,y="degrees"))+theme(legend.position="none")+
  theme(panel.background=element_rect(fill=NA),
        panel.border=element_rect(fill=NA,color="black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        title=element_text(size=14)
  )
box_score<-ggplot(subset(combo3,trait %in% c("score")))+geom_boxplot(aes(factor(time), data, fill=factor(treatment)))+labs(list(title="(D) Plot roll score",x=NULL,y="score"))+theme(legend.position="none")+
  theme(panel.background=element_rect(fill=NA),
        panel.border=element_rect(fill=NA,color="black"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        title=element_text(size=14)
  )

box_GSF
box_inclination
box_roll
box_score

grid.arrange(box_inclination, box_roll, box_GSF, box_score)




###OUTPUT###
#output .Rdata for population level score data and subset evaluation data 

save(combo, file="./data/clean_data/data_subset_LR.Rdata")
save(visual_score1, file="./data/clean_data/data_population_score.Rdata")











