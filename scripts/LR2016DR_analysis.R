#LR2016DR_analysis 


library(plyr)
library(ggplot2)
library(reshape2)
library(GGally)
library(lsmeans)



save(combo, file="./data/clean_data/data_subset_LR.Rdata")
save(visual_score1, file="./data/clean_data/data_population_score.Rdata")
save(combine, file="./data/clean_data/data_population_plant.Rdata")



load("./data/clean_data/data_subset_LR.Rdata")
load("./data/clean_data/data_population_score.Rdata")
load("./data/clean_data/data_population_plant.Rdata")

#analyze in order of population-plant/canopy-level to subset-leaf/canopy-level 
#communicate in reverse order (subset-leaf/canopy to population-plant/canopy)
#

#####POPULATION#####


#trim population-plant level data to just traits of interest 
interest<-c("basal_width","branch_number","culm_height","tiller_height","tiller_number",
            "per_plant_leaf_mass","per_plant_panicle_mass","per_plant_stem_mass","per_plant_vegetative_mass","per_plant_total_mass",
            "reproductive_vegetative_ratio","panicle_emerge","leaf_number_total",
            "height","tiller_count")
combine1<-subset(combine, trait %in% interest)

#trait abbreviations 
abb<-c("BW","BN","CH","TH","TN","LM","PM","SM","VM","TM","RV","PE","LN","MH","MT")
key<-cbind(interest, abb)
colnames(key)<-c("trait","abb")
combine1<-merge(combine1, key, by="trait")

#genotype averages (harvest and score)
combineg<-ddply(combine1, c("genotype", "treatment", "trait", "abb"), summarise, average=mean(data))
visualg<-ddply(visual_score, c("genotype"), summarise, score=mean(data))

png(file="./results/score_hist.png")
hist(visualg$score)
dev.off()


#make wide by treatment (to calculate treatment differences and identify(?) absolutes)
combinew<-dcast(combineg, genotype+trait+abb~treatment, value.var="average")

#calculate differences 
combinew$dw_diff<-combinew$dry-combinew$wet
combinew$abs_diff<-abs(combinew$dw_diff)
combinew$rel_diff<-combinew$dw_diff/combinew$wet
combinew$dw_percent<-combinew$dry/combinew$wet
combinew$log_percent<-log(combinew$dw_percent)
combinew$mid_diff<-combinew$dw_diff/((combinew$dry+combinew$wet)/2)

  
#long by data type (dry, wet, difference, etc) 
combinel<-melt(combinew, id.vars=c("genotype","trait","abb"), measure.vars=c("dry","wet","dw_diff","abs_diff","rel_diff","dw_percent","log_percent","mid_diff"), variable.name="type",value.name="data")
#trait_type combination 
combinel$abb_type<-paste(combinel$abb, combinel$type, sep="_")

#go wide again by trait+type combination 
combineww<-dcast(combinel, genotype~abb_type, value.var="data")

#join with score data
combineww<-merge(combineww, visualg, by=c("genotype"))


#show some correlations 
pairs(combineww[,c("score","CH_abs_diff","CH_dry","CH_dw_diff","CH_dw_percent","CH_log_percent","CH_mid_diff","CH_rel_diff","CH_wet")])

pairs(combineww[,c("score","TM_abs_diff","TM_dry","TM_dw_diff","TM_dw_percent","TM_log_percent","TM_mid_diff","TM_rel_diff","TM_wet")])

pairs(combineww[,c("score","LM_abs_diff","LM_dry","LM_dw_diff","LM_dw_percent","LM_log_percent","LM_mid_diff","LM_rel_diff","LM_wet")])

pairs(combineww[,c("score",
                   "TM_dry",
                   "VM_dry",
                   "SM_dry",
                   "LM_dry",
                   "PM_dry",
                   "CH_dry",
                   "TN_dry",
                   "BN_dry",
                   "LN_dry",
                   "RV_dry",
                   "PE_dry"
                   )])

pairs(combineww[,c("score",
                   "TM_wet",
                   "VM_wet",
                   "SM_wet",
                   "LM_wet",
                   "PM_wet",
                   "CH_wet",
                   "TN_wet",
                   "BN_wet",
                   "LN_wet",
                   "RV_wet",
                   "PE_wet"
)])

pairs(combineww[,c("score",
                   "TM_log_percent",
                   "VM_log_percent",
                   "SM_log_percent",
                   "LM_log_percent",
                   "PM_log_percent",
                   "CH_log_percent",
                   "TN_log_percent",
                   "BN_log_percent",
                   "LN_log_percent",
                   "RV_log_percent",
                   "PE_log_percent"
)])





ggplot(data=combineww, aes(x=score, y=TM_dry))+geom_point()+geom_smooth(method="lm")
fit<-lm(formula=combineww$TM_dry~combineww$score)
summary(fit)




ggplot(data=combineww, aes(x=score, y=TM_dry))+geom_point()+geom_smooth(method="lm")
ggplot(data=combineww, aes(x=score, y=VM_dry))+geom_point()+geom_smooth(method="lm")
ggplot(data=combineww, aes(x=score, y=LM_dry))+geom_point()+geom_smooth(method="lm")
ggplot(data=combineww, aes(x=score, y=SM_dry))+geom_point()+geom_smooth(method="lm")


ggplot(data=combineww, aes(x=score, y=BN_dry))+geom_point()+geom_smooth(method="lm")

ggplot(data=combineww, aes(x=factor(score), y=TM_dry))+
  geom_jitter(width=0.12, aes(size=RV_dry))

ggplot(data=combineww, aes(x=factor(score), y=CH_dry))+
  geom_jitter(width=0.12)




fit<-lm(formula=combineww$VM_dry~combineww$score)
summary(fit)




#report means and number of observations 
check<-aggregate(VM_dry~score, combineww, mean)
number<-aggregate(VM_dry~score, combineww, length)

#multiple comparisons, try least square means, more robust
model = lm(VM_dry ~ factor(score),
           data=combineww)

leastsquare = lsmeans(model,
                      pairwise ~ factor(score),
                      adjust = "tukey")

grouping<-cld(leastsquare,
    alpha   = 0.05,
    Letters = letters,
    adjust="tukey")

#place grouping letters above boxes 
ggplot(data=combineww, aes(factor(score), VM_dry))+
  #geom_boxplot()+
  geom_jitter(width=0.125)+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(VM_dry, 2), y=VM_dry), size=3)+
  geom_text(data=number, aes(label=VM_dry, y=0))+
  geom_text(data=grouping, aes(label=.group, x=factor(grouping$score), y=Inf))+
  theme_classic()

png(file="./results/TMdry_v_score.png")
ggplot(data=combineww, aes(factor(score), VM_dry))+
  #geom_boxplot()+
  geom_jitter(width=0.125)+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(VM_dry, 2), y=VM_dry), size=3)+
  geom_text(data=number, aes(label=VM_dry, y=0))+
  geom_text(data=grouping, aes(label=.group, x=factor(grouping$score), y=Inf))+
  theme_classic()
dev.off()












ggplot(data=combinel, aes(x=data))+
  geom_histogram()+
  facet_wrap(~abb_type, scales="free")


trythis<-subset(combinel, type==c("dry","wet","dw_diff"))


ggplot(data=trythis, aes(x=data, color=type))+
  geom_freqpoly()+
  facet_wrap(~trait, scales="free")








#subset population level 
#subplot_id's of interest 
subset15<-unique(combo2$subplot_id)







#correlations of population plant level traits 
#is leaf rolling limited to any specific architectural or developmental "ideotype"? 
#doesn't really look like it :( 




#ideal model: leaf_rolling = genotype * water * time (where time is seasonal)
#true model: leaf_rolling = genotype? 
#^but does this really make any sense? only have leaf rolling score for dry plots 
#two biological replicates for most traits 
#have ~207 genotypes 


#15 genotype subset leaf/canopy time and treatment dynamics aka interactions 

#model: leaf_rolling = genotype * water * time (where time is diurnal) 
#use psudeo replication? 
#one true biological replicate (treatment applied to only one true biological replicate: a single 6*5 1cm square grid. drought by water withholding and well-watered by )
#have multiple measurements (subsamples) for some traits (roll angle, inclination angle)
#have 15 geno 

ggplot()+geom_boxplot(data=combo2, aes(factor(time), data, fill=factor(treatment)))+facet_wrap(~trait, scale="free")+theme(legend.position="none")

#calculate genotype averages from subset15 data 
combog<-ddply(combo2, c("genotype","treatment","time","trait"), summarise, average=mean(data))
#visualize 
ggplot()+
  geom_boxplot(data=combog, aes(factor(time), average, fill=factor(treatment)))+
  facet_wrap(~trait, scale="free")+
  theme(legend.position="none")


ggplot()+
  geom_jitter(data=combog, aes(interaction(treatment, time), average, color=treatment))+
  facet_wrap(~trait, scales="free")

png(file="./results/subset_boxplot.png")
ggplot()+
  geom_boxplot(data=combog, aes(interaction(treatment, time), average, fill=factor(treatment)))+
  #geom_jitter(data=combog, aes(interaction(treatment, time), average, color=treatment))+
  facet_wrap(~trait, scale="free")
dev.off()

#need to do some sort of ?t-test? ?ANOVA? here to tell if distributions are different...
#model is rolling = genotype * water * time 
# or is genotype our observation??? 


#raw data trait correlations? (color points by time and treatment)
#wide by trait 
combot<-dcast(combog, genotype+treatment+time~trait, value.var="average")
ggpairs(combot, columns=(4:8), aes(color=interaction(treatment, time), shape=time))

combot1<-subset(combot, time=="midday")
combot1<-subset(combot1, treatment=="dry")

ggpairs(combot1, columns=(4:8), aes(shape=time))



#correlation of time and treatment differences 
#calculate time change (inclination, roll, GSF in dry plots), wide by time 
combow_time<-dcast(combog, trait+genotype+treatment~time, value.var="average")
combow_time$diff_time<-combow_time$midday-combow_time$dawn
combow_time$per_time<-combow_time$midday/combow_time$dawn

#wide by trait (try percents)
combow_timew<-dcast(combow_time, genotype+treatment~trait, value.var="per_time")

combow_timew1<-subset(combow_timew, treatment=="dry")
ggpairs(combow_timew1, columns=(3:7), lower=list(continuous="smooth"))

png(file="./results/subset_percent_corr.png")
ggpairs(combow_timew1, columns=(3:7), lower=list(continuous="smooth"))
dev.off()

#remove score column (full of NA and Inf)
combow_timew1$score<-NULL


#join subset 15 time differences with score and LAI treatment differences and harvest TM, height etc 
plant15<-subset(combine1, subplot_id %in% subset15)
score15<-subset(visual_score, subplot_id %in% subset15)

#new score column merge to wide diurnal change dataset
score15<-score15[,c(2,5)]
colnames(score15)<-c("genotype","score")
combow_timew1<-merge(combow_timew1, score15, by=c("genotype"))







combol_time<-melt(combow_time, id.vars=c("genotype","treatment"), measure.vars=c("dawn","midday","diff_time","per_time"), variable.name="type",value.name="data")

combol_time$difference<-paste(combol_time$treatment, combol_time$type, sep="_")


#calculate treatment change (PAI), wide by treatment 
combow_trt<-dcast(combog, trait+genotype+time~treatment, value.var="average")
combow_trt$diff_trt<-combow_trt$dry-combow_trt$wet
combow_trt$per_trt<-combow_trt$dry/combow_trt$wet

#trim to treatment and rename traits, columns 
combow_trt1<-subset(combow_trt, time=="dawn")
combow_trt1<-subset(combow_trt1, trait=="LAI")


combow_time1<-combow_time[,c(1,2,6)]
combow_trt1<-combow_trt[,c(1,2,6)]

combow_both<-merge()




#correlate time change with treatment change 



