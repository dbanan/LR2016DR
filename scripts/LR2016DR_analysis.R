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
load("./data/clean_data/data_subset_CT.Rdata")



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




ggplot(data=combineww, aes(x=score, y=TM_dry))+geom_point()+geom_smooth(method="lm")
ggplot(data=combineww, aes(x=score, y=VM_dry))+geom_point()+geom_smooth(method="lm")
ggplot(data=combineww, aes(x=score, y=LM_dry))+geom_point()+geom_smooth(method="lm")
ggplot(data=combineww, aes(x=score, y=SM_dry))+geom_point()+geom_smooth(method="lm")
ggplot(data=combineww, aes(x=score, y=PM_dry))+geom_point()+geom_smooth(method="lm")



#score vs height (can make a fit line), poor r-squared but the line is significantly non-zero
combineww1<-combineww
combineww1$score[combineww1$score==3]<-2.5
ggplot(data=combineww1, aes(x=score, y=CH_dry))+geom_point()+geom_smooth(method="lm")
fit<-lm(formula=combineww1$CH_dry~combineww1$score)
summary(fit)

png(file="./results/score_v_CHdry.png")
ggplot(data=combineww1, aes(x=score, y=CH_dry))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black")+
  geom_text(aes(x=2, y=810, label="y = 47.76x + 337.30"))+
  geom_text(aes(x=2, y=790, label="p-value = 0.002"))+
  geom_text(aes(x=2, y=770, label="r-squared = 0.04"))+
  theme_classic()
dev.off()


#score vs biomass
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

png(file="./results/score_v_TMdry.png")
ggplot(data=combineww, aes(factor(score), TM_dry))+
  #geom_boxplot()+
  geom_jitter(width=0.125)+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(TM_dry, 2), y=TM_dry), size=3)+
  geom_text(data=number, aes(label=TM_dry, y=0))+
  geom_text(data=grouping, aes(label=.group, x=factor(grouping$score), y=Inf))+
  theme_classic()
dev.off()

png(file="./results/score_v_VMdry.png")
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

#try some sort of response plot 
png(file="./results/subset_response.png")
ggplot(data=combog, aes(interaction(time, treatment), average, group=interaction(genotype, treatment)))+
  geom_point(aes(color=treatment))+
  geom_line(alpha=0.1)+
  #geom_text(aes(label=ifelse(time=="midday", as.character(genotype), '')), hjust=0, size=2)+
  facet_wrap(~trait, scale="free")
dev.off()


ggplot(data=subset(combo2, time=="midday"), aes(interaction(treatment, genotype), data))+
  geom_boxplot(aes(fill=treatment))+
  facet_wrap(~trait, scales="free")
ggplot(data=subset(combo2, treatment=="dry"), aes(interaction(time, genotype), data))+
  geom_boxplot(aes(fill=time))+
  scale_fill_manual(values=c("purple","orange"))+
  facet_wrap(~trait, scales="free")


#need to do some sort of ?t-test? ?ANOVA? here to tell if distributions are different...
#model is rolling = genotype * water * time 
# or is genotype our observation??? 


#raw data trait correlations? (color points by time and treatment)
#wide by trait 
combot<-dcast(combog, genotype+treatment+time~trait, value.var="average")
ggpairs(combot, columns=(4:8), aes(color=interaction(treatment, time), shape=time), lower=list(continuous="smooth"))


combot1<-subset(combot, time=="midday")
combot1<-subset(combot1, treatment=="dry")

ggpairs(combot1, columns=(4:8), aes(shape=time))



#correlation of time and treatment differences 
#calculate time change (inclination, roll, GSF in dry plots), wide by time 
combow_time<-dcast(combog, trait+genotype+treatment~time, value.var="average")
combow_time$time_diff<-combow_time$midday-combow_time$dawn
combow_time$time_rel<-combow_time$time_diff/combow_time$dawn
combow_time$time_per<-combow_time$midday/combow_time$dawn
combow_time$time_logper<-log(combow_time$time_per)

#remove score (no sense in calculating difference)
combow_time<-combow_time[!(combow_time$trait=="score"),]

#combine trait and treatment 
combow_time$trait_trt<-paste(combow_time$treatment, combow_time$trait, sep="_time_")

#plot of differences 
#barplot by genotype individually 
ggplot(data=combow_time, aes(genotype, time_logper, fill=treatment))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~trait, scales="free")
#genotypes grouped as boxplot
ggplot(data=combow_time, aes(factor(treatment), time_diff, fill=treatment))+
  geom_boxplot()+
  facet_wrap(~trait, scales="free")
#genotypes individually as response plot 
ggplot(data=combow_time, aes(treatment, time_diff, group=genotype))+
  geom_point(aes(color=treatment))+
  geom_line(alpha=0.1)+
  facet_wrap(~trait, scale="free")




#productivity treatment differences 
#join subset 15 time differences with score and LAI treatment differences and harvest TM, height etc 
plant15<-subset(combine1, subplot_id %in% subset15)
score15<-subset(visual_score, subplot_id %in% subset15)

#new score column merge to wide diurnal change dataset
score15<-score15[,c(2,5)]
colnames(score15)<-c("genotype","score")


#calculate treatment change (productivity), wide by treatment 
#fetch LAI
laig<-subset(combog, time=="dawn" & trait=="LAI")
laig<-laig[,c(1,2,4,5)]
colnames(laig)[4]<-"data"
#fetch biomass
bmg<-subset(plant15, trait=="per_plant_total_mass")
bmg<-bmg[,c(1,4,5,6)]
#fetch canopy temperature
CT16dr2<-CT16dr1[,c(2,3,4,5)]

#join productivity and CT for calculating treatment differences 
prodg<-rbind(laig, bmg, CT16dr2)

#make productivity wide by treatment and calculate treatment differences 
prodw_trt<-dcast(prodg, genotype+trait~treatment, value.var="data")
prodw_trt$trt_diff<-prodw_trt$dry-prodw_trt$wet
prodw_trt$trt_per<-prodw_trt$dry/prodw_trt$wet





#join difference calculations with absolute data into a very looong data set 
#use "type" column as treatment / difference column
#absolute (all traits)
prodg1<-prodg
colnames(prodg1)[2]<-"type"
prodg1$trait<-as.character(prodg1$trait)
prodg1$trait[prodg1$trait=="LAI"]<-"LAIdawn"

combog1<-subset(combog, time=="midday")
colnames(combog1)[2]<-"type"
colnames(combog1)[5]<-"data"
combog1<-combog1[,c(1,2,4,5)]
combog1$trait<-as.character(combog1$trait)
combog1$trait[combog1$trait=="LAI"]<-"LAImidday"

CT16dr2<-CT16dr1
CT16dr2$data[CT16dr2$data>44]<-NA
colnames(CT16dr2)[2]<-"type"
CT16dr2<-CT16dr2[,c(2,3,4,5)]

abs_all<-rbind(prodg1, combog1, CT16dr2)

#time difference (rolling traits)
combow_time1<-combow_time
combow_time1$difference<-"time_diff"
combow_time1$type<-paste(combow_time1$difference, combow_time1$treatment, sep="_")
combow_time2<-combow_time1[,c(1,2,12,8)]
colnames(combow_time2)[4]<-"data"

#treatment difference (productivity traits)
prodw_trt1<-prodw_trt
prodw_trt1$type<-"trt_diff"
prodw_trt2<-prodw_trt1[,c(1,2,7,6)]
colnames(prodw_trt2)[4]<-"data"

#stack abs, time_diff, trt_diff 
subset_stack<-rbind(abs_all, combow_time2, prodw_trt2)
#paste type and trait into special trait to go wide with 
subset_stack$special<-paste(subset_stack$type, subset_stack$trait, sep="_")




#canopy temperature versus leaf rolling 
#leaf rolling traits expressed as a diurnal difference
rolling<-combow_time1[,c(1:3,8)]
colnames(rolling)[4]<-"data"
temp<-CT16dr1[,c(2:5)]
temp$data[temp$data>44]<-NA
compare<-rbind(rolling, temp)
#wide on trait 
comparew<-dcast(compare, genotype+treatment~trait, value.var="data")

ggpairs(comparew, columns=c(3:8), lower=list(continuous="smooth"))

ggpairs(data=subset(comparew, treatment=="dry"), columns=c(3:8), lower=list(continuous="smooth"))
ggpairs(data=subset(comparew, treatment=="wet"), columns=c(3:8), lower=list(continuous="smooth"))


ggpairs(comparew, columns=c(3:8), lower=list(continuous="smooth"), aes(color=treatment))

#roll change and canopy temp are correlated 
fitct<-lm(formula=comparew$CT30DAS~comparew$roll)
summary(fitct)
png(file="./results/roll_v_CT.png")
ggplot(data=comparew, aes(x=roll, y=CT30DAS))+
  geom_point(aes(color=treatment), size=2)+
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  geom_text(label="y=-5.87x+40.10", aes(x=0.3, y=34.5))+
  geom_text(label="***", aes(x=0.3, y=34))+
  geom_text(label="r-squared=0.45", aes(x=0.3, y=33.5))+
  theme_classic()
dev.off()




#correlations 
#wide by special 
specialw<-dcast(subset_stack, genotype~special, value.var="data")

#remove some crazy outliers 
specialw1<-specialw
specialw1$time_diff_dry_inclination[specialw1$time_diff_dry_inclination>2]<-NA

#huge matrix (can't see at all what is going on)
ggpairs(specialw, columns=(2:31), lower=list(continuous="smooth"))

#zoom in
ggpairs(specialw, columns=c(2:10), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(23:31), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(11:14,19:22,2,3,6,8,10), lower=list(continuous="smooth"))


ggpairs(specialw, columns=c(11:18), lower=list(continuous="smooth"))

ggpairs(specialw, columns=c(2,3,11:14), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(2,3,11:14,19,20), lower=list(continuous="smooth"))

ggpairs(specialw, columns=c(11:18,2,3,23,24), lower=list(continuous="smooth"))

ggpairs(specialw, columns=c(2,3,23,24,6,8,27,29,19:22), lower=list(continuous="smooth"))


ggpairs(specialw, columns=c(2,3,21,22,10:14), lower=list(continuous="smooth"))


ggplot(specialw, aes(x=time_diff_dry_roll, y=trt_diff_LAI))+
  geom_point()

ggplot(specialw, aes(x=time_diff_dry_roll, y=trt_diff_per_plant_total_mass))+
  geom_point()


ggplot(specialw)+
  geom_point(aes(x=time_diff_dry_roll, y=trt_diff_LAI), color="blue")+
  #geom_smooth(method="lm", aes(x=time_diff_dry_roll, y=trt_diff_LAI), color="black")+
  geom_point(aes(x=time_diff_dry_roll, y=trt_diff_per_plant_total_mass), color="red")
  #geom_smooth(method="lm", aes(x=time_diff_dry_roll, y=trt_diff_per_plant_total_mass), color="grey")


#dry rolling time deltas correlate with each other and with LRS 
png(file="./results/time_delta_dry_corr.png")
ggpairs(specialw, columns=c(10,14,13,11,12), lower=list(continuous="smooth"))
dev.off()

#are wet and dry time deltas correlated?...not super strongly
png(file="./results/time_delta_wetdry_corr.png")
ggpairs(specialw1, columns=c(11:18), lower=list(continuous="smooth"))
dev.off()

#dry time deltas are more correlted with biomass and biomass response 
ggpairs(specialw, columns=c(10:14,2,3), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(10:14,19,20), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(10:14,21,22), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(10:14,7,8), lower=list(continuous="smooth"))

png(file="./results/time_delta_biomass_corr.png", width=800, height=800)
ggpairs(specialw, columns=c(21,22,10,14,13,11,12), lower=list(continuous="smooth"))
dev.off()


#wet time deltas not really correlated with much
ggpairs(specialw, columns=c(15:18,19,20), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(15:18,23,24), lower=list(continuous="smooth"))
ggpairs(specialw, columns=c(15:18,21,22,28,29), lower=list(continuous="smooth"))

#try to find images where canopy changed alot due to leaf adjustment 
ggplot(data=specialw, aes(x=time_diff_dry_roll, y=time_diff_dry_LAI))+
  geom_point()+
  geom_text(aes(label=genotype))
  
  
ggplot(data=specialw, aes(x=time_diff_wet_inclination, y=time_diff_wet_LAI))+
  geom_point()+
  geom_text(aes(label=genotype))




#wide by just trait 
traitw<-dcast(subset_stack, genotype+type~trait, value.var="data")

ggpairs(traitw, columns=(3:12), lower=list(continuous="smooth"), aes(color=type))



#correlate time change with treatment change 
#wide by trait to enable correlations 
#time differences within treatment 
combow_timew<-dcast(combow_time, genotype~trait_trt, value.var="diff_time")
#treatment differences and absolute values? 
prodw_trtw<-dcast(prodw_trt, genotype~trait, value.var="diff_trt")

#merge time and treatment differences and score
bothw<-merge(combow_timew, prodw_trtw, by="genotype")
bothw<-merge(bothw, score15, by="genotype")








#correlation of wet and dry time differences 
ggpairs(combow_timew, columns=(2:9), lower=list(continuous="smooth"))

ggpairs(bothw, columns=(2:12), lower=list(continuous="smooth"))

ggpairs(bothw, columns=c(2:5,10:12), lower=list(continuous="smooth"))


