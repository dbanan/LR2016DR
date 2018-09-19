#LR2016DR_analysis_180503

#this analysis will bring together and analyze the various 16DR datasets in preparation for a leaf rolling manuscript
#then builds on notes from 5/2/18 meeting with Andrew 
#...
#biomass/time calculation 
#relative difference transformations 
#non-parametric tests 

#AND organize work workspace 

#6/9/18 started version control

library(tools)
library(plyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(corrplot)

install.packages("car")
install.packages("lsmeans")
install.packages("multcompView")
install.packages("dunn.test")
library(lsmeans)
library(multcompView)
library(dunn.test)
library(car)

setwd("/rsync/box/Darshi work/STORY leaf movement")
save.image("LR2016DR_analysis_180511.Rdata")
load("LR2016DR_analysis_180511.Rdata")

#leaf rolling 
#panicle emergence 
#harvest weights 
#harvest architecture 
#midseason architecture? (height, tiller n)
#midseason PAI 

#figure goals 
#do leaves roll?
#how are leaf level and canopy level traits correlated? 
#does visual score capture these traits? 
#how does visual score vary across population? 
#is leaf rolling correlated with productivity? productivity response to drought? (subset and whole pop)

#####LOAD#####
#bring in the equivalent of step2 data (formatted, duplicates removed, outliers flagged)

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




#####COMBINE#####
#combine all population level datasets 
combine<-rbind(harvest_traits1, harvest_weights, panicle_emergence, clean1, visual_score1)

#remove B100, BLANK, and TB_0430? 
combine<-combine[!combine$genotype=="B100",]
combine<-combine[!combine$genotype=="BLANK",]
combine<-combine[!combine$genotype=="TB_0430",]

#remove extra panicle emergence data (stick with orignal)
combine<-combine[!(combine$trait=="panicle_emerge_DAS"),]


#simple boxplots of all 
ggplot()+geom_boxplot(data=combine, aes(factor(treatment), data, fill=factor(treatment)))+facet_wrap(~trait, scale="free")



#####TRAIT KEY#####
#trait name abbreviation key
trait<-as.character(unique(combine$trait))
abb<-c("BW","BN","CH","TH","TN","LM","IP","IC","IT","PM","SM","VM","TM","RV","HDAS","GN","LN","DN","GP",
       "H_LM","H_PM","H_SM","H_VM","H_TM","H_GN","H_LN","H_DN",
       "P_LM","P_PM","P_SM","P_VM","P_TM","P_GN","P_LN","P_DN,", 
       "PE","MH","MT","S0","S2","S3","S4")
key<-cbind(trait, abb)

combine<-merge(combine, key, by="trait")

#####SUBSET#####
#trim to 15 genotype subset (subset15)
combine_subset<-subset(combine, subplot_id %in% subset15)

#wide by treatment 
subsetw<-dcast(combine_subset, genotype+abb+trait~treatment, value.var="data")
#calculate differences 
subsetw$abs_diff<-subsetw$dry-subsetw$wet
subsetw$rel_diff<-subsetw$abs_diff/subsetw$wet

#try transformations of relative treatment difference 
#log 
subsetw$log_diff<-log(subsetw$rel_diff)
#log doesnt work with zero and negative values 
#sign*log(abs(rel_diff)
subsetw$slog_diff<-sign(subsetw$abs_diff)*(log(abs(subsetw$rel_diff)))

#cube root (complicated expression b/c of how R uses complex numbers? idk)
subsetw$cube_diff<-sign(subsetw$rel_diff)*(abs(subsetw$rel_diff))^(1/3)

#OK, now try alternaticves to our dogmatic "relative treatment difference" 
#percent tells you direction on either side of 100% rather than sign
#percent itself is right skewed 
subsetw$percent<-subsetw$dry/subsetw$wet
#log(percent)
subsetw$log_per<-log(subsetw$percent)

#what about a different way of calculating relative treatment difference? 
subsetw$mid_diff<-subsetw$abs_diff/((subsetw$dry+subsetw$wet)/2)
#but in this case, what is the significance of the mid-value (average)?


#wide by trait 
subset_dry<-subsetw[,c(1,2,3,4)]
subset_dry$name<-paste("d", subset_dry$abb, sep="_")
subset_dryw<-dcast(subset_dry, genotype~name, value.var="dry")

subset_wet<-subsetw[,c(1,2,3,5)]
subset_wet$name<-paste("w", subset_wet$abb, sep="_")
subset_wetw<-dcast(subset_wet, genotype~name, value.var="wet")

subset_rel<-subsetw[,c(1,2,3,7)]
subset_rel$name<-paste("r", subset_rel$abb, sep="_")
subset_relw<-dcast(subset_rel, genotype~name, value.var="rel_diff")

subset_mid<-subsetw[,c(1,2,3,13)]
subset_mid$name<-paste("m", subset_mid$abb, sep="_")
subset_midw<-dcast(subset_mid, genotype~name, value.var="mid_diff")

subset_lp<-subsetw[,c(1,2,3,12)]
subset_lp$name<-paste("l", subset_lp$abb, sep="_")
subset_lpw<-dcast(subset_lp, genotype~name, value.var="log_per")

subset_wide<-merge(subset_dryw, subset_wetw, by=c("genotype"))
subset_wide<-merge(subset_wide, subset_relw, by=c("genotype"))
subset_wide<-merge(subset_wide, subset_midw, by=c("genotype"))
subset_wide<-merge(subset_wide, subset_lpw, by=c("genotype"))




#need to join these subset harvest traits with subset rolling specific traits 
#combotime down to just dry values 
combotime1<-subset(combotime, treatment=="dry")

#add abbreviations
trait1<-unique(combotime1$trait)
abb1<-c("DSF", "GSF", "IA", "ISF", "LAI", "RA", "S0")
key1<-cbind(trait1, abb1)
colnames(key1)<-c("trait", "abb")

combotime1<-merge(combotime1, key1, by=c("trait"))

#cut traits not of interest 
interest<-c("GSF", "IA", "RA")
combotime1<-subset(combotime1, abb %in% interest)
  
#separate by data type, wide by trait 
#also pull out just midday and dawn values?

time_rel<-combotime1[,c(1,3,12,8)]
time_rel$name<-paste("r", time_rel$abb, sep="_")
time_relw<-dcast(time_rel, genotype~name, value.var="rel_diff")

time_mid<-combotime1[,c(1,3,12,9)]
time_mid$name<-paste("m", time_mid$abb, sep="_")
time_midw<-dcast(time_mid, genotype~name, value.var="mid_diff")

time_per<-combotime1[,c(1,3,12,10)]
time_per$name<-paste("p", time_mid$abb, sep="_")
time_perw<-dcast(time_per, genotype~name, value.var="percent")

time_lp<-combotime1[,c(1,3,12,11)]
time_lp$name<-paste("lp", time_lp$abb, sep="_")
time_lpw<-dcast(time_lp, genotype~name, value.var="log_per")

time_wide<-merge(time_relw, time_midw, by=c("genotype"))
time_wide<-merge(time_wide, time_perw, by=c("genotype"))
time_wide<-merge(time_wide, time_lpw, by=c("genotype"))

#also bring in dawn PAI values and treatment responses
time_wide<-merge(time_wide, PAI, by=c("genotype"))

#join with rest of harvest traits 
subset_time<-merge(time_wide, subset_wide, by=c("genotype"))

#LOOK AT IT 

ggplot(data=time_wide, aes(x=r_RA, y=PAI.dry))+geom_point()+geom_smooth(method="lm")
fit<-lm(formula=time_wide$PAI.trt.dawn~time_wide$m_GSF)
summary(fit)

ggplot(data=subset_time, aes(x=m_GSF, y=PAI.trt.dawn))+geom_point()+geom_smooth(method="lm")



ggplot(data=subset_time, aes(x=r_RA, y=l_TM))+geom_point()+geom_smooth(method="lm")

ggplot(data=subset_time, aes(x=d_S4, y=l_TM))+geom_point()+geom_smooth(method="lm")



ggplot(data=subset_time, aes(x=r_GSF, y=l_SM))+geom_point()+geom_smooth(method="lm")

fit<-lm(formula=subset_time$l_SM~subset_time$r_GSF)
summary(fit)



ggplot(data=subset_time, aes(x=PAI.trt.dawn, y=r_TM))+geom_point()+geom_smooth(method="lm")

fit<-lm(formula=subset_time$r_TM~subset_time$PAI.trt.dawn)
summary(fit)



ggplot(data=subset_time, aes(x=PAI.dry, y=d_TM))+geom_point()+geom_smooth(method="lm")

fit<-lm(formula=subset_time$d_TM~subset_time$PAI.dry)
summary(fit)


ggplot(data=subset_time, aes(x=PAI.wet, y=w_TM))+geom_point()+geom_smooth(method="lm")

fit<-lm(formula=subset_time$w_TM~subset_time$PAI.wet)
summary(fit)

#midseason PAI vs harvest biomass
ggplot(data=subset_time, aes(x=PAI.wet, y=w_TM))+geom_point()+geom_smooth(method="lm")
ggplot(data=subset_time, aes(x=PAI.dry, y=d_TM))+geom_point()+geom_smooth(method="lm")
ggplot(data=subset_time, aes(x=PAI.trt.dawn, y=r_TM))+geom_point()+geom_smooth(method="lm")

fit<-lm(formula=subset_time$d_TM~subset_time$PAI.dry)
summary(fit)

#rolling traits vs midseason PAI and harvest traits 
ggplot(data=subset_time, aes(x=PAI.dry, y=d_BW))+geom_point()+geom_smooth(method="lm")
ggplot(data=subset_time, aes(x=d_S4, y=l_BW))+geom_point()+geom_smooth(method="lm")


ggplot(data=subset_time, aes(x=PAI.dry, y=d_BW))+geom_point()+geom_smooth(method="lm")



ggplot(data=subset_time, aes(x=d_S4, y=r_RA))+geom_point()+geom_smooth(method="lm")

ggplot(data=subset_time, aes(x=r_GSF, y=r_RA))+geom_point()+geom_smooth(method="lm")
ggplot(data=subset_time, aes(x=l_P_TM, y=l_TM))+geom_point()+geom_smooth(method="lm")
ggplot(data=subset_time, aes(x=w_P_TM, y=w_TM))+geom_point()+geom_smooth(method="lm")


pairs(subset_time[,c("d_S4", "r_GSF", "r_RA", "r_IA", "PAI.trt.dawn", "l_P_TM", "l_TM")])

pairs(subset_time[,c("d_S4", "r_GSF", "r_RA", "r_IA", 
                     "PAI.trt.dawn", "PAI.wet", "PAI.dry",
                     "l_TM", "w_TM", "d_TM"
                     )])







#visualize score versus response 
png(file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180502 leaf rolling tests/subset_rTM.png")
ggplot(data=subset_wide, aes(x=d_S4, y=r_TM))+geom_point()+geom_smooth(method="lm")
dev.off()

fit<-lm(formula=subset_wide$r_VM~subset_wide$d_S0)
summary(fit)


ggplot(data=subset_wide, aes(x=d_S4, y=l_TH))+geom_point()+geom_smooth(method="lm")

fit<-lm(formula=subset_wide$l_TM~subset_wide$d_S4)
summary(fit)



check<-aggregate(l_CH~d_S4, subset_wide, mean)
number<-aggregate(l_CH~d_S4, subset_wide, length)

ggplot(subset_wide, aes(factor(d_S4), l_VM))+geom_boxplot()


#####GENOTYPE AVERAGES#####
#average by genotype*treatment (n=2)
combine_g<-ddply(combine, c("genotype", "treatment", "trait", "abb"), summarise, average=mean(data))
#wide genotype by treatment 
combine_gw<-dcast(combine_g, genotype+trait+abb~treatment, value.var="average")

#calculate relative treatment differences 
combine_gw$abs_diff<-combine_gw$dry-combine_gw$wet
combine_gw$rel_diff<-combine_gw$abs_diff/combine_gw$wet 

#try transformations of relative treatment difference 
#log 
combine_gw$log_diff<-log(combine_gw$rel_diff)
#log doesnt work with zero and negative values 
#sign*log(abs(rel_diff)
combine_gw$slog_diff<-sign(combine_gw$abs_diff)*(log(abs(combine_gw$rel_diff)))

#cube root (complicated expression b/c of how R uses complex numbers? idk)
combine_gw$cube_diff<-sign(combine_gw$rel_diff)*(abs(combine_gw$rel_diff))^(1/3)

#OK, now try alternaticves to our dogmatic "relative treatment difference" 
#percent tells you direction on either side of 100% rather than sign
#percent itself is right skewed 
combine_gw$percent<-combine_gw$dry/combine_gw$wet
#log(percent)
combine_gw$log_per<-log(combine_gw$percent)

#what about a different way of calculating relative treatment difference? 
combine_gw$mid_diff<-combine_gw$abs_diff/((combine_gw$dry+combine_gw$wet)/2)
#but in this case, what is the significance of the mid-value (average)?

#get rid of relative difference extreme outliers 
combine_gw$rel_diff[combine_gw$rel_diff>2.75]<-NA
#boxplot for relative differences 
ggplot()+geom_boxplot(data=combine_gw, aes(y=rel_diff, x=abb))

#stretch wide on everything so every trait has wet, dry, abs_diff, and rel_diff
combine_dry<-combine_gw[,c(1,2,3,4)]
combine_dry$name<-paste("d", combine_dry$abb, sep="_")
combine_dryw<-dcast(combine_dry, genotype~name, value.var="dry")

combine_wet<-combine_gw[,c(1,2,3,5)]
combine_wet$name<-paste("w", combine_wet$abb, sep="_")
combine_wetw<-dcast(combine_wet, genotype~name, value.var="wet")

combine_abs<-combine_gw[,c(1,2,3,6)]
combine_abs$name<-paste("a", combine_abs$abb, sep="_")
combine_absw<-dcast(combine_abs, genotype~name, value.var="abs_diff")

combine_rel<-combine_gw[,c(1,2,3,7)]
combine_rel$name<-paste("r", combine_rel$abb, sep="_")
combine_relw<-dcast(combine_rel, genotype~name, value.var="rel_diff")

combine_mid<-combine_gw[,c(1,2,3,13)]
combine_mid$name<-paste("m", combine_mid$abb, sep="_")
combine_midw<-dcast(combine_mid, genotype~name, value.var="mid_diff")

combine_lp<-combine_gw[,c(1,2,3,12)]
combine_lp$name<-paste("l", combine_lp$abb, sep="_")
combine_lpw<-dcast(combine_lp, genotype~name, value.var="log_per")

combine_wide<-merge(combine_dryw, combine_wetw, by=c("genotype"))
combine_wide<-merge(combine_wide, combine_absw, by=c("genotype"))
combine_wide<-merge(combine_wide, combine_relw, by=c("genotype"))
combine_wide<-merge(combine_wide, combine_midw, by=c("genotype"))
combine_wide<-merge(combine_wide, combine_lpw, by=c("genotype"))





ggplot(data=combine_wide, aes(x=d_BN, y=l_TM))+geom_point()+geom_smooth(method="lm")

ggplot(data=combine_wide, aes(x=d_S4, y=m_VM))+geom_point()+geom_smooth(method="lm")
ggplot(data=combine_wide, aes(x=d_S4, y=m_BW))+geom_point()+geom_smooth(method="lm")


ggplot(combine_wide, aes(factor(d_S4), d_VM))+geom_boxplot()
#report means and number of observations 
check<-aggregate(d_LN~d_S4, combine_wide, mean)
number<-aggregate(d_LN~d_S4, combine_wide, length)
#boxplot by score with mean and observation number reported
ggplot(data=combine_wide, aes(factor(d_S4), d_CH))+
  geom_boxplot()+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(d_CH, 2), y=d_CH+40), size=3)+
  geom_text(data=number, aes(label=d_CH), y=100)+
  theme_classic()



fit<-lm(formula=combine_wide$l_TM~combine_wide$d_TN)
summary(fit)


model=lm(r_TM~factor(d_S4), data=combine_wide)

kruskal.test(d_VM ~ factor(d_S4), data = combine_wide)

attach(combine_wide)

dunn.test(d_VM, d_S4, method="none", list=TRUE, kw=TRUE)

detach(combine_wide)






#####GENOTYPE REP AVERAGES#####
#average by genotype*rep*treatment (same as subplot) (individual observations) (n=1?) 
combine_gr<-ddply(combine, c("genotype","rep","treatment", "trait","abb"), summarise, average=mean(data))
#wide genotype*rep by treatment 
combine_grw<-dcast(combine_gr, genotype+rep+trait+abb~treatment, value.var="average")

#calculate relative treatment differences 
combine_grw$abs_diff<-combine_grw$dry-combine_grw$wet
combine_grw$rel_diff<-combine_grw$abs_diff/combine_grw$wet 

#percent 
combine_grw$percent<-combine_grw$dry/combine_grw$wet
#log(percent)
combine_grw$log_per<-log(combine_grw$percent)

#what about a different way of calculating relative treatment difference? 
combine_grw$mid_diff<-combine_grw$abs_diff/((combine_grw$dry+combine_grw$wet)/2)


#get rid of relative difference extreme outliers 
combine_grw$rel_diff[combine_grw$rel_diff>2.75]<-NA
#boxplot for relative differences 
ggplot()+geom_boxplot(data=combine_grw, aes(y=log_per, x=abb))

#stretch wide on everything so every trait has wet, dry, abs_diff, and rel_diff
combinegr_dry<-combine_grw[,c(1:4,5)]
combinegr_dry$name<-paste("d", combinegr_dry$abb, sep="_")
combinegr_dryw<-dcast(combinegr_dry, genotype+rep~name, value.var="dry")

combinegr_wet<-combine_grw[,c(1:4,6)]
combinegr_wet$name<-paste("w", combinegr_wet$abb, sep="_")
combinegr_wetw<-dcast(combinegr_wet, genotype+rep~name, value.var="wet")

combinegr_abs<-combine_grw[,c(1:4,7)]
combinegr_abs$name<-paste("a", combinegr_abs$abb, sep="_")
combinegr_absw<-dcast(combinegr_abs, genotype+rep~name, value.var="abs_diff")

combinegr_rel<-combine_grw[,c(1:4,8)]
combinegr_rel$name<-paste("r", combinegr_rel$abb, sep="_")
combinegr_relw<-dcast(combinegr_rel, genotype+rep~name, value.var="rel_diff")

combinegr_mid<-combine_grw[,c(1:4,11)]

combinegr_lp<-combine_grw[,c(1:4,10)]


combinegr_wide<-merge(combinegr_dryw, combinegr_wetw, by=c("genotype", "rep"))
combinegr_wide<-merge(combinegr_wide, combinegr_absw, by=c("genotype", "rep"))
combinegr_wide<-merge(combinegr_wide, combinegr_relw, by=c("genotype", "rep"))


#####SCORE AS FACTOR?#####

#####GT#####
#various scores as boxplot categories 
ggplot(combine_wide, aes(factor(d_S0), d_TM))+geom_boxplot()
ggplot(combine_wide, aes(factor(d_S4), d_TM))+geom_boxplot()
ggplot(combine_wide, aes(factor(d_S3), d_TM))+geom_boxplot()
ggplot(combine_wide, aes(factor(d_S2), d_TM))+geom_boxplot()
#where do NA's come from? (there are observations with no height or biomass data i think)


#report means and number of observations 
check<-aggregate(d_CH~d_S4, combine_wide, mean)
number<-aggregate(d_CH~d_S4, combine_wide, length)
#boxplot by score with mean and observation number reported
ggplot(data=combine_wide, aes(factor(d_S4), d_CH))+
  geom_boxplot()+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(d_CH, 2), y=d_CH+40), size=3)+
  geom_text(data=number, aes(label=d_CH), y=20)+
  theme_classic()

#do t-test to see if high score is truly different from low score 
x<-subset(combine_wide$d_CH, combine_wide$d_S4==0)
y<-subset(combine_wide$d_CH, combine_wide$d_S4==2)

t.test(x, y)

#describe the trend (linear) between score and biomass (or height, or trait response)
fit<-lm(formula=combine_wide$d_PE~combine_wide$d_S4)
summary(fit)

ggplot(data=combine_wide, aes(x=d_S4, y=d_PE))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()

ggplot(data=combine_wide, aes(x=d_PE, y=d_CH))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()



#####GRT######
#visualize, test genotype*rep*treatment version 
ggplot(combinegr_wide, aes(factor(d_S0), d_TM))+geom_boxplot()
ggplot(combinegr_wide, aes(factor(d_S4), d_TM))+geom_boxplot()
ggplot(combinegr_wide, aes(factor(d_S3), d_TM))+geom_boxplot()
ggplot(combinegr_wide, aes(factor(d_S2), d_TM))+geom_boxplot()


#do t-test to see if high score is truly different from low score 
x<-subset(combinegr_wide$d_CH, combinegr_wide$d_S4==0)
y<-subset(combinegr_wide$d_CH, combinegr_wide$d_S4==3)

t.test(x, y)

#one way ANOVA
results=aov(d_TM~factor(d_S4), data=combinegr_wide)
summary(results)

results<-aov(model)
summary(results)

#multiple comparisons (not really accounting for extreme unbalanced data)
TukeyHSD(results, conf.level = 0.95)

#try least square means, more robust
model = lm(r_CH ~ factor(d_S4),
           data=combinegr_wide)


leastsquare = lsmeans(model,
                      pairwise ~ factor(d_S4),
                      adjust = "tukey")

cld(leastsquare,
    alpha   = 0.05,
    Letters = letters,
    adjust="tukey")


png(file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180502 leaf rolling tests/scoreVsdCH.png")
ggplot(data=subset(combinegr_wide, !is.na(d_S4)), aes(factor(d_S4), d_CH))+
  geom_boxplot()+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(d_CH, 2), y=d_CH+40), size=3)+
  geom_text(data=number, aes(label=d_CH), y=20)+
  geom_text(data=abc_group, aes(label=abc, y=925))+
  theme_classic()
dev.off()

capture.output(summary(results), file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180502 leaf rolling tests/r_CHtest.csv")

png(file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180502 leaf rolling tests/scoreVsrCH.png")
ggplot(data=subset(combinegr_wide, !is.na(d_S4)), aes(factor(d_S4), r_CH))+
  geom_boxplot()
dev.off()


#report means and number of observations 
check<-aggregate(d_CH~d_S4, combinegr_wide, mean)
number<-aggregate(d_CH~d_S4, combinegr_wide, length)
#boxplot by score with mean and observation number reported
ggplot(data=combinegr_wide, aes(factor(d_S4), d_CH))+
  geom_boxplot()+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(d_CH, 2), y=d_CH+40), size=3)+
  geom_text(data=number, aes(label=d_CH), y=20)+
  theme_classic()

ggplot(data=combinegr_wide, aes(factor(d_S4), d_CH))+
  geom_boxplot()+
  stat_summary(fun.y=mean, color="red", geom="point")+
  geom_text(data=check, aes(label=round(d_CH, 2), y=d_CH+40), size=3)+
  geom_text(data=number, aes(label=d_CH), y=20)+
  geom_text(data=abc_group, aes(label=abc, y=925))+
  theme_classic()

d_S4<-c(0,1,2,3)
abc<-c("ab","a","bc","c")
abc_group<-data.frame(d_S4, abc)



#describe the trend (linear) between score and biomass (or height, or trait response)
fit<-lm(formula=combinegr_wide$d_TM~combinegr_wide$d_S4)
summary(fit)

ggplot(data=combinegr_wide, aes(x=d_S4, y=d_PE))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()

ggplot(data=combinegr_wide, aes(x=d_CH, y=d_S4))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()





#split data into size quantiles 
tableOne<-within(combinegr_wide, quartile<-as.integer(cut(d_CH, quantile(d_CH, probs=0:4/4, na.rm=TRUE), include.lowest = TRUE)))

ggplot(data=tableOne, aes(x=factor(quartile), y=d_S4))+geom_point()+geom_jitter()


#test score v response correlation just with upper most quantile(s)
png(file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180502 leaf rolling tests/scoreQuantile.png")
ggplot(data=subset(tableOne, quartile=="4"), aes(x=d_S4, y=r_TM))+geom_point()+geom_smooth(method="lm")
dev.off()






#####MIDSEASON####
#score versus midseason height and tiller number 
png(file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180502 leaf rolling tests/scoreVsrMH.png")
ggplot(data=subset(combinegr_wide, !is.na(d_S4)), aes(factor(d_S4), y=r_MH))+geom_boxplot()
dev.off()

results=aov(r_MH~factor(d_S4), data=combinegr_wide)
summary(results)


capture.output(summary(results), file="/rsync/box/Darshi work/LAB MEETINGS/BANAN 180502 leaf rolling tests/r_MHtest.csv")









#####CORRELATIONS#####

#dry architecture 
pairs(combine_wide[,c("d_S0", "d_PE", "d_CH", "d_TH", "d_TN", "d_BN"
                      , "d_BW", "d_RV")])
#architecture response 
pairs(combine_wide[,c("d_RS", "d_PE", "r_CH", "r_TH", "r_TN", "r_BN", "r_BW", "r_RV")])

pairs(combine_wide[,c("d_RS", "d_PE", "d_MH", "d_CH", "d_TM", "r_MH", "r_CH", "r_TM")])

#dry productivity 
pairs(combine_wide[,c("d_RS", "d_PE", "d_LM", "d_PM", "d_SM", "d_VM", "d_TM")])

#productivity response 
pairs(combine_wide[,c("d_RS", "d_PE", "r_LM", "r_PM", "r_SM", "r_VM", "r_TM")])

#dry productivty vs productivity response 
pairs(combine_wide[,c("d_RS", "d_PE", "d_LM", "d_PM", "d_SM", "d_VM", "d_TM", "r_LM", "r_PM", "r_SM", "r_VM", "r_TM")])

























