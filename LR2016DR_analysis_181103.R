#LR2016DR_analysis 


library(plyr)
library(ggplot2)
library(reshape2)
library(GGally)
library(lsmeans)
library(gridExtra)


save(combo, file="./data/clean_data/data_subset_LR.Rdata")
save(visual_score1, file="./data/clean_data/data_population_score.Rdata")
save(combine, file="./data/clean_data/data_population_plant.Rdata")



load("./data/clean_data/data_subset_LR.Rdata")
load("./data/clean_data/data_population_score.Rdata")
load("./data/clean_data/data_population_plant.Rdata")
load("./data/clean_data/data_subset_CT.Rdata")


save.image(file="tooMuch.Rdata")



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
scoreVch<-ggplot(data=combineww1, aes(x=score, y=CH_dry))+
  geom_point(color=combineww1$sibling)+
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  ylim(0,NA)+
  geom_text(aes(x=2, y=810, label="y = 47.76x + 337.30"))+
  geom_text(aes(x=2, y=790, label="p-value = 0.002"))+
  geom_text(aes(x=2, y=770, label="r-squared = 0.04"))+
  theme_classic()
  #theme(panel.background = element_rect(fill = "white", colour = "black"))
scoreVch
dev.off()

coefs <- coef(lm(CH_dry ~ score, data = combineww1))
scoreVch<-ggplot(combineww1, aes(score, CH_dry, group=score))+
  geom_jitter(width=0.075, color=combineww1$sibling)+
  geom_abline(intercept=coefs[1], slope=coefs[2], size=0.5)+
  ylim(0,NA)+
  ylab("DS culm height (mm)")+
  xlab("")+
  #geom_text(aes(x=2, y=810, label="y = 47.76x + 337.30"))+
  #geom_text(aes(x=2, y=790, label="p-value = 0.002"))+
  #geom_text(aes(x=2, y=770, label="r-squared = 0.04"))+
  annotate("text",x=2,y=830,label="y = 47.76x + 337.30")+
  annotate("text",x=2,y=790,label="p-value = 0.002")+
  annotate("text",x=2,y=750,label="r-squared = 0.04")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),axis.ticks.x=element_blank(),axis.text.x=element_blank())
  #theme_classic()+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank())
scoreVch







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
scoreVvm<-ggplot(data=combineww1, aes(factor(score), VM_dry))+
  #geom_boxplot()+
  geom_jitter(width=0.125, color=combineww1$sibling)+
  #stat_summary(fun.y=mean, color="blue", geom="point")+
  ylab("DS vegetative mass (g)")+xlab("leaf rolling score")+
  #geom_text(data=check, aes(label=round(VM_dry, 2), y=VM_dry))+
  #geom_text(data=number, aes(label=VM_dry, y=0))+
  geom_text(data=grouping, aes(label=.group, x=factor(grouping$score), y=6))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
scoreVvm
dev.off()


png("./results/scoreVchvm.png", width=700, height=400)
grid.arrange(scoreVch, scoreVvm, ncol=2)
dev.off()



#try the sibling thing 
#closest zero
sibling<-c("TB_0612",
           "TB_0368",
           "TB_0462",
           "TB_0335",
           "TB_0239",
           "TB_0003",
           "TB_0446",
           "TB_0581",
           "TB_0436",
           "TB_0433",
           "TB_0386",
           "TB_0382",
           "TB_0095",
           "TB_0508",
           "TB_0476",
           "TB_0509",
           "TB_0037",
           "TB_0122",
           "TB_0208",
           "TB_0254",
           "TB_0269",
           "TB_0260",
           "TB_0129",
           "TB_0146",
           "TB_0255",
           "TB_0437",
           "TB_0344",
           "TB_0015",
           "TB_0605",
           "TB_0050",
           "TB_0469",
           "TB_0540",
           "TB_0534",
           "TB_0480"
)
#closest 0-05
sibling<-c("TB_0612","TB_0551",
            "TB_0368","TB_0125",
            "TB_0462",
            "TB_0335","TB_0130","TB_0208","TB_0254","TB_0269",
            "TB_0239","TB_0260",
            "TB_0003","TB_0129",
            "TB_0446","TB_0343",
            "TB_0581","TB_0255",
            "TB_0436",
            "TB_0433","TB_0205",
            "TB_0386","TB_0015",
            "TB_0382","TB_0610",
            "TB_0095","TB_0050",
            "TB_0508","TB_0486",
            "TB_0476","TB_0477",
            "TB_0509","TB_0513",
            "TB_0037","TB_0503")
            
#closest 0-0.5 with full data available 
sibling<-c()



combineww2<-subset(combineww, genotype %in% sibling)
combineww1$sibling[combineww1$genotype %in% sibling]<-"black"
combineww1$sibling[is.na(combineww1$sibling)]<-"grey"

png(file="./results/scoreVMdiff.png", width=700, height=400)
grid.arrange(
ggplot(data=combineww1, aes(factor(score), VM_rel_diff))+
  #geom_boxplot()+
  geom_jitter(width=0.125, color=combineww1$sibling)+
  geom_hline(yintercept=0)+
  theme(panel.background = element_rect(fill = "white", colour = "black")),
ggplot(data=combineww1, aes(factor(score), CH_rel_diff))+
  #geom_boxplot()+
  geom_jitter(width=0.125, color=combineww1$sibling)+
  geom_hline(yintercept=0)+
  theme(panel.background = element_rect(fill = "white", colour = "black")),ncol=2)
dev.off()



#calculate size and mass differences between high and low rolling pairs/sets 
#pull out size and mass
combineww3<-combineww1[,c("genotype","score","VM_dry","VM_wet","VM_dw_percent","CH_dry","CH_wet","CH_dw_percent")]
#merge with potential lo-rollers 
potLo<-merge(combineww3, coph3, by="genotype")
#ditch those with missing data (NAs in CH or VM)
newdf<-potLo[!rowSums(is.na(potLo))>0,]

write.csv(newdf, file="./results/rollhilofull.csv")

traitlook<-subset(combineww1, genotype %in% rollhi)
traitlook<-traitlook[,c("genotype","score","VM_dry","VM_wet","VM_dw_percent","CH_dry","CH_wet","CH_dw_percent")]

#re-arranged it by hand in excel to make matches 
hilo<-read.csv("./data/clean_data/rollHiLoComparison.csv", header=T)
hilo$hilo_VM<-hilo$roll_hi_VM-hilo$roll_lo_VM
hilo$hilo_CH<-hilo$roll_hi_CH-hilo$roll_lo_CH

#take the hi genotype average 
hilo1<-ddply(hilo, c("roll_hi"), summarise, hiloVM=mean(hilo_VM), hiloCH=mean(hilo_CH), hilodist=mean(distance))


#plot bars #note to self aes(x=reorder(roll_hi, hilodist))
barVM<-ggplot(hilo1, aes(x=roll_hi, y=hiloVM))+
  geom_bar(stat="identity", fill="grey")+
  geom_hline(yintercept=0)+
  #geom_text(aes(label=sprintf("%0.2f", round(hilodist, digits = 2)),vjust=ifelse(hiloVM>=0,-0.3,1.3)))+
  xlab(expression("accessions with leaf roll score">=2))+ylab("")+
  scale_y_continuous("vegetative mass difference (g)", position="right")+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
barVM

barCH<-ggplot(hilo1, aes(x=roll_hi, y=hiloCH))+
  geom_bar(stat="identity", fill="grey")+
  geom_hline(yintercept=0)+
  labs(title="severe versus no-to-low rolling sibling comparison")+
  #geom_text(aes(label=sprintf("%0.2f", round(hilodist, digits = 2)),y=-130), size=3.5, angle=45)+
  xlab("")+ylab("")+
  scale_y_continuous("culm height difference (mm)", position="right", limits=c(-150,300))+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_blank())
barCH





#put the size and mass stuff all together
png("./results/lrsChVmSiblings.png",width=700, height=400)
grid.arrange(scoreVch, barCH, scoreVvm, barVM, ncol=2)
dev.off()






#trying this again with as little missing data as possible and with wet and difference data included
hilo<-read.csv("./data/clean_data/rollHiLoComparison_noMissing.csv", header=T)
#merge with phentypes data 
#los
los<-subset(combineww3, genotype %in% hilo$rolllo)
colnames(los)<-c("rolllo","lo_score","lo_VM_dry","lo_VM_wet","lo_VM_dw_percent","lo_CH_dry","lo_CH_wet","lo_CH_dw_percent")
#his 
his<-subset(combineww3, genotype %in% hilo$rollhi)
colnames(his)<-c("rollhi","hi_score","hi_VM_dry","hi_VM_wet","hi_VM_dw_percent","hi_CH_dry","hi_CH_wet","hi_CH_dw_percent")
#merge
hilo1<-merge(hilo, los, by="rolllo")
hilo2<-merge(hilo1, his, by="rollhi")





#make this long for easier manipulation 
hilo2$hilo_VM_dry<-hilo2$hi_VM_dry-hilo2$lo_VM_dry
hilo2$hilo_CH_dry<-hilo2$hi_CH_dry-hilo2$lo_CH_dry
hilo2$hilo_VM_wet<-hilo2$hi_VM_wet-hilo2$lo_VM_wet
hilo2$hilo_CH_wet<-hilo2$hi_CH_wet-hilo2$lo_CH_wet
hilo2$hilo_VM_dif<-hilo2$hi_VM_dw_percent-hilo2$lo_VM_dw_percent
hilo2$hilo_CH_dif<-hilo2$hi_CH_dw_percent-hilo2$lo_CH_dw_percent

hilo3<-ddply(hilo2, c("rollhi"), summarise, 
             hiloVMd=mean(hilo_VM_dry), 
             hiloCHd=mean(hilo_CH_dry), 
             hiloVMw=mean(hilo_VM_wet), 
             hiloCHw=mean(hilo_CH_wet),
             hiloVMr=mean(hilo_VM_dif), 
             hiloCHr=mean(hilo_CH_dif),
             hilodist=mean(distance))


boxplot(hilo2$lo_VM_dry, hilo2$hi_VM_dry, names=c("lo","hi"), main="VM dry 0.06")
boxplot(hilo2$lo_VM_wet, hilo2$hi_VM_wet, names=c("lo","hi"), main="VM wet 0.02")
boxplot(hilo2$lo_VM_dw_percent, hilo2$hi_VM_dw_percent, names=c("lo","hi"), main="VM difference ns")

t.test(hilo2$lo_VM_dry, hilo2$hi_VM_dry)
t.test(hilo2$lo_VM_wet, hilo2$hi_VM_wet)
t.test(hilo2$lo_VM_dw_percent, hilo2$hi_VM_dw_percent)


boxplot(hilo2$lo_CH_dry, hilo2$hi_CH_dry, names=c("lo","hi"), main="CH dry 0.06")
boxplot(hilo2$lo_CH_wet, hilo2$hi_CH_wet, names=c("lo","hi"), main="CH wet 0.0009")
boxplot(hilo2$lo_CH_dw_percent, hilo2$hi_CH_dw_percent, names=c("lo","hi"), main="CH difference 0.0006")

t.test(hilo2$lo_CH_dry, hilo2$hi_CH_dry)
t.test(hilo2$lo_CH_wet, hilo2$hi_CH_wet)
t.test(hilo2$lo_CH_dw_percent, hilo2$hi_CH_dw_percent)


grid.arrange(
  boxplot(hilo2$lo_VM_dry, hilo2$hi_VM_dry, names=c("lo","hi"), top="VM dry 0.06"),
  boxplot(hilo2$lo_VM_wet, hilo2$hi_VM_wet, names=c("lo","hi"), top="VM wet 0.02"),
  boxplot(hilo2$lo_VM_dw_percent, hilo2$hi_VM_dw_percent, names=c("lo","hi"), top="VM difference ns"),
  
  boxplot(hilo2$lo_CH_dry, hilo2$hi_CH_dry, names=c("lo","hi"), top="CH dry 0.06"),
  boxplot(hilo2$lo_CH_wet, hilo2$hi_CH_wet, names=c("lo","hi"), top="CH wet 0.0009"),
  boxplot(hilo2$lo_CH_dw_percent, hilo2$hi_CH_dw_percent, names=c("lo","hi"), top="CH difference 0.0006"),
  ncol=2
)




barVM<-ggplot(hilo3, aes(x=rollhi, y=hiloVMd))+
  geom_bar(stat="identity", fill="grey")+
  geom_hline(yintercept=0)+
  #geom_text(aes(label=sprintf("%0.2f", round(hilodist, digits = 2)),vjust=ifelse(hiloVM>=0,-0.3,1.3)))+
  xlab(expression("accessions with leaf roll score">=2))+ylab("")+
  scale_y_continuous("vegetative mass difference (g)", position="right")+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
barVM

barCH<-ggplot(hilo3, aes(x=rollhi, y=hiloCHd))+
  geom_bar(stat="identity", fill="grey")+
  geom_hline(yintercept=0)+
  #geom_text(aes(label=sprintf("%0.2f", round(hilodist, digits = 2)),vjust=ifelse(hiloVM>=0,-0.3,1.3)))+
  xlab("")+ylab("")+
  scale_y_continuous("culm height difference (mm)", position="right")+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_blank())
barCH





hilol<-melt(hilo3, id.vars=c("rollhi","hilodist"), measure.vars=c("hiloVMd","hiloCHd","hiloVMw","hiloCHw","hiloVMr","hiloCHr"), variable.name="name",value.name="data")

hilol$treatment[hilol$name %in% c("hiloVMd","hiloCHd")]<-"dry"
hilol$treatment[hilol$name %in% c("hiloVMw","hiloCHw")]<-"wet"
hilol$treatment[hilol$name %in% c("hiloVMr","hiloCHr")]<-"response"
hilol$trait[hilol$name %in% c("hiloVMd","hiloVMw","hiloVMr")]<-"VM"
hilol$trait[hilol$name %in% c("hiloCHd","hiloCHw","hiloCHr")]<-"CH"




ggplot(subset(hilol, treatment %in% c("dry")), aes(x=rollhi, y=data))+
  geom_bar(stat="identity", position="dodge", fill="grey")+
  geom_hline(yintercept=0)+
  facet_wrap(~trait, scales="free")+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                       panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank(),
                     axis.text.x=element_text(angle=45,vjust = 1, hjust=1))

png("./results/lrsChVmSiblingsWD.png",width=700, height=400)
ggplot(subset(hilol, treatment %in% c("dry","wet")), aes(x=rollhi, y=data, fill=treatment))+
  geom_bar(stat="identity", position="dodge")+
  geom_hline(yintercept=0)+
  facet_wrap(~trait, scales="free")+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_text(angle=45,vjust = 1, hjust=1))
dev.off()
  
  
  #geom_text(aes(label=sprintf("%0.2f", round(hilodist, digits = 2)),vjust=ifelse(hiloVM>=0,-0.3,1.3)))+
  #xlab(expression("accessions with leaf roll score">=2))+ylab("")+
  #scale_y_continuous("vegetative mass difference (g)", position="right")+
  #theme_minimal()+theme(panel.grid.minor.y=element_blank(),
   #                     panel.grid.minor.x=element_blank(),
    #                    panel.grid.major.x=element_blank(),
     #                   axis.text.x=element_text(angle=45,vjust = 1, hjust=1))



###PULL OUT ROOT SAMPLES FOR SIBLINGS 
#print out list of hi siblings and potential lo siblings to pull root samples for processing
sibList<-read.csv("./data/clean_data/rollHiLoSiblingsList.csv", header=T)

#infile experimental design
design16dr<-read.csv("/rsync/box/Setaria/2016 Setaria/16DR_subplot_genotype_list_FINAL.csv",header=T, stringsAsFactors=F,na.strings=".")
#format column names
colnames(design16dr)<-c("subplot_id","rep","awning","treatment","genotype")

sibList1<-merge(sibList, design16dr, by=c("genotype"))

write.csv(sibList1, "./results/HiLoSiblings.csv")






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

#combo2 is the infiled data 
ggplot(data=combo2, aes(factor(time), data, fill=factor(treatment)))+
  geom_boxplot()+
  facet_wrap(~trait, scales="free")+
  theme(legend.position="none")

#look for outliers
look<-subset(combo2, trait %in% c("GSF", "LAI", "inclination", "roll"))
look1<-subset(look, genotype=="TB_0348")
ggplot(data=look1, aes(factor(interaction(subsample, time, treatment)), color=(interaction(time, treatment)), data))+
  geom_point()+
  geom_line(aes(group=genotype), alpha=0.25)+
  facet_wrap(~trait, scale="free")+
  theme(legend.position="none")

#remove very high LAI values 
combo3<-combo2
combo3$data[combo3$trait=="LAI" & combo3$data>2.5]<-NA
combo3$data[combo3$trait=="inclination" & combo3$data>70]<-NA




#calculate genotype averages from subset15 data 
combog<-ddply(combo3, c("genotype","treatment","time","trait"), summarise, average=mean(data, na.rm=TRUE))
#visualize 
ggplot()+
  geom_boxplot(data=combog, aes(factor(time), average, fill=factor(treatment)))+
  facet_wrap(~trait, scale="free")+
  theme(legend.position="none")


#remove the dry inclination outlier 
#TB_0581
combog$average[combog$treatment=="dry"&combog$trait=="inclination"&combog$genotype=="TB_0581"]<-NA



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
#add to this the mean for each trait*time*treatment combo 
png(file="./results/subset_rxn.png")
rxn_raw<-ggplot(data=subset(combog, trait %in% c("GSF","LAI","inclination","roll")), aes(interaction(time, treatment), average, group=interaction(genotype, treatment)))+
  geom_point(aes(color=treatment))+
  geom_line(alpha=0.2)+
  scale_x_discrete(labels = c("dawn","midday","dawn","midday"))+
  ylim(0,NA)+
  stat_summary(aes(interaction(time, treatment), average, label=round(..y..,2)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
  #geom_text(aes(label=ifelse(time=="midday", as.character(genotype), '')), hjust=0, size=2)+
  facet_wrap(~trait, scale="free", ncol=4)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())
rxn_raw
dev.off()

#t-tests of raw observations dawn vs midday within each treatment 

#two factor anova
do_aov<-aov(average~treatment+time+treatment:time, data=subset(combog, trait=="roll"))
summary(do_aov)
TukeyHSD(do_aov)





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

png("./results/rolling_raw_corr.png")
ggpairs(combot1, columns=(4:8), aes(shape=time))
dev.off()

ggplot(data=combot1, aes(x=roll, y=inclination))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))


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

#use percent (time_per) as difference metric going forward 
combow_time1<-combow_time[,c(1:3,8)]

#remove outliers 
combow_time2<-combow_time1
#combow_time2$time_per[combow_time2$trait=="LAI" & combow_time2$time_per>1.26]<-NA
combow_time2$time_per[combow_time2$trait=="inclination" & combow_time2$time_per>2]<-NA
combow_time2$time_per[combow_time2$treatment=="wet" & combow_time2$trait=="GSF" & combow_time2$time_per>1.4]<-NA
combow_time2$time_per[combow_time2$trait=="GSF" & combow_time2$time_per<0.75]<-NA

#boxplot of differences
box_diff<-ggplot(data=combow_time2, aes(factor(treatment), time_per, fill=treatment))+
  geom_boxplot()+
  ylim(0,1.6)+
  stat_summary(aes(treatment, time_per), fun.y=mean, color="black", geom="point", inherit.aes = FALSE)+
  facet_wrap(~trait, ncol=4)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none",axis.title.x=element_blank())
box_diff

png("./results/rxnRaw_boxDiff.png", width=700, height=400)
grid.arrange(rxn_raw, box_diff)
dev.off()

#need to do a test to show the time differences are indeed different
do_aov<-aov(time_per~treatment, data=subset(combow_time2, trait=="GSF"))
summary(do_aov)

t.test(time_per~treatment, data=subset(combow_time2, trait=="GSF"))






#plot of differences 
#barplot by genotype individually 
ggplot(data=combow_time, aes(genotype, time_logper, fill=treatment))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~trait, scales="free")
#genotypes grouped as boxplot
png(file="./results/time_delta_boxplot.png")
ggplot(data=combow_time, aes(factor(treatment), time_per, fill=treatment))+
  geom_boxplot()+
  facet_wrap(~trait, scales="free")
dev.off()
#genotypes individually as response plot 
ggplot(data=combow_time2, aes(treatment, time_per, group=genotype))+
  geom_point(aes(color=treatment))+
  geom_line(alpha=0.1)+
  facet_wrap(~trait, scale="free")





#combine diurnal time change with canopy temperature 
rolling<-combow_time2
colnames(rolling)[4]<-"data"
temp<-CT16dr1[,c(2:5)]
temp$data[temp$data>44]<-NA
compare<-rbind(rolling, temp)
#wide on trait 
comparew<-dcast(compare, genotype+treatment~trait, value.var="data")

#add score
comparew2<-merge(comparew, score15, by="genotype")
comparew2$score[comparew2$treatment=="wet"]<-0

#delta rolling wet and dry combined
ggpairs(comparew, columns=c(3:6), lower=list(continuous="smooth"))
ggpairs(comparew2, columns=c(9,3:8), lower=list(continuous="smooth"), aes(color=treatment))


#dry and wet deltas for leaf and canopy traits 
png(file="./results/delta_wd_corr.png", width=700, height=700)
ggpairs(comparew2, columns=c(3:6), 
        lower=list(continuous="smooth"), 
        upper=list(continuous=wrap("cor")),
        axisLabels="internal",
        aes(color=treatment), se=FALSE)
dev.off()

#calculate corr values and pvalues 
comparew2d<-subset(comparew2, treatment=="dry")
cor.test(comparew2d$GSF, comparew2d$roll)

comparew2w<-subset(comparew2, treatment=="wet")
cor.test(comparew2w$inclination, comparew2w$roll)





#delta rolling wet and dry with CT 
png(file="./results/delta_wd_ct_corr.png")
ggpairs(comparew2, columns=c(9,3:8), lower=list(continuous="smooth"), aes(color=treatment))
dev.off()

png(file="./results/delta_ct_corr.png")
ggpairs(comparew, columns=c(3:8), lower=list(continuous="smooth"))
dev.off()


comparewd<-subset(comparew2, treatment=="dry")
compareww<-subset(comparew2, treatment=="wet")

comparewd<-comparewd[,c(1,3:6)]
compareww<-compareww[,c(1,3:6)]

colnames(comparewd)<-c("genotype","diff_dry_GSF","diff_dry_LAI","diff_dry_inclination","diff_dry_roll")
colnames(compareww)<-c("genotype","diff_wet_GSF","diff_wet_LAI","diff_wet_inclination","diff_wet_roll")

compare_wd<-merge(comparewd, compareww, by="genotype")


png(file="./results/time_delta_wetdry_corr2.png", width=700, height=700)
ggpairs(compare_wd, columns=c(2:9), lower=list(continuous="smooth"))
dev.off()


#maybe just focus on a couple "oddities"? 
fit<-lm(formula=compare_wd$diff_wet_GSF~compare_wd$diff_dry_roll)
summary(fit)


png(file="./results/wetDryMismatch.png", width=500, height=270)
grid.arrange( 
ggplot(data=compare_wd, aes(x=diff_dry_roll, y=diff_wet_GSF))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  geom_text(label="y=0.26x+0.73 ** r^2=0.53", aes(0.7, 0.75))+
  theme(panel.background = element_rect(fill = "white", colour = "black")),
ggplot(data=compare_wd, aes(x=diff_dry_roll, y=diff_wet_inclination))+
  geom_point()+
  geom_hline(yintercept=1)+
  geom_vline(xintercept=1)+
  theme(panel.background = element_rect(fill = "white", colour = "black")),
ncol=2)
dev.off()
  
png(file="./results/wetDryMismatch.png", width=300, height=270)
ggplot(data=compare_wd, aes(x=diff_dry_roll, y=diff_wet_inclination))+
  geom_point()+
  geom_hline(yintercept=1)+
  geom_vline(xintercept=1)+
  xlab("DS roll angle midday:dawn")+
  ylab("WW inclination angle midday:dawn")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


#relate rolling time delta to canopy temperature and biomass treatment delta 
#comparew2
#calculate dawn LAI treatment difference 

laig<-subset(combog, time=="dawn" & trait=="LAI")
laig<-laig[,c(1,2,4,5)]
colnames(laig)[4]<-"data"


laiw<-dcast(laig, genotype~treatment, value.var="data")
laiw$trait<-"dawnLAI"
laiw$trt_per<-laiw$dry/laiw$wet
laiw1<-laiw[,c(1,5)]
colnames(laiw1)[2]<-"dLAI_diff_trt"

#merge LAI trt delta with rolling time delta and dry score 
compare_big<-merge(compare_wd, laiw1, by="genotype")
compare_big<-merge(compare_big, score15, by="genotype")


#rolling + dLAI 
rolling1<-subset(rolling, treatment=="dry")
rolling1<-rolling1[,c(1,2,4)]
score151<-score15
score151$trait<-"score"
colnames(score151)[2]<-"data"
rolling2<-rbind(rolling1, score151)

rollLAI<-merge(rolling2, laiw1, by="genotype")
#relate dry diurnal change to treatment response 
rollLAIcorr<-ggplot(data=subset(rollLAI, trait %in% c("GSF","roll","score")), aes(x=data, y=dLAI_diff_trt))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black")+
  facet_wrap(~trait, scales="free_x")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
rollLAIcorr

#regressions for dry diurnal change to treatment response 
#embed these stats in figure 4
fit<-lm(data=rollLAI, subset=(trait=="score"), dLAI_diff_trt~data)
summary(fit)




#relate wet diurnal change to wet productiivity 
rollingw<-subset(rolling, treatment=="wet")
laiww<-laiw[,c(1,3)]
rollinglaiw<-merge(rollingw, laiww, by="genotype")
png("./results/diurnalW_prodW.png")
ggplot(rollinglaiw, aes(x=data, y=wet))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trait, scales="free_x")
dev.off()
fit<-lm(data=rollinglaiw, subset=(trait=="LAI"), wet~data)
summary(fit)

#relate wet diurnal change to treatment response 
ggplot(data=subset(rollLAI, trait %in% c("LAI","roll","inclination")), aes(x=data, y=dLAI_diff_trt))+
  geom_point()+
  #geom_smooth(method="lm", se=FALSE, color="black")+
  facet_wrap(~trait, scales="free_x")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

rollLAIwd<-merge(rolling, laiw1, by="genotype")

png("./results/rollWDvLAId.png", width=900, height=900)
ggplot(data=subset(rollLAIwd, trait %in% c("GSF","LAI","roll","inclination")), aes(data, dLAI_diff_trt, group=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(aes(color=treatment), method="lm", se=FALSE)+
  facet_wrap(~trait, scales="free_x")
dev.off()

fit<-lm(formula=compare_wd$diff_wet_GSF~compare_wd$diff_dry_roll)
summary(fit)

rollLAIw<-subset(rollLAIwd, treatment=="dry")
rollLAIw1<-subset(rollLAIw, trait=="roll")

fit<-lm(formula=rollLAIw1$dLAI_diff_trt~rollLAIw1$data)
summary(fit)



#rolling + CT
justroll<-subset(rolling, trait=="roll")
justroll<-justroll[,c(2:4)]
colnames(justroll)[3]<-"time_diff_roll"
CT16dr3<-CT16dr2
colnames(CT16dr3)<-c("treatment","genotype","trait","CT")
rollCT<-merge(justroll, CT16dr3, by=c("genotype", "treatment"))

rollCTcorr<-ggplot(data=rollCT, aes(x=time_diff_roll, y=CT, color=treatment))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE)+
  geom_smooth(method="lm", se=FALSE, inherit.aes=FALSE, linetype="dashed", aes(x=time_diff_roll, y=CT), color="black")+
  facet_wrap(~trait)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), legend.position="none")
rollCTcorr

#CT v rolling correlations 
fit<-lm(data=rollCT, subset=(trait=="CT40DAS"&treatment=="dry"), CT~time_diff_roll)
summary(fit)



png("./results/roll_CT_LAI.png", width=700, height=400)
grid.arrange(rollLAIcorr, rollCTcorr)
dev.off()



#score capture time delta rolling 
rollingScore<-merge(rolling1, score15, by="genotype")
png("./results/scoreRolling.png", width=700, height=400)
ggplot(data=rollingScore, aes(x=data, y=score))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  facet_wrap(~trait, scales="free_x")
dev.off()


#correlation values for LRS vs adjustments
fit<-lm(data=rollingScore, subset=(trait=="inclination"), score~data)
summary(fit)



#some outliers to remove from GSF and LAI...take these out of full dataset next times
rollingScore<-merge(rolling1, score15, by="genotype")
rollingScore<-rollingScore[!(rollingScore$data>1.40 & rollingScore$score==0 & rollingScore$trait=="GSF"),]
rollingScore<-rollingScore[!(rollingScore$data<0.9 & rollingScore$score==0 & rollingScore$trait=="LAI"),]
  




#delta rolling vs dawn inclination 
dawninc<-subset(combow_time, trait=="inclination")
dawninc<-dawninc[,c(2,3,4)]

dryroll<-subset(combow_time, trait=="roll")
dryroll<-dryroll[,c(2,3,8)]

rollinc<-merge(dryroll, dawninc, by=c("genotype","treatment"))

ggplot(rollinc, aes(x=time_per, y=dawn))+
  geom_point(aes(color=treatment))+
  geom_smooth(method="lm",aes(group=treatment))

ggplot(subset(rollinc, treatment=="dry"), aes(x=time_per, y=dawn))+
  geom_point()+
  geom_smooth(method="lm")

rollinc1<-subset(rollinc, treatment=="dry")
summary(lm(formula=rollinc1$dawn~rollinc1$time_per))



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

combow_time3<-combow_time2
colnames(combow_time3)[4]<-"data"
colnames(combow_time3)[3]<-"type"
#stack abs, time_diff, trt_diff 
subset_stack<-rbind(abs_all, combow_time3, prodw_trt2)
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

#some outliers to consider taking out 
comparew1<-comparew
comparew1$LAI[comparew1$LAI>1.26]<-NA
comparew1$inclination[comparew1$inclination>2]<-NA
comparew1$GSF[comparew1$treatment=="wet"&comparew1$GSF>1.4]<-NA
#add score
comparew2<-merge(comparew1, score15, by="genotype")
comparew2$score[comparew2$treatment=="wet"]<-0

#delta rolling wet and dry combined
ggpairs(comparew, columns=c(3:6), lower=list(continuous="smooth"))
ggpairs(comparew2, columns=c(9,3:8), lower=list(continuous="smooth"), aes(color=treatment))


#delta rolling wet and dry with CT 
png(file="./results/delta_wd_ct_corr.png")
ggpairs(comparew2, columns=c(9,3:8), lower=list(continuous="smooth"), aes(color=treatment))
dev.off()

png(file="./results/delta_ct_corr.png")
ggpairs(comparew, columns=c(3:8), lower=list(continuous="smooth"))
dev.off()


comparewd<-subset(comparew2, treatment=="dry")
compareww<-subset(comparew2, treatment=="wet")

comparewd<-comparewd[,c(1,3:6)]
compareww<-compareww[,c(1,3:6)]

colnames(comparewd)<-c("genotype","diff_dry_GSF","diff_dry_LAI","diff_dry_inclination","diff_dry_roll")
colnames(compareww)<-c("genotype","diff_wet_GSF","diff_wet_LAI","diff_wet_inclination","diff_wet_roll")

compare_wd<-merge(comparewd, compareww, by="genotype")


png(file="./results/time_delta_wetdry_corr2.png")
ggpairs(compare_wd, columns=c(2:9), lower=list(continuous="smooth"))
dev.off()




specialw1$trt_diff_CT30DAS[specialw1$trt_diff_CT30DAS>1.3]<-NA
specialw1$trt_diff_CT40DAS[specialw1$trt_diff_CT40DAS>1.3]<-NA

#teratment diff CT vs rolling and productivity 
png(file="./results/CTtrt_roll.png")
ggpairs(specialw1, columns=c(19:21,11:14), lower=list(continuous="smooth"))
dev.off()

png(file="./results/dryroll_wetincline.png")
ggplot(specialw1, aes(x=time_diff_dry_roll, y=time_diff_wet_inclination))+
  geom_point()+
  geom_hline(yintercept=1)+
  geom_vline(xintercept=1)
dev.off()


#roll change and canopy temp are correlated 
fitct<-lm(formula=comparew$CT40DAS~comparew$roll)
summary(fitct)
png(file="./results/roll_v_CT30.png")
ggplot(data=comparew, aes(x=roll, y=CT30DAS))+
  geom_point(aes(color=treatment), size=2)+
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  geom_text(label="y=-5.87x+40.10", aes(x=0.3, y=34.5))+
  geom_text(label="***", aes(x=0.3, y=34))+
  geom_text(label="r-squared=0.45", aes(x=0.3, y=33.5))+
  theme_classic()
dev.off()

png(file="./results/roll_v_CT40.png")
ggplot(data=comparew, aes(x=roll, y=CT40DAS))+
  geom_point(aes(color=treatment), size=2)+
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  geom_text(label="y=-2.65x+38.65", aes(x=0.3, y=36))+
  geom_text(label="***", aes(x=0.3, y=35.5))+
  geom_text(label="r-squared=0.17", aes(x=0.3, y=35))+
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
png(file="./results/time_delta_wet_corr.png")
ggpairs(specialw, columns=c(15:18), lower=list(continuous="smooth"))
dev.off()

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




#check if raw values of GSF and LAI agree with each other 
combo2w<-dcast(combo2, subplot_id+time+subsample~trait, value.var="data")

ggplot(data=combo2w, aes(x=GSF, y=LAI))+geom_point()

