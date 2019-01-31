#once again revise LR2016DR manuscript code 
#analyses and figures 

#figure list 
#Fig 1 reaction norms and percent change boxplots 
#Fig 2 correlation of leaf and canopy diurnal changes 
#Fig 3 dawn PAI correlation with leaf and canopy diurnal adjustment 
#Fig 4 high rolling versus low rolling sibling productivity 

#S1 weather 
#S2 LRS correlation with leaf and canopy adjustment 
#S3 ? lack of wet dry coordination ? maybe maybe not try facet_grid to produce 4x4 matrix that shows that wet and dry don't have the same type of leaf movement response to water availability
#S4 LRS histogram 
#S5 LRS lack of pattern of productivity percent change 
#S6 dendrogram of SNP distances and identification of high and low rolling siblings 

#libraries used 
library(plyr)
library(ggplot2)
library(reshape2)
library(GGally)
library(lsmeans)
library(gridExtra)
library(cowplot)

#tick marks in found here
#https://stackoverflow.com/questions/26367296/how-do-i-make-my-axis-ticks-face-inwards-in-ggplot2


#cleaned data load 
load("./data/clean_data/data_subset_LR.Rdata")
load("./data/clean_data/data_population_score.Rdata")
load("./data/clean_data/data_population_plant.Rdata")
load("./data/clean_data/data_subset_CT.Rdata")

#and go 

###########
###SETUP###
###########

#list of 15 genotype subset 
subset15<-unique(combo2$subplot_id)



##############
###FIGURE 1###
##############
#figure 1 reaction norms and adjustment boxplots 
#well-watered on left, water deficit on right 
#leaf traits (roll, inclination) on left, canopy traits (GSF, LAI) on right 
#try grey scale version 

combo2$treatment[combo2$treatment=="dry"]<-"water deficit"
combo2$treatment[combo2$treatment=="wet"]<-"well watered"
combo2$trait<-as.character(as.factor(combo2$trait))
combo2$trait[combo2$trait=="LAI"]<-"PAI"
combo2$trait[combo2$trait=="GSF"]<-"CLP" #canopy light penetration 

combo2$trt_abb[combo2$treatment=="water deficit"]<-"WD"
combo2$trt_abb[combo2$treatment=="well watered"]<-"WW"



#first vis of raw leaf and canopy data 
ggplot(data=combo2, aes(factor(time), data, fill=factor(treatment)))+
  geom_boxplot()+
  facet_wrap(~trait, scales="free")+
  theme(legend.position="none")

#look for outliers
look<-subset(combo2, trait %in% c("CLP", "PAI", "inclination", "roll"))
look1<-subset(look, genotype=="TB_0348")
#vis to check for outliers 
ggplot(data=look1, aes(factor(interaction(subsample, time, treatment)), color=(interaction(time, treatment)), data))+
  geom_point()+
  geom_line(aes(group=genotype), alpha=0.25)+
  facet_wrap(~trait, scale="free")+
  theme(legend.position="none")

#remove outliers 
combo3<-combo2
combo3$data[combo3$trait=="PAI" & combo3$data>2.5]<-NA
combo3$data[combo3$trait=="inclination" & combo3$data>70]<-NA

#calculate genotype averages from subset15 data 
combog<-ddply(combo3, c("genotype","trt_abb","time","trait"), summarise, average=mean(data, na.rm=TRUE))
#visualize genotype averages 
ggplot()+
  geom_boxplot(data=combog, aes(factor(time), average, fill=factor(trt_abb)))+
  facet_wrap(~trait, scale="free")+
  theme(legend.position="none")

#remove a genotype average outlier point 
combog$average[combog$treatment=="dry"&combog$trait=="inclination"&combog$genotype=="TB_0581"]<-NA







#order treatment (wet first, then dry)
#combog$treatment<-factor(combog$trt_abb, levels=c("WW", "WD"), ordered=TRUE)

#visualize reaction norms 
combog1<-combog[!(combog$trait=="score"),]
combog1$trt_abb<-factor(combog1$trt_abb, levels=c("WW", "WD"), ordered=TRUE)

ggplot(transform(combog1, trait=factor(trait, levels=c("roll","inclination","PAI","CLP"))), 
       aes(interaction(time, trt_abb), average, group=interaction(genotype, trt_abb)))+
  geom_point(aes(color=trt_abb))+
  geom_line(alpha=0.2)+
  #scale_color_grey()+
  scale_x_discrete(labels = c("dawn","midday","dawn","midday"))+
  scale_color_manual(values=c("grey","black"))+
  ylim(0,NA)+
  stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,2)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
  #geom_text(aes(label=ifelse(time=="midday", as.character(genotype), '')), hjust=0, size=2)+
  facet_wrap(~trait, scale="free", ncol=4)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none", axis.text.x =element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())




#plotted separately and then together
f1a<-ggplot(subset(combog1, trait=="roll"), 
       aes(interaction(time, trt_abb), average, group=interaction(genotype, trt_abb)))+
  geom_point(aes(color=trt_abb))+
  geom_line(alpha=0.2)+
  scale_color_manual(values=c("grey","black"))+
  scale_x_discrete(labels = c("dawn","midday","dawn","midday"))+
  ylim(0,NA)+
  #labs(title="leaf roll angle")+
  ylab(expression("leaf roll angle "(degree)))+
  stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,0)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none", 
        axis.text.x =element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )


f1b<-ggplot(subset(combog1, trait=="inclination"), 
       aes(interaction(time, trt_abb), average, group=interaction(genotype, trt_abb)))+
  geom_point(aes(color=trt_abb))+
  geom_line(alpha=0.2)+
  scale_color_manual(values=c("grey","black"))+
  scale_x_discrete(labels = c("dawn","midday","dawn","midday"))+
  ylim(0,NA)+
  #labs(title="leaf inclination angle")+
  ylab(expression("leaf inclination angle " (degree)))+
  stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,0)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none", 
        axis.text.x =element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )


f1c<-ggplot(subset(combog1, trait=="PAI"), 
       aes(interaction(time, trt_abb), average, group=interaction(genotype, trt_abb)))+
  geom_point(aes(color=trt_abb))+
  geom_line(alpha=0.2)+
  scale_color_manual(values=c("grey","black"))+
  scale_x_discrete(labels = c("dawn","midday","dawn","midday"))+
  ylim(0,NA)+
  #labs(title="PAI")+
  ylab(expression(paste("PAI ",(m^{2}/m^{2}))))+
  stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,2)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none", 
        axis.text.x =element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )


f1d<-ggplot(subset(combog1, trait=="CLP"), 
       aes(interaction(time, trt_abb), average, group=interaction(genotype, trt_abb)))+
  geom_point(aes(color=trt_abb))+
  geom_line(alpha=0.2)+
  scale_color_manual(values=c("grey","black"))+
  scale_x_discrete(labels = c("dawn","midday","dawn","midday"))+
  ylim(0,NA)+
  #labs(title="CLP")+
  ylab("CLP (proportion)")+
  stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,2)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none", 
        #plot.title=element_text(hjust=0),
        #legend.title=element_blank(),
        #legend.position=c(0.9,0.5),
        axis.text.x =element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )


#diurnal adjustment 
#calculate time change (inclination, roll, GSF in dry plots), wide by time 
combow_time<-dcast(combog, trait+genotype+trt_abb~time, value.var="average")
combow_time$time_diff<-combow_time$midday-combow_time$dawn
combow_time$time_rel<-combow_time$time_diff/combow_time$dawn
combow_time$time_per<-combow_time$midday/combow_time$dawn
combow_time$time_logper<-log(combow_time$time_per)
#remove score (no sense in calculating difference)
combow_time<-combow_time[!(combow_time$trait=="score"),]
#combine trait and treatment 
#combow_time$trait_trt<-paste(combow_time$trt_abb, combow_time$trait, sep="_time_")
#use percent (time_per) as difference metric going forward 
combow_time1<-combow_time[,c(1:3,8)]
#remove outliers 
combow_time2<-combow_time1
#combow_time2$time_per[combow_time2$trait=="LAI" & combow_time2$time_per>1.26]<-NA
combow_time2$time_per[combow_time2$trait=="inclination" & combow_time2$time_per>2]<-NA
combow_time2$time_per[combow_time2$trt_abb=="WW" & combow_time2$trait=="CLP" & combow_time2$time_per>1.4]<-NA
combow_time2$time_per[combow_time2$trait=="CLP" & combow_time2$time_per<0.75]<-NA

#boxplot of differences
ggplot(data=combow_time2, aes(factor(trt_abb), time_per, fill=trt_abb))+
  geom_boxplot()+
  ylim(0,1.6)+
  stat_summary(aes(trt_abb, time_per), fun.y=mean, color="black", geom="point", inherit.aes = FALSE)+
  facet_wrap(~trait, ncol=4)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none",axis.title.x=element_blank())

#order treatment 
combow_time2$trt_abb<-factor(combow_time2$trt_abb, levels=c("WW","WD"), ordered=TRUE)

#order traits and plot 
#change axis labels 
ggplot(transform(combow_time2, trait=factor(trait, levels=c("roll","inclination","PAI","CLP"))), 
       aes(factor(trt_abb), time_per, fill=trt_abb))+
  geom_boxplot()+
  #scale_fill_grey()+
  scale_fill_manual(values=c("grey","white"))+
  ylim(0,1.6)+
  stat_summary(aes(trt_abb, time_per), fun.y=mean, color="black", geom="point", inherit.aes = FALSE)+
  facet_wrap(~trait, ncol=4)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none",axis.title.x=element_blank())

#plot these separately for formatting and then grid them together
f1e<-ggplot(subset(combow_time2, trait=="roll"), aes(factor(trt_abb), time_per, fill=trt_abb))+
  geom_boxplot()+
  scale_fill_manual(values=c("grey","white"))+
  ylim(0,1.6)+
  ylab("leaf roll angle (midday/dawn)")+
  stat_summary(aes(trt_abb, time_per), fun.y=mean, color="black", geom="point", inherit.aes = FALSE)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none",axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

f1f<-ggplot(subset(combow_time2, trait=="inclination"), aes(factor(trt_abb), time_per, fill=trt_abb))+
  geom_boxplot()+
  scale_fill_manual(values=c("grey","white"))+
  ylim(0,1.6)+
  ylab("leaf inclination angle (midday/dawn)")+
  stat_summary(aes(trt_abb, time_per), fun.y=mean, color="black", geom="point", inherit.aes = FALSE)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none",axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

f1g<-ggplot(subset(combow_time2, trait=="PAI"), aes(factor(trt_abb), time_per, fill=trt_abb))+
  geom_boxplot()+
  scale_fill_manual(values=c("grey","white"))+
  ylim(0,1.6)+
  ylab("PAI (midday/dawn)")+
  stat_summary(aes(trt_abb, time_per), fun.y=mean, color="black", geom="point", inherit.aes = FALSE)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none",axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

f1h<-ggplot(subset(combow_time2, trait=="CLP"), aes(factor(trt_abb), time_per, fill=trt_abb))+
  geom_boxplot()+
  scale_fill_manual(values=c("grey","white"))+
  ylim(0,1.6)+
  ylab("CLP (midday/dawn)")+
  stat_summary(aes(trt_abb, time_per), fun.y=mean, color="black", geom="point", inherit.aes = FALSE)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(legend.position="none",
    #legend.title=element_blank(),
    axis.title.x=element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

#and now grid reaction norms with boxplots of percent difference 
plot_grid(f1a, f1b, f1c, f1d, f1e, f1f, f1g, f1h, ncol=4, labels=c("A","B","C","D","E","F","G","H"), align="v")

png(file="./results/f1rev.png", width=1200, height=500)
plot_grid(f1a, f1b, f1c, f1d, f1e, f1f, f1g, f1h, ncol=4, labels=c("A","B","C","D","E","F","G","H"), align="v")
dev.off()

#re-run statisitical analyses (t-tests)

##############
###FIGURE 2###
##############

#correlation of leaf and canopy diurnal adjustments 
comparew<-dcast(combow_time2, genotype+trt_abb~trait, value.var="time_per")
ggpairs(comparew, columns=c(3:6), 
        lower=list(continuous="smooth"), 
        upper=list(continuous=wrap("cor")),
        axisLabels="internal",
        aes(color=trt_abb), se=FALSE)

#need a way to plot this in an aesthetic way
#try each pair one by one and then grid them out in a 3x3 with blanks? 


f2a<-ggplot(comparew, aes(x=roll, y=inclination, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  theme(legend.position="none")
f2b<-ggplot(comparew, aes(x=roll, y=CLP, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  theme(legend.position="none")
f2c<-ggplot(comparew, aes(x=roll, y=PAI, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  theme(legend.position="none")
f2d<-ggplot(comparew, aes(x=inclination, y=CLP, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  theme(legend.position="none")
f2e<-ggplot(comparew, aes(x=inclination, y=PAI, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  theme(legend.position="none")
f2f<-ggplot(comparew, aes(x=CLP, y=PAI, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  theme(legend.position="none")

plot_grid(f2a, NULL, NULL, 
          f2b, f2d, NULL,
          f2c, f2e, f2f,
          ncol=3)

##############
###FIGURE 3###
##############

#correlate midseason productivity response to diurnal adjustments 
#calculate PAI treatment response, merge with dirunal adjustment 

#pull out scores 
score15<-subset(visual_score, subplot_id %in% subset15)
score15<-score15[,c(2,5)]
colnames(score15)[2]<-"LRS"

#pull out dawn PAI and calculate treatment response 
dawnPAI<-subset(combog, time=="dawn" & trait=="PAI")
dawnPAIw<-dcast(dawnPAI, genotype~trt_abb, value.var="average")
dawnPAIw$trt_per<-dawnPAIw$WD/dawnPAIw$WW
dawnPAIw<-dawnPAIw[,c(1,4)]
colnames(dawnPAIw)[2]<-"dPAI_trt_per"



#diurnal adjustment correlation with dawn PAI treatment response 
#join dawn PAI treatment response with leaf and canopy diurnal adjustment 
dPAI_roll<-merge(dawnPAIw, comparew, by="genotype")
dPAI_roll<-merge(dPAI_roll, score15, by="genotype")

f3a<-ggplot(data=subset(dPAI_roll, trt_abb=="WD"), aes(x=roll, y=dPAI_trt_per))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black")+
  ylab("dawn PAI (WD/WW)")+
  xlab("WD leaf roll angle (midday/dawn)")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

f3b<-ggplot(data=subset(dPAI_roll, trt_abb=="WD"), aes(x=CLP, y=dPAI_trt_per))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black")+
  ylab("dawn PAI (WD/WW)")+
  xlab("WD CLP (midday/dawn)")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
        

f3c<-ggplot(data=subset(dPAI_roll, trt_abb=="WD"), aes(x=LRS, y=dPAI_trt_per))+
  geom_point()+
  geom_smooth(method="lm", se=FALSE, color="black")+
  ylab("dawn PAI (WD/WW)")+
  xlab("LRS")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
        

plot_grid(f3a, f3b, f3c, ncol=3, labels=c("A","B","C"), align="h", rel_widths=c(1.2,1,1))



##############
###FIGURE 4###
##############

#relate productivity to score 
#compare siblings 

###################
###SUPPLEMENTALS###
###################

#S1 weather, soil moisture 

#S2 score correlation with diurnal rolling 




#S3 lack of correlation between wet and dry diurnal rolling? 

#S4 score histogram 

#S5 LRS not correlated with biomass treatment response 

#S6 dendrogram and sibling identification 







