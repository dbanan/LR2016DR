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
library(emmeans)


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
  #stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,0)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
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
  #stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,0)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
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
  #stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,2)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
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
  #stat_summary(aes(interaction(time, trt_abb), average, label=round(..y..,2)), fun.y=mean, color="black", geom="text", inherit.aes = FALSE, size=3)+
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
coefs_f2b <- coef(lm(CLP ~ roll, data = subset(comparew, trt_abb=="WD")))
coefs_f2c <- coef(lm(PAI ~ roll, data = subset(comparew, trt_abb=="WD")))

#coefs_f2fww <- coef(lm(PAI ~ CLP, data = subset(comparew, trt_abb=="WW")))
#coefs_f2fwd <- coef(lm(PAI ~ CLP, data = subset(comparew, trt_abb=="WD")))

coefs_f2f <- coef(lm(PAI ~ CLP, data = comparew))
cor.test(comparew$PAI, comparew$CLP)


f2a<-ggplot(comparew, aes(x=roll, y=inclination, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  xlim(0.3, 1.5)+ylim(0.3, 1.5)+
  xlab("leaf roll angle (midday/dawn)")+ylab("leaf inclination angle (midday/dawn)")+
  theme(legend.position="none")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
f2b<-ggplot(comparew, aes(x=roll, y=CLP, color=trt_abb))+
  geom_point(size=3)+
  geom_abline(intercept=coefs_f2b[1], slope=coefs_f2b[2], size=0.5)+
  annotate("text",x=0.8,y=0.63,label="y = -0.48x + 1.44 *", size=5)+
  annotate("text",x=0.8,y=0.57,label="r = -0.55", size=5)+
  scale_color_manual(values=c("grey", "black"))+
  xlim(0.3, 1.5)+ylim(0.3, 1.5)+
  xlab("leaf roll angle (midday/dawn)")+ylab("CLP (midday/dawn)")+
  theme(legend.position="none")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
f2c<-ggplot(comparew, aes(x=roll, y=PAI, color=trt_abb))+
  geom_point(size=3)+
  geom_abline(intercept=coefs_f2c[1], slope=coefs_f2c[2], size=0.5)+
  annotate("text",x=0.8,y=0.63,label="y = 0.21x + 0.81 *", size=5)+
  annotate("text",x=0.8,y=0.57,label="r = 0.52", size=5)+
  scale_color_manual(values=c("grey", "black"))+
  xlim(0.3, 1.5)+ylim(0.3, 1.5)+
  xlab("leaf roll angle (midday/dawn)")+ylab("PAI (midday/dawn)")+
  theme(legend.position="none")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
f2d<-ggplot(comparew, aes(x=inclination, y=CLP, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  xlim(0.3, 1.5)+ylim(0.3, 1.5)+
  xlab("leaf inclination angle (midday/dawn)")+ylab("CLP (midday/dawn)")+
  theme(legend.position="none")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
f2e<-ggplot(comparew, aes(x=inclination, y=PAI, color=trt_abb))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey", "black"))+
  xlim(0.3, 1.5)+ylim(0.3, 1.5)+
  xlab("leaf inclination angle (midday/dawn)")+ylab("PAI (midday/dawn)")+
  theme(legend.position="none")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
f2f<-ggplot(comparew, aes(x=CLP, y=PAI, color=trt_abb))+
  geom_point(size=3)+
  #geom_abline(intercept=coefs_f2fww[1], slope=coefs_f2fww[2], size=0.5, color="grey")+
  #annotate("text",x=1,y=1.35,label="y = -1.04x + 2.05 **")+
  #annotate("text",x=1,y=1.25,label="r = -0.71")+
  #geom_abline(intercept=coefs_f2fwd[1], slope=coefs_f2fwd[2], size=0.5)+
  #annotate("text",x=1,y=1.35,label="y = -0.36x + 1.35 ***")+
  #annotate("text",x=1,y=1.25,label="r = -0.80")+
  geom_abline(intercept=coefs_f2f[1], slope=coefs_f2f[2], size=0.5, linetype=2)+
  annotate("text",x=0.8,y=0.63,label="y = -0.49x + 1.52 ***", size=5)+
  annotate("text",x=0.8,y=0.57,label="r = -0.78", size=5)+
  scale_color_manual(values=c("grey", "black"))+
  xlim(0.3, 1.5)+ylim(0.3, 1.5)+
  xlab("CLP (midday/dawn)")+ylab("PAI (midday/dawn)")+
  theme(legend.position="none")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )



plot_grid(f2a, NULL, NULL, 
          f2b, f2d, NULL,
          f2c, f2e, f2f,
          ncol=3,
          labels=c("A","","","B","C","","D","E","F"))

plot_grid(f2a, f2f,
          f2b, f2d, 
          f2c, f2e, 
          ncol=2, 
          labels=c("A","B","C","D","E","F"))

png(file="./results/f2rev.png", width=660, height=990)
plot_grid(f2a, f2f,
          f2b, f2d,
          f2c, f2e, 
          ncol=2, 
          labels=c("A","B","C","D","E","F"))
dev.off()




#calculate corr values and pvalues and coefficients, input to supplemental table (by hand:|) 

compareww<-subset(comparew, trt_abb=="WW")
cor.test(compareww$PAI, compareww$CLP)

coef(lm(PAI ~ CLP, data = subset(comparew, trt_abb=="WW")))

comparewd<-subset(comparew, trt_abb=="WD")
cor.test(comparewd$PAI, comparewd$CLP)

coef(lm(PAI ~ CLP, data = subset(comparew, trt_abb=="WD")))


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

dPAI_rollwd<-subset(dPAI_roll, trt_abb=="WD")

f3a<-ggplot(dPAI_rollwd, aes(x=roll, y=dPAI_trt_per))+
  geom_point()+
  geom_abline(intercept=coefs_f3a[1], slope=coefs_f3a[2], size=0.5)+
  annotate("text",x=0.87,y=1.05,label="y = -0.35x + 1.21", size=5)+
  annotate("text",x=0.87,y=1.02,label=expression(paste(r^{2},"=0.45 **")), size=5)+
  ylab("dawn PAI (WD/WW)")+
  xlab("WD leaf roll angle (midday/dawn)")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

f3b<-ggplot(dPAI_rollwd, aes(x=CLP, y=dPAI_trt_per))+
  geom_point()+
  geom_abline(intercept=coefs_f3b[1], slope=coefs_f3b[2], size=0.5)+
  annotate("text",x=1.35,y=0.8,label="y = 0.08x + 0.81", size=5)+
  annotate("text",x=1.35,y=0.77,label=expression(paste(r^{2},"=0.27 *")), size=5)+
  ylab("dawn PAI (WD/WW)")+
  xlab("WD CLP (midday/dawn)")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )#+
  #theme(axis.title.y=element_blank(),
  #      axis.text.y=element_blank())
        

f3c<-ggplot(data=subset(dPAI_roll, trt_abb=="WD"), aes(x=LRS, y=dPAI_trt_per))+
  geom_point()+
  geom_abline(intercept=coefs_f3c[1], slope=coefs_f3c[2], size=0.5)+
  annotate("text",x=1,y=1.05,label="y = 0.08x + 0.81", size=5)+
  annotate("text",x=1,y=1.02,label=expression(paste(r^{2},"=0.48 **")), size=5)+
  ylab("dawn PAI (WD/WW)")+
  xlab("leaf rolling score")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )#+
  #theme(axis.title.y=element_blank(),
  #      axis.text.y=element_blank())
        

coefs_f3a<-coef(lm(dPAI_trt_per~roll, data=dPAI_rollwd))
cor.test(dPAI_rollwd$dPAI_trt_per, dPAI_rollwd$roll)

coefs_f3b<-coef(lm(dPAI_trt_per~CLP, data=dPAI_roll))
cor.test(dPAI_rollwd$dPAI_trt_per, dPAI_rollwd$CLP)

coefs_f3c<-coef(lm(dPAI_trt_per~LRS, data=dPAI_roll))
cor.test(dPAI_rollwd$dPAI_trt_per, dPAI_rollwd$LRS)

plot_grid(f3a, f3b, f3c, ncol=3, labels=c("A","B","C"))

png(file="./results/f3rev.png", width=900, height=500)
plot_grid(f3a, f3b, f3c, ncol=3, labels=c("A","B","C"))
dev.off()

##############
###FIGURE 4###
##############

#format whole plant trait dataset 
#select traits of interest 
interest<-c("basal_width","branch_number","culm_height","tiller_height","tiller_number",
            "per_plant_leaf_mass","per_plant_panicle_mass","per_plant_stem_mass","per_plant_vegetative_mass","per_plant_total_mass",
            "reproductive_vegetative_ratio","panicle_emerge","leaf_number_total",
            "height","tiller_count")
allplant1<-subset(allplant, trait %in% interest)

#trait abbreviations 
abb<-c("BW","BN","CH","TH","TN","LM","PM","SM","VM","TM","RV","PE","LN","MH","MT")
key<-cbind(interest, abb)
colnames(key)<-c("trait","abb")
allplant1<-merge(allplant1, key, by="trait")

#genotype averages for plant and score 
allplantg<-ddply(allplant1, c("genotype", "treatment", "trait", "abb"), summarise, average=mean(data))
visualg<-ddply(visual_score, c("genotype"), summarise, score=mean(data))

#make wide by treatment 
allplantw<-dcast(allplantg, genotype+trait+abb~treatment, value.var="average")

#calculate treatment differences 
allplantw$dw_diff<-allplantw$dry-allplantw$wet
allplantw$abs_diff<-abs(allplantw$dw_diff)
allplantw$rel_diff<-allplantw$dw_diff/allplantw$wet
allplantw$dw_percent<-allplantw$dry/allplantw$wet
allplantw$log_percent<-log(allplantw$dw_percent)
allplantw$mid_diff<-allplantw$dw_diff/((allplantw$dry+allplantw$wet)/2)

#long by data type 
allplantl<-melt(allplantw, id.vars=c("genotype","trait","abb"), measure.vars=c("dry","wet","dw_diff","abs_diff","rel_diff","dw_percent","log_percent","mid_diff"), variable.name="type",value.name="data")
#trait_type combination 
allplantl$abb_type<-paste(allplantl$abb, allplantl$type, sep="_")

#go wide again by trait+type combination 
allplantww<-dcast(allplantl, genotype~abb_type, value.var="data")

#join with score data
allplantww<-merge(allplantww, visualg, by=c("genotype"))

#score vs height (can make a fit line), poor r-squared but the line is significantly non-zero
allplantww1<-allplantww
allplantww1$score[allplantww1$score==3]<-2.5

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

allplantww2<-subset(allplantww, genotype %in% sibling)
allplantww1$sibling[allplantww1$genotype %in% sibling]<-"black"
allplantww1$sibling[is.na(allplantww1$sibling)]<-"grey"

#calculate size and mass differences between high and low rolling pairs/sets 
#pull out size and mass
allplantww3<-allplantww1[,c("genotype","score","VM_dry","VM_wet","VM_dw_percent","CH_dry","CH_wet","CH_dw_percent")]
#merge CH and VM with CT40
allplantww3<-merge(allplantww3, CT40w, by=c("genotype"))
#merge with SD data (P1, raw)
allplantww3<-merge(allplantww3, SDw, by=c("genotype"))

#merge with potential lo-rollers 
potLo<-merge(allplantww3, coph3, by="genotype")
#ditch those with missing data (NAs in CH or VM)
newdf<-potLo[!rowSums(is.na(potLo))>0,]

###keep an eye on this, maybe later try just VM?###




#special graph plant size versus leaf hydraulic conductance

hilo_pair<-hilo
colnames(hilo_pair)<-c("hi","lo","distance")
hilo_pair$pair<-paste(hilo_pair$hi, hilo_pair$lo, sep="_")
hilo_pairl<-melt(hilo_pair, id.vars=c("distance", "pair"), measure.vars=c("hi","lo"), variable.name="roll", value.name="genotype")

pheno_pair<-merge(hilo_pairl, allplantww3, by="genotype")

ggplot(pheno_pair, aes(x=SD_dry, y=VM_dry))+
  geom_point(size=5, aes(color=roll))

cor.test(pheno_pair$VM_dry, pheno_pair$SD_dry)

pheno_pairhi<-subset(pheno_pair, roll=="hi")
cor.test(pheno_pairhi$VM_dry, pheno_pairhi$SD_dry)

pheno_pairlo<-subset(pheno_pair, roll=="lo")
cor.test(pheno_pairlo$VM_dry, pheno_pairlo$SD_dry)



ggplot(allplantww3, aes(x=SD_dry, y=VM_dry))+
  geom_point()+
  geom_smooth(method="lm")




png(file="./results/SD_VM.png")
ggplot(pheno_pair, aes(x=SD_dry, y=VM_dry))+
  geom_point(size=5, aes(color=roll))
dev.off()



#trying this again with as little missing data as possible and with wet and difference data included
hilo<-read.csv("./data/clean_data/rollHiLoComparison_noMissing.csv", header=T)
#merge with phentypes data 
#los
los<-subset(allplantww3, genotype %in% hilo$rolllo)
colnames(los)<-c("rolllo","lo_score","lo_VM_dry","lo_VM_wet","lo_VM_dw_percent","lo_CH_dry","lo_CH_wet","lo_CH_dw_percent")
#colnames(los)<-c("rolllo","lo_score","lo_VM_dry","lo_VM_wet","lo_VM_dw_percent","lo_CH_dry","lo_CH_wet","lo_CH_dw_percent","lo_CT_dry","lo_CT_wet","lo_CT_dw_percent","lo_SD_dry","lo_SD_wet","lo_SD_dw_percent")
#his 
his<-subset(allplantww3, genotype %in% hilo$rollhi)
colnames(his)<-c("rollhi","hi_score","hi_VM_dry","hi_VM_wet","hi_VM_dw_percent","hi_CH_dry","hi_CH_wet","hi_CH_dw_percent")
#colnames(his)<-c("rollhi","hi_score","hi_VM_dry","hi_VM_wet","hi_VM_dw_percent","hi_CH_dry","hi_CH_wet","hi_CH_dw_percent","hi_CT_dry","hi_CT_wet","hi_CT_dw_percent","hi_SD_dry","hi_SD_wet","hi_SD_dw_percent")
#merge
hilo1<-merge(hilo, los, by="rolllo")
hilo2<-merge(hilo1, his, by="rollhi")

#calculate sibling pair differences 
hilo2$hilo_VM_dry<-hilo2$hi_VM_dry-hilo2$lo_VM_dry
hilo2$hilo_CH_dry<-hilo2$hi_CH_dry-hilo2$lo_CH_dry
hilo2$hilo_CT_dry<-hilo2$hi_CT_dry-hilo2$lo_CT_dry
hilo2$hilo_SD_dry<-hilo2$hi_SD_dry-hilo2$lo_SD_dry

hilo2$hilo_VM_wet<-hilo2$hi_VM_wet-hilo2$lo_VM_wet
hilo2$hilo_CH_wet<-hilo2$hi_CH_wet-hilo2$lo_CH_wet
hilo2$hilo_CT_wet<-hilo2$hi_CT_wet-hilo2$lo_CT_wet
hilo2$hilo_SD_wet<-hilo2$hi_SD_wet-hilo2$lo_SD_wet

hilo2$hilo_VM_dif<-hilo2$hi_VM_dw_percent-hilo2$lo_VM_dw_percent
hilo2$hilo_CH_dif<-hilo2$hi_CH_dw_percent-hilo2$lo_CH_dw_percent
hilo2$hilo_CT_dif<-hilo2$hi_CT_dw_percent-hilo2$lo_CT_dw_percent
hilo2$hilo_SD_dif<-hilo2$hi_SD_dw_percent-hilo2$lo_SD_dw_percent


#average out to take care of hi rollers with multiple close siblings 
hilo3<-ddply(hilo2, c("rollhi"), summarise, 
             hiloVMd=mean(hilo_VM_dry), 
             hiloCHd=mean(hilo_CH_dry),
             hiloCTd=mean(hilo_CT_dry),
             hiloSDd=mean(hilo_SD_dry),
             hiloVMw=mean(hilo_VM_wet), 
             hiloCHw=mean(hilo_CH_wet),
             hiloCTw=mean(hilo_CT_wet),
             hiloSDw=mean(hilo_SD_wet),
             hiloVMr=mean(hilo_VM_dif), 
             hiloCHr=mean(hilo_CH_dif),
             hiloCTr=mean(hilo_CT_dif),
             hiloSDr=mean(hilo_SD_dif),
             hilodist=mean(distance))
hilo3<-ddply(hilo2, c("rollhi"), summarise, 
             hiloVMd=mean(hilo_VM_dry), 
             hiloCHd=mean(hilo_CH_dry),
             #hiloCTd=mean(hilo_CT_dry),
             #hiloSDd=mean(hilo_SD_dry),
             hiloVMw=mean(hilo_VM_wet), 
             hiloCHw=mean(hilo_CH_wet),
             #hiloCTw=mean(hilo_CT_wet),
             #hiloSDw=mean(hilo_SD_wet),
             hiloVMr=mean(hilo_VM_dif), 
             hiloCHr=mean(hilo_CH_dif),
             #hiloCTr=mean(hilo_CT_dif),
             #hiloSDr=mean(hilo_SD_dif),
             hilodist=mean(distance))

hilo4<-hilo3[complete.cases(hilo3),]
hilo21<-hilo2[complete.cases(hilo2),]
hilo22<-ddply(hilo21, c("rollhi"), summarise, 
              lo_VM=mean(lo_VM_dry),
              hi_VM=mean(hi_VM_dry),
              lo_CH=mean(lo_CH_dry),
              hi_CH=mean(hi_CH_dry))


t.test(hilo2$lo_VM_dry, hilo2$hi_VM_dry)
boxplot(hilo2$lo_VM_dry, hilo2$hi_VM_dry, names=c("low LRS","high LRS"), main="WD vegetative mass (g) p-value=0.0589")

t.test(hilo2$lo_CH_dry, hilo2$hi_CH_dry)
boxplot(hilo2$lo_CH_dry, hilo2$hi_CH_dry, names=c("low LRS","high LRS"), main="WD culm height (mm) p-value=0.0641")

t.test(hilo4$hiloVMd)
t.test(hilo4$hiloCHd)


png(file="./results/boxplot_hiloVM.png")
boxplot(hilo2$lo_VM_dry, hilo2$hi_VM_dry, names=c("low LRS","high LRS"), main="WD vegetative mass (g) p-value=0.0589")
dev.off()

png(file="./results/boxplot_hiloCH.png")
boxplot(hilo2$lo_CH_dry, hilo2$hi_CH_dry, names=c("low LRS","high LRS"), main="WD culm height (mm) p-value=0.0641")
dev.off()  




f4d<-ggplot(hilo4, aes(x=rollhi, y=hiloVMd))+
  geom_bar(stat="identity", fill="grey")+
  geom_hline(yintercept=0)+
  #geom_text(aes(label=sprintf("%0.2f", round(hilodist, digits = 2)),vjust=ifelse(hiloVM>=0,-0.3,1.3)))+
  xlab(expression("accessions with leaf roll score">=2))+ylab("")+
  scale_y_continuous("vegetative mass difference (g)", position="right")+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_text(angle=45,vjust = 1, hjust=1))

f4c<-ggplot(hilo4, aes(x=rollhi, y=hiloCHd))+
  geom_bar(stat="identity", fill="grey")+
  geom_hline(yintercept=0)+
  #geom_text(aes(label=sprintf("%0.2f", round(hilodist, digits = 2)),vjust=ifelse(hiloVM>=0,-0.3,1.3)))+
  xlab(expression("accessions with leaf roll score">=2))+ylab("")+
  scale_y_continuous("culm height difference (mm)", position="right")+
  theme_minimal()+theme(panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_text(angle=45,vjust = 1, hjust=1))




ggplot(data=allplantww1, aes(x=score, y=CH_dry))+geom_point()+geom_smooth(method="lm")

fit<-lm(formula=allplantww1$CH_dry~allplantww1$score)
summary(fit)

coefs <- coef(lm(CH_dry ~ score, data = allplantww1))






modelLRS_CH = lm(CH_dry ~ factor(score),
                 data=allplantww1)

leastsquareLRS_CH = lsmeans(modelLRS_CH,
                            pairwise ~ factor(score),
                            adjust = "tukey")

groupingLRS_CH<-cld(leastsquareLRS_CH,
                    alpha   = 0.05,
                    Letters = letters,
                    adjust="tukey")


f4a<-ggplot(allplantww1, aes(factor(score), CH_dry))+
  geom_jitter(width=0.075, color=allplantww1$sibling)+
  #geom_abline(intercept=coefs[1], slope=coefs[2], size=0.5)+
  ylim(0,NA)+
  ylab("WD culm height (mm)")+
  xlab("")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

#LRS vs WD vegetative mass 
modelLRS_VM = lm(VM_dry ~ factor(score),
           data=allplantww1)

leastsquareLRS_VM = lsmeans(modelLRS_VM,
                      pairwise ~ factor(score),
                      adjust = "tukey")

groupingLRS_VM<-cld(leastsquareLRS_VM,
              alpha   = 0.05,
              Letters = letters,
              adjust="tukey")

f4b<-ggplot(data=allplantww1, aes(factor(score), VM_dry))+
  #geom_boxplot()+
  geom_jitter(width=0.125, color=allplantww1$sibling)+
  #stat_summary(fun.y=mean, color="blue", geom="point")+
  ylab("WD vegetative mass (g)")+xlab("leaf rolling score")+
  #geom_text(data=check, aes(label=round(VM_dry, 2), y=VM_dry))+
  #geom_text(data=number, aes(label=VM_dry, y=0))+
  geom_text(data=groupingLRS_VM, aes(label=.group, x=factor(groupingLRS_VM$score), y=6))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )


png(file="./results/f4rev.png", width=900, height=500)
plot_grid(f4a, f4c, f4b, f4d, ncol=2, labels=c("A","B","C","D"))
dev.off()

plot_grid(f4a, f4c, f4b, f4d, ncol=2, labels=c("A","B","C","D"))






##############
###FIGURE 5###
##############

#join score and canopy temperature 
LRS_CT<-merge(CTfullg, visualg, by=c("genotype"))

LRS_CT<-subset(LRS_CT, trait=="CT40DAS")

LRS_CT$sibling[LRS_CT$genotype %in% sibling]<-"black"
LRS_CT$sibling[is.na(LRS_CT$sibling)]<-"grey"

LRS_CT$score[LRS_CT$score==3]<-2.5

ggplot(LRS_CT1, aes(x=factor(score), y=average))+
  geom_jitter(width=0.125, color=LRS_CT1$sibling)
  

#do canopy temperature wide by treatment and calculate a difference and then merge with the sibling hi lo comparison 
CT40<-subset(CTfullg, trait=="CT40DAS")
CT40w<-dcast(CT40, genotype~treatment, value.var="average")
colnames(CT40w)<-c("genotype","CT40_dry","CT40_wet")
  
CT40w$CT40_dw_percent<-CT40w$CT40_dry/CT40w$CT40_wet



##############
###FIGURE 6###
##############

#compare LRS and SD 
LRS_SD<-merge(SDg, visualg, by="genotype")

LRS_SD$score[LRS_SD$score==3]<-2.5

ggplot(LRS_SD, aes(x=factor(score), y=average))+
  geom_jitter(width=0.125)

LRS_SD1<-subset(LRS_SD, treatment=="dry")
LRS_SD1<-subset(LRS_SD1, collection=="P1")
LRS_SD1<-subset(LRS_SD1, trait=="SDraw")

LRS_SD1$sibling[LRS_SD1$genotype %in% sibling]<-"black"
LRS_SD1$sibling[is.na(LRS_SD1$sibling)]<-"grey"

ggplot(LRS_SD1, aes(x=factor(score), y=average))+
  geom_jitter(width=0.125, color=LRS_SD1$sibling)

modelLRS_SD = lm(average ~ factor(score),
                 data=LRS_SD1)

leastsquareLRS_SD = lsmeans(modelLRS_SD,
                            pairwise ~ factor(score),
                            adjust = "tukey")

groupingLRS_SD<-cld(leastsquareLRS_SD,
                    alpha   = 0.05,
                    Letters = letters,
                    adjust="tukey")






#stomatal density wide by treatment and calculate difference to merge with other traits for sibling comparison 
SDg1<-subset(SDg, collection=="P1")
SDg1<-subset(SDg1, trait=="SDraw")

SDw<-dcast(SDg1, genotype~treatment, value.var="average")
colnames(SDw)<-c("genotype","SD_dry","SD_wet")
SDw$SD_dw_percent<-SDw$SD_dry/SDw$SD_wet



#try stomatal density with subset 15
SDsub<-subset(SDfull, subplot_id %in% subset15)
SDsub1<-SDsub
SDsub1$trt_abb[SDsub1$treatment=="wet"]<-"WW"
SDsub1$trt_abb[SDsub1$treatment=="dry"]<-"WD"
SDsub2<-SDsub1[,c(3,7,4)]

comparewSD<-merge(comparew, SDsub2, by=c("genotype","trt_abb"))


ggplot(comparewSD, aes(x=SDraw, y=roll, color=trt_abb))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("grey", "black"))
#:( nothing here

###################
###SUPPLEMENTALS###
###################

#S1 weather, soil moisture 

#S2 score correlation with diurnal rolling 
#remove some select points? outliers in CLP and PAI 

s2a<-ggplot(subset(dPAI_roll, trt_abb=="WD"), aes(x=roll, y=LRS))+
  geom_point()+geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  annotate("text",x=0.45,y=1,label=expression(paste(r^{2},"=0.51 **")), size=5)+
  xlab("WD leaf roll angle (midday/dawn)")+ylab("leaf rolling score")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
s2b<-ggplot(subset(dPAI_roll, trt_abb=="WD"), aes(x=inclination, y=LRS))+
  geom_point()+
  xlab("WD leaf inclination angle (midday/dawn)")+ylab("leaf rolling score")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
s2c<-ggplot(subset(dPAI_roll, trt_abb=="WD"), aes(x=CLP, y=LRS))+
  geom_point()+geom_smooth(method="lm", se=FALSE, color="black", size=0.5)+
  annotate("text",x=1.3,y=1,label=expression(paste(r^{2},"=0.47 **")), size=5)+
  xlab("WD CLP (midday/dawn)")+ylab("leaf rolling score")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
s2d<-ggplot(subset(dPAI_roll, trt_abb=="WD"), aes(x=PAI, y=LRS))+
  geom_point()+geom_smooth(method="lm", se=FALSE, color="black",size=0.5)+
  annotate("text",x=0.85,y=1,label=expression(paste(r^{2},"=0.34 *")), size=5)+
  xlab("WD PAI (midday/dawn)")+ylab("leaf rolling score")+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

plot_grid(s2a, s2b, s2c, s2d, labels=c("A","B","C","D"))

png(file="./results/s2rev.png", width=900, height=500)
plot_grid(s2a, s2b, s2c, s2d, labels=c("A","B","C","D"))
dev.off()


#S3 lack of correlation between wet and dry diurnal rolling? 

#S4 score histogram 

#S5 LRS not correlated with biomass treatment response 

#S6 dendrogram and sibling identification 







