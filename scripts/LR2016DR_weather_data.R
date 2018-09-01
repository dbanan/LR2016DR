#LR2016DR_weather_data

library(plyr)
library(reshape2)

#infile, analyze, visualize weather data from 2016 field experiment 

#weather data goals 
#build experiment timeline (experiment events, development, measurement sequence)
#average AM and PM during heat dome diurnal from 10min average data 
#season averages 


#infile weather data 
#infile 24 hour data
ave24header<-read.csv("./data/raw_data/24hr_weather_16sept2016.csv", skip=1, header=F, nrows=1, as.is=T)
ave24header<-read.csv("./data/raw_data/24hr_weather_16sept2016.csv", skip=1, header=F, nrows=1, as.is=T)
ave24<-read.csv("./data/raw_data/24hr_weather_16sept2016.csv", header=F, skip=4)
colnames(ave24)<-ave24header
#convert timestamp format 
ave24$timestamp<-as.POSIXct(ave24$TIMESTAMP, format="%m/%d/%Y %H:%M")
#constrain time series of interest 
ave241<-subset(ave24, timestamp>="2016-06-01" & timestamp<="2016-08-31")

#infile 10 minute data 
ave10header<-read.csv("./data/raw_data/10min_weather_16sept2016.csv", skip=1, header=F, nrows=1, as.is=T)
ave10<-read.csv("./data/raw_data/10min_weather_16sept2016.csv", header=F, skip=4)
colnames(ave10)<-ave10header
#convert timestamp format 
ave10$timestamp<-as.POSIXct(ave10$TIMESTAMP, format="%m/%d/%Y %H:%M")
#constrain time series of interest 
ave101<-subset(ave10, timestamp>="2016-07-20 00:00:00" & timestamp<="2016-07-24 23:50:00")


#diurnal average AM and PM conditions 
#Temp, RH, wind, PAR

#when were diurnal measurements taken? 
#midday
#7/22/18 
#13:35:25
#16:57:41

#dawn 
#7/23/18
#07:55:20
#08:55:03


AM<-subset(ave101, timestamp>="2016-07-23 07:50:00" & timestamp<="2016-07-23 09:00:00")
PM<-subset(ave101, timestamp>="2016-07-22 13:30:00" & timestamp<="2016-07-22 17:00:00")
AM$period<-"morning"
PM$period<-"afternoon"
AMPM<-rbind(AM, PM)

AMPMl<-melt(AMPM, id.vars=c("timestamp", "period"), 
            measure.vars=c("RH_Avg","AirTemp_Avg","PAR_APOGEE_Avg"),
            variable.name="variable",
            value.name="data")

AMPMave<-ddply(AMPMl, c("period","variable"), summarise, average=mean(data))












