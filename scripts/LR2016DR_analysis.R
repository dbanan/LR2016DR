#LR2016DR_analysis 


save(combo, file="./data/clean_data/data_subset_LR.Rdata")
save(visual_score1, file="./data/clean_data/data_population_score.Rdata")
save(combine, file="./data/clean_data/data_population_plant.Rdata")



load("./data/clean_data/data_subset_LR.Rdata")
load("./data/clean_data/data_population_score.Rdata")
load("./data/clean_data/data_population_plant.Rdata")

#analyze in order of population-plant/canopy-level to subset-leaf/canopy-level 
#communicate in reverse order (subset-leaf/canopy to population-plant/canopy)
#


#subset population-plant level data to just traits of interest 
interest<-c("basal_width","branch_number","culm_height","tiller_height","tiller_number",
            "per_plant_leaf_mass","per_plant_panicle_mass","per_plant_stem_mass","per_plant_vegetative_mass","per_plant_total_mass",
            "reproductive_vegetative_mass_ratio","panicle_emerge","leaf_number_total",
            "height","tiller_count")
combine1<-subset(combine, trait %in% interest)
#make wide by treatment (to calculate treatment differences and identify(?) absolutes)














#join plant data and score data
plant<-rbind(visual_score, combine1)


#genotype averages 
plantAve<-ddply(plant, c("genotype","treatment","trait"), summarise, average=mean(data))

#wide by trait for pairs plots 
plantw<-dcast(plantAve )
subsetw<-dcast(combine_subset, genotype+abb+trait~treatment, value.var="data")





#subset population level 
#subplot_id's of interest 
subset15<-unique(combo$subplot_id)







#correlations of population plant level traits 
#is leaf rolling limited to any specific architectural or developmental "ideotype"? 





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


subset15<-unique(combo$subplot_id)








