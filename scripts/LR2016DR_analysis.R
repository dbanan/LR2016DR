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
#


#subset population level 




#combine population level data 

#correlations of population plant level traits 

#model: leaf_rolling = genotype * water * time (where time is seasonal)
#^but does this really make any sense? only have leaf rolling score for dry plots 
#two biological replicates for most traits 
#have ~207 genotypes 


#15 genotype subset leaf/canopy time and treatment dynamics aka interactions 

#model: leaf_rolling = genotype * water * time (where time is diurnal) 
#use psudeo replication? 
#one true biological replicate (treatment applied to only one true biological replicate: a single 6*5 1cm square grid. drought by water withholding and well-watered by )
#have multiple measurements (subsamples) for some traits (roll angle, inclination angle)
#have 15 geno 











