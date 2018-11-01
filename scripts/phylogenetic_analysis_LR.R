#phylogenetic_analysis

library(plyr)

library(circlize)
library(ape)
library(dendextend)
library(colorspace)

install.packages('circlize')
install.packages("ape")
install.packages("dendextend")
install.packages("colorspace")


#this script will assess leaf rolling across the accession's phylogeny (?) 

#load phenotype data
load("data_population_score.Rdata")


#list of genotypes used in experiment 
genotypes<-unique(visual_score1$genotype)

#trim phenotypes to just those of interest (score4, 0-3)
pheno<-subset(visual_score1, visual_score1$trait=="score4")
pheno<-ddply(pheno, c("genotype"), summarise, score4=mean(data))

#condense score further to three classes (no rolling, some rolling, severe rolling)
pheno1<-pheno
pheno1$score4[pheno1$score4>=0.5 & pheno1$score4<=1.5]<-1
pheno1$score4[pheno1$score4>=2]<-2



#load genetic map 
#work with filtered map 
snp_012<-read.table("~/Desktop/snp_100k/snp_100k_phylo.012", header=F, sep="\t", stringsAsFactors=FALSE)
snp_012<-snp_012[-1]

snp_pos<-read.table("~/Desktop/snp_100k/snp_100k_phylo.012.pos", header=F, sep="\t", stringsAsFactors=FALSE)
snp_pos$marker<-paste(snp_pos$V1, snp_pos$V2, sep="_")
snp_pos<-as.vector(snp_pos[,3])

snp_indv<-read.table("~/Desktop/snp_100k/snp_100k_phylo.012.indv", header=F, sep="\t", stringsAsFactors=FALSE)



colnames(snp_012)<-snp_pos  
rownames(snp_012)<-snp_indv$V1

snp_012<-subset(snp_012, rownames(snp_012) %in% genotypes)

snp_pheno<-snp_012
snp_pheno$genotype<-rownames(snp_pheno)
snp_pheno<-merge(snp_pheno, pheno1, by=c("genotype"))

rownames(snp_pheno)<-snp_pheno$genotype
snp_pheno$genotype<-NULL

#does the data need to be rotated (transposed)? ...actually i dont think so, so maybe skip this
#snp_map<-t(snp_012)
#snp_map<-as.data.frame(snp_map)

#code label color as a column 
snp_pheno1<-snp_pheno
snp_pheno1$color[snp_pheno1$score4==0]<-"red"
snp_pheno1$color[snp_pheno1$score4==1]<-"orange"
snp_pheno1$color[snp_pheno1$score4==2]<-"blue"

#calculate distances and clusters 
d_snp<-dist(snp_pheno1[,1:1316])
hc_snp<-hclust(d_snp)

#simple visualization 
plot(hc_snp, hang=-1)

dend<-as.dendrogram(hc_snp)
plot(dend)


#two ways to show denodrogram 
#one way
#lab.color=matrix("black",nrow(snp_pheno))
#rownames(lab.color)=rownames(snp_pheno)
#lab.color[which(score4==0)]="red"
#lab.color[which(score4==1)]="orange"
#lab.color[which(score4==2)]="blue"

#plot(as.phylo(hc_snp), cex = 0.6, tip.color = lab.color, label.offset = .75, type="fan")

png("./results/phylo_score.png", height=900, width=900)
plot(as.phylo(hc_snp), tip.color=snp_pheno1[,c(1318)], type="fan")
dev.off()





#another way
dend=as.dendrogram(hc_snp)
labels_colors(dend)<-snp_pheno1$color[order.dendrogram(dend)]
#labels_colors(dend)=lab.color[order.dendrogram(dend)]
dend<-color_branches(dend, k=5)
plot(dend)

png("./results/phylo_score_circle.png", height=1200, width=1200)
circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
dev.off()












#save these for now 
circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)

plot(as.phylo(hc.map), type = "fan", cex = 0.6, tip.color = lab.color, label.offset = .75)

dend<-as.dendrogram(hc_snp)

dend<-rotate(dend, 1:203)
dend<-color_branches(dend, k=5)

labels_colors(dend) <-
  rainbow_hcl(3)[sort_levels_values(
    as.numeric(snp_pheno[,1317])[order.dendrogram(dend)]
  )]

labels(dend)<-as.character(snp_pheno[,1])[order.dendrogram(dend)]

plot(dend)
plot(dend, horiz =  TRUE)
legend(legend = score4, fill = rainbow_hcl(3))

#end save these for now







#archive attempts 
#output dendrograms 
png("~/Desktop/meeting_dendrogram.png")
plot(dend)
dev.off()




png(file="./results/LRS_dendro.png")
plot(as.phylo(hc_snp), tip.color=snp_pheno[,c(1317)])
dev.off()

plot(as.phylo(hc_snp), tip.color=snp_pheno[,c(1317)], type="fan")











#so far this works 
d_snp<-dist(snp_012)
hc_snp<-hclust(d_snp)

plot(hc_snp)

plot(hc_snp, hang=-1, cex=0.6)

hcd = as.dendrogram(hc_snp)                             
plot(hcd, cex=0.6)   



plot(as.phylo(hc_snp), type = "fan")


