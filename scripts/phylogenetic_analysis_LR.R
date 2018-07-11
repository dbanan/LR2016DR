#phylogenetic_analysis


library(circlize)
library(ape)
library(dendextend)

install.packages("circlize")
install.packages("ape")
install.packages("dendextend")


#this script will assess leaf rolling across the accession's phylogeny (?) 

#load phenotype data
load("data_population_score.Rdata")


#list of genotypes used in experiment 
genotypes<-unique(visual_score1$genotype)

#trim phenotypes to just those of interest 
pheno<-subset(visual_score1, visual_score1$trait=="score4")
pheno<-ddply(pheno, c("genotype"), summarise, score4=max(data))



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
snp_pheno<-merge(snp_pheno, pheno, by=c("genotype"))


#does the data need to be rotated (transposed)? ...actually i dont think so, so maybe skip this
#snp_map<-t(snp_012)
#snp_map<-as.data.frame(snp_map)

d_snp<-dist(snp_pheno[,2:1317])
hc_snp<-hclust(d_snp)


plot(hc_snp)

dend<-as.dendrogram(hc_snp)

plot(dend)

colors_to_use<-as.numeric(snp_pheno[,1318])
colors_to_use <- colors_to_use[order.dendrogram(dend)]


labels_colors(dend)

labels_colors(dend)<-colors_to_use

#so far this works 
d_snp<-dist(snp_012)
hc_snp<-hclust(d_snp)

plot(hc_snp)

plot(hc_snp, hang=-1, cex=0.6)

hcd = as.dendrogram(hc_snp)                             
plot(hcd, cex=0.6)   



plot(as.phylo(hc_snp), type = "fan")


