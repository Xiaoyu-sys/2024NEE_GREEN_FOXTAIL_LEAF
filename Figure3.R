#############################################################################################################################
## A-B ##

source("VarPartDist.txt") # this script is from Jiao et al., (2020): Balance between community assembly processes mediates species coexistence in agricultural soil microbiomes across eastern China

## Dist_Matrix
otu <- read.delim('leaf_16s/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- as.data.frame (t(otu))

Dist_Matrix<-as.matrix(vegdist(otu,method="bray"))

# import Environmental factors data
env1=read.delim("leaf_16s/env.txt",row.names=1,check.name=F,sep="\t")
colnames(env1)=c("pH","NO3","NH4","AK","AP","TK","TP","TN","SOC","CN","Altitude","Longitude","Latitude",'MAT','AT','RH','MAP','AS','AR')
Env=select(env1,-CN,-Longitude,-Latitude)

# Geo: Longitude,Latitude
Geo=select(env1,Longitude,Latitude)

VPD_bac=VarPartDist(Dist_Matrix,Env,Geo,Geo_Co=TRUE,Number_Permutations=999)

VPD_fun=VarPartDist(Dist_Matrix,Env,Geo,Geo_Co=TRUE,Number_Permutations=999)

write.table(VPD_bac,"leaf_16s/VPD_bac.tsv",sep="\t",row.names=TRUE)
write.table(VPD_fun,"leaf_its/VPD_fun.tsv",sep="\t",row.names=TRUE)


#############################################################################################################################
## C ##
source("Neutral model.R")

# bacteria
otu <- read.delim('leaf_16s/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- as.data.frame (t(otu))
nm_bac=Neutral.fit(otu)

# fungi
otu <- read.delim('leaf_its/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- as.data.frame (t(otu))
nm_fun=Neutral.fit(otu)


#############################################################################################################################
## D ##
## A community table with samples as rows and taxa as columns.
library(spaa)
library(plyr)
library(tidyverse)

source("habitat niche breadth analysis.R")

otu <- read.delim('leaf_16s/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)
niche_bac=as.data.frame(Com.niche(otu));names(niche_bac)="Bacteria"


#############################################################################################################################
## E ##
load("D:/Xiaoyu Zai/PhD research project/Domestication of foxtail millet/My Research Project/data/leaf_16s/icamp/16s.iCAMP.Confidence.detail.rda")

tax<- read.delim("leaf_16s/taxonomy.tsv",row.names=1,check.names = F,header=T, stringsAsFactors = FALSE)
sp.bin=res$detail$taxabin$sp.bin

tax %>% separate(Taxon,c("kingdom","phylum","class","order","family","genus","species"),sep="[;]")->tax
tax=tax[rownames(sp.bin),];sp.bin=sp.bin[match(rownames(tax),rownames(tax)),]
sp.bin.tax=cbind(sp.bin,tax)
sp.bin.tax$bin.id.new=factor(sp.bin.tax$bin.id.new)
write.table(sp.bin.tax,"leaf_16s/icamp/sp.bin.tax.tsv",sep="\t")

## The relative importance of each bin to the assembly process
b=read.delim("leaf_16s/icamp/16s.BinContributeToProcess_EachGroup.csv",header=T, stringsAsFactors = FALSE,sep=",")
b=b[,-1:-3]
colnames(b)=gsub('[bin]','',colnames(b))
a=order(-colMeans(b[,2:ncol(b)]))
b_long <- as.data.frame(pivot_longer(b, cols = 2:ncol(b),names_to = "Bin_ID", values_to = "Relative_importance"))
b_long$Process=factor(b_long$Process,levels=c("HeS","HoS","DL","HD","DR"))
b_long$Bin_ID=c(paste("bin",b_long$Bin_ID,sep=""))

colors5=c( "#B1BBEF","#6DB7C2","#BB8C9E","#EC8C82","#FFE68E")

p1_1=ggplot(b_long,aes(x = reorder(Bin_ID,-Relative_importance), y = Relative_importance*100,fill=Process))+geom_bar(position="stack",stat="identity")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
  labs(title="",x="",y="The relative importance (%) of community assembly process")+
  theme(axis.text.x=element_text(size=10,angle=270,vjust=0.5,hjust=0,color = "black"),
        axis.text.y=element_text(size=10,color = "black"), legend.text=element_text(size=10),legend.position="right",legend.key.size=unit(0.3,"cm"))+
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+#去除刻度线
  guides(fill=guide_legend(title=NULL))+theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  scale_fill_manual(values = colors5) #+coord_polar(theta = "x")
p1_1















