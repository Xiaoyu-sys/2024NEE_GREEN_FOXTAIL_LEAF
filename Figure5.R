################################################################################################################################
## A-B #
conda activate qiime2-2021.2

qiime tools import \
  --input-path 16S_isolates/isolates.fasta \
  --output-path 16S_isolates/seq.qza \
  --type 'FeatureData[Sequence]'

ime qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 16S_isolates/seq.qza \
  --o-alignment 16S_isolates/aligned-seq.qza \
  --o-masked-alignment 16S_isolates/masked-aligned-seq.qza \
  --o-tree 16S_isolates/unrooted-tree.qza \
  --o-rooted-tree 16S_isolates/rooted-tree.qza \
  --p-n-threads 70
  
qiime tools export \
  --input-path 16S_isolates/rooted-tree.qza \
  --output-path 16S_isolates
  
  
################################################################################################################################
## C ##
library(EasyAovWlxPlot)
library(ggplot2)

data <- read.delim("isolates/phyllosphere/BXL_infection/metadata2.txt", sep = "\t",  check.names = FALSE)

result = MuiKwWlx(data = data,num = 2)
result

result1 = FacetMuiPlotresultBar(data = data,num = 2,result = result,sig_show ="abc",ncol = 1 )
gsd=result1[[2]]

gsd$group=factor(gsd$group,levels = c("Infection without BXL","CK","BacSynCom","YeastSynCom","BacYeastSynCom"))
ggplot(gsd,aes(group,mean,fill=group))+
    geom_errorbar(aes(ymin=mean, ymax=(mean+SD)),width=0.5,position=position_dodge(.9))+
    geom_col(position="dodge")+ 
    theme(strip.text=element_text(size = rel(0.7),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
    theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
    labs(title="",x="",y="Lesion area of each leaf (cm2)",fill="")+
    theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black',))+
    theme(axis.text.x=element_text(size=8,angle=0,vjust=0,hjust=0.5,color = "black",face = "bold"),
          axis.text.y=element_text(size=8,color = "black",face = "bold"),legend.position = "right")+
    #theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
    geom_text(aes(group,mean,label = groups), position = position_dodge(0.9),stat="identity")+
	coord_flip()
	
################################################################################################################################
## D ##
library(tidyverse)
library(reshapes2)

df=read.table("isolates/phyllosphere/bacteria/final_isolates_taxonomy.tsv",sep="\t",row.names=1,header=T)
df=melt(df)

colors3=c("#DF7469", "#00AED1", "#6DBB60")
df$Genus=factor(df$Genus,levels=c("Bacillus","Pantoea","Methylobacterium-Methylorubrum"))

ggplot(df,aes(reorder(Feature_ID,-value),value,fill=Genus))+
    geom_bar(stat="identity")+
	facet_wrap(.~variable,1)+
    theme(strip.text=element_text(size = rel(0.7),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
    theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
    labs(title="Bacterial isolates",x="",y="Percentage inhibition (100%)",fill="")+
    theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black',))+
    theme(axis.text.x=element_text(size=8,angle=0,vjust=0,hjust=0.5,color = "black",face = "bold"),
          axis.text.y=element_text(size=8,color = "black",face = "bold"),legend.position = "right")+
    #theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
	coord_flip()





