coreMicro=function(otu,tax){
  library(vegan)
  library(tidyverse)
  
  otu_r<- otu/rowSums(otu)# to calculate the relative abundance
  otu_r=na.omit(otu_r)
  
  tax %>% separate(Taxon,c("domain","phylum","class","order","family","genus","species"),sep="[;]")->tax 
  tax<-subset(tax,rownames(tax) %in% colnames(otu_r))
  tax<-tax[match(colnames(otu_r),rownames(tax)),]
  
  g<- aggregate(t(otu_r), by=list(tax$genus), sum);g=g[!(g$Group.1==""),] # to obtain the genus table
  rownames(g)<-g[,1];g<-g[,-1]
  
  #Genera with relative abundance more than 0.01% present in 100% samples
  g_r_filter1 <- data.frame(g[which(apply(g, 1, function(x){length(which(x>0.001))/length(x)}) >= 1),])# g is a data frame with taxa as rows and samples as columns  
  g_r_filter1=as.data.frame (t(g_r_filter1))
  rownames(g_r_filter1)=gsub('[.]','-',rownames(g_r_filter1))
  results=list(a=g,b=g_r_filter1)
  return(results)
}

# to obtain core bacterial genera
otu <- as.data.frame(t(read.delim("leaf_16s/ASV_table_CSS.tsv", row.names=1,sep="\t",check.names = F, stringsAsFactors = FALSE)))#samples as rows and taxa as columns
tax<- read.delim("leaf_16s/taxonomy.tsv",row.names=1,check.names = F,header=T, stringsAsFactors = FALSE)
core_bac=coreMicro(otu,tax)

otu <- as.data.frame(t(read.delim("leaf_its/ASV_table_CSS.tsv", row.names=1,sep="\t",check.names = F, stringsAsFactors = FALSE)))
tax<- read.delim("leaf_its/taxonomy.tsv",row.names=1,check.names = F,header=T, stringsAsFactors = FALSE)
core_fun=coreMicro(otu,tax)

# Display as a histgram
tax_histgram<-function(tax,group){#tax is a data frame with sample names in rows and groups in columns, and group is a metadata file with sample names in rows and groups in columns.
  library(plyr)
  library(rio)
  library(ggpubr)
  library(reshape2)
  library(ggplot2)
  
  tax=tax[match(rownames(group),rownames(tax)),]
  tax$sample <- rownames(tax);group$sample <- rownames(group)
  tax<-melt(tax,ID="sample")
  tax<-merge(tax,group,by="sample",all=TRUE)
  tax$group2 <- factor(tax$group2)
  
  source("summarySE.R")
  gsd<-summarySE(tax, measurevar="value",groupvars=c("group2","variable"))
  gsd$group2<-factor(gsd$group2)
  idx=order(gsd$value,decreasing = F)
  gsd=gsd[idx,]
  
  Wilcoxon<-compare_means(value ~ group2,tax,group.by =c("variable") ,p.adjust.method = "fdr")
  Kruskal<-compare_means(value ~ group2, tax,group.by =c("variable"),method = "kruskal.test",p.adjust.method = "fdr")
  
  results<-list(a=Wilcoxon,b=Kruskal,c=gsd)
  return(results)
}

group=read.delim("leaf_16s/metadata.tsv",row.names = 1,check.names = F)

result1=tax_histgram(core_bac$b,group)
result2=tax_histgram(core_fun$b,group)

gsd1=result1$c;gsd1$variable=factor(gsd1$variable,levels=c("g__Bacillus","g__Massilia","g__Pseudomonas","g__Pantoea","g__Methylobacterium-Methylorubrum","g__Arthrobacter","g__Curtobacterium"))
p1<-ggplot(gsd1,aes(variable,value,fill=group2))+geom_bar(position="stack",stat="identity",width=0.5)+scale_fill_manual(values=colors7)+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
  labs(title="Bacteria",x="",y="Relative abundance")+
  theme(axis.text.x=element_text(size=10,angle=0,vjust=0,hjust=0.5,color = "black"),
        axis.text.y=element_text(size=10,color = "black"), legend.text=element_text(size=10),legend.key.size=unit(0.2,"inches"))+
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+#去除刻度线
  guides(fill=guide_legend(title="Sites"))+theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+coord_flip()
p1

gsd2=result2$c;gsd2$variable=factor(gsd2$variable,levels=c("g__Vishniacozyma","g__Mycosphaerella","g__Alternaria","g__Filobasidium","g__Nothophoma","g__Sporobolomyces","g__Cladosporium","g__Epicoccum","g__Neosetophoma","g__Septoria","g__Selenophoma","g__Phaeosphaeria","g__Edenia"))

p2<-ggplot(gsd2,aes(variable,value,fill=group2))+geom_bar(position="stack",stat="identity",width=0.9)+scale_fill_manual(values=colors7)+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
  labs(title="Fungi",x="",y="Relative abundance")+
  theme(axis.text.x=element_text(size=10,angle=0,vjust=0,hjust=0.5,color = "black"),
        axis.text.y=element_text(size=10,color = "black"), legend.text=element_text(size=10),legend.key.size=unit(0.2,"inches"))+
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+#去除刻度线
  guides(fill=guide_legend(title="Sites"))+theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+coord_flip()
p2

plot_grid(p1,p2, ncol =1, align = "vh",labels = "AUTO",label_size = 12)

## the correlations of core bacterial genera and fungal genera
library(pheatmap)
library(psych)
library(ggplotify)
library(cowplot)

core_bac$b=core_bac$b[match(rownames(core_fun$b),rownames(core_bac$b)),]

colnames(core_bac$b)=c(paste("B",colnames(core_bac$b),sep=": "))
colnames(core_fun$b)=c(paste("F",colnames(core_fun$b),sep=": "))
df=cbind(core_bac$b,core_fun$b)
colnames(df)=gsub('[g__]','',colnames(df))


phheat<-function(comm1,comm2){

library(pheatmap)
library(psych)
library(ggplotify)
res <- corr.test(comm1,comm2,method = "spearman",adjust = "fdr")

p=pheatmap(res$r, fontsize_number=10,fontsize = 10,cluster_rows = FALSE,
display_numbers = matrix(ifelse(res$p.adj <= 0.01, "**",
ifelse(res$p.adj <= 0.05 ,"*"," ")), nrow(res$p.adj)),
cluster_cols = FALSE,main="",cellheight=10,cellwidth = 10)
g=as.ggplot(p)

return(g)
}

phheat(df,df)

##Linear regression between specific genera
df=as.data.frame(t(read.delim("leaf_16s/core_bac_fun_genera.tsv",sep="\t",row.names=1,check.names=F)))
df_log=log10(df)

p1=ggplot(df_log,aes(x=`B: Bacillus`, y=`F: Alternaria`) )+ geom_point(shape=19)+geom_smooth(method = 'lm',level=0.95)+
labs(title = NULL, x = "Bacillus: log10 (relative abundance)", y = 'Alternaria: log10 (relative abundance)')
summary(lm(df_log$`B: Bacillus`~df_log$`F: Alternaria`,data=df_log))

p1=ggdraw(add_sub(p1,"y=1.5677x+0.3586, R2=0.3064, p < 0.001", vpadding=grid::unit(0, "lines"),y = 10, x = 0.3, hjust = 0,size=10))
p1


p2=ggplot(df_log,aes(x=`B: Methylobacterium-Methylorubrum`, y=`F: Alternaria`) )+ geom_point(shape=19)+geom_smooth(method = 'lm',level=0.95)+
labs(title = NULL, x = "Methylobacterium: log10 (relative abundance)", y = 'Alternaria: log10 (relative abundance)')
summary(lm(df_log$`B: Methylobacterium-Methylorubrum`~df_log$`F: Alternaria`,data=df_log))

p2=ggdraw(add_sub(p2,"y=-0.5472x+0.1789, R2=0.1692, p < 0.01", vpadding=grid::unit(0, "lines"),y = 10, x = 0.3, hjust = 0,size=10))
p2

p3=ggplot(df_log,aes(x=`F: Vishniacozyma`, y=`F: Alternaria`) )+ geom_point(shape=19)+geom_smooth(method = 'lm',level=0.95)+
labs(title = NULL, x = "Vishniacozyma: log10 (relative abundance)", y = 'Alternaria: log10 (relative abundance)')
summary(lm(df_log$`F: Vishniacozyma`~df_log$`F: Alternaria`,data=df_log))

p3=ggdraw(add_sub(p3,"y=-0.7622x+0.2606, R2=0.1556, p < 0.01", vpadding=grid::unit(0, "lines"),y = 10, x = 0.3, hjust = 0,size=10))
p3

plot_grid(p1,p2,p3)


