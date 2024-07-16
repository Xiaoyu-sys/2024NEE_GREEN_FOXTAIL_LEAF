#Set the working directroy  to your desired directory 
setwd("D:/demo")

## A ## map

library(readr) # import data
library(rio) # import data
library(ggplot2 # Visualization
        
        encodings <- readr::guess_encoding("data/china_mapdata.csv") # to confirm the encoding types of your data
        chinamap <- read_csv("data/china_mapdata.csv", locale = locale(encoding = encodings$encoding[1])) # import your data according to the encoding types
        chinamap$prov_name <- factor(chinamap$prov_name)
        
        
        df <- subset(chinamap, prov_name == "河南省" |prov_name == "山西省" |prov_name == "陕西省" |prov_name == "河北省"|prov_name == "北京市") # to select the province corresponding to your sampling location
        
        
        p <- ggplot(df, aes(x = prov_long, y = prov_lat,  
                            group = group))+ geom_polygon(fill="white" ) +
          geom_path(colour = "darkgrey") +  
          geom_point(x=112.432521,y=37.68130,color="#BEB8DC",size=5)+ 
          annotate("text",x=112.432521,y=37.38130,label="Taiyuan",size=5,color="#BEB8DC")
        p
        
        p <- p+geom_point(x=108.019417,y=34.25004,color="#E7DAD2",size=5)+ 
          annotate("text",x=108.000,y=34.00004,label="Yangling",size=5,color="#E7DAD2")
        
        p <- p+geom_point(x=113.534877,y=34.825208,color="#96C37D",size=5)+ 
          annotate("text",x=113.334877,y=34.625208,label="Zhengzhou",size=5,color="#96C37D")
        
        p <- p+geom_point(x=116.277502,y=40.024059,color="#8ECFC9",size=5)+ 
          annotate("text",x=116.1785,y=39.77369,label="Beijing",size=5,color="#8ECFC9")
        
        p <- p+geom_point(x=114.689282,y=37.891583,color="#FFBE7A",size=5)+ 
          annotate("text",x=114.989282,y=37.51583,label="Luancheng",size=5,color="#FFBE7A")
        
        p <- p+geom_point(x=114.529025,y=38.026944,color="#FA7F6F",size=5)+ 
          annotate("text",x=114.529025,y=38.026944,label="Shijiazhuang",size=5,color="#FA7F6F")
        
        p <- p+geom_point(x=114.249027,y=37.904325,color="#82B0D2",size=5)+ 
          annotate("text",x=114.049027,y=38.054325,label="Taihang",size=5,color="#82B0D2")
        
        p+theme(axis.text.x=element_text(size=15,color = "black"),
                axis.text.y=element_text(size=15,color = "black"))
        
#############################################################################################################################
## B and F ## Phylogeny based on SNP data; the correlation of microbes and plant genotypic SNP distance via mantel test
        
        ## 1. To conduct mantel test, a correlation based analysis between two data sets of interest,
        ## i.e., microbial abundance distance matrix and host plant genotypic SNP distance matrix in this case; and
        ## 2. Plot a Tanglegram, a side-by side dendrogram visualization of two data sets to reveal 
        ## whether the microbiome profile and host genotypic distance are congruently aligned 
        ## or vary based on the distance matrix data provided respectively.
        ## this R script is from Sewunet Abera, and i modified it for my data
        
        #Set the working directroy  to your desired directory 
        setwd("D:/demo")
        
        #check if the session is routed to the set workind directory
        getwd()
        
        #check your if files and directories are present and accessible in the working directory
        list.files()
        
        #install the necessary packages for the analysis uisng biocmanager
        BiocManager::install(c("vegan", "ape", "ggplot2", "CCA",
                               "tidyverse", "cluster", "factoextra", "dendextend"))
        
        
        ##If you want to do the installation one by one follow the examples below:
        #BiocManager::install("vegan")
        #BiocManager::install("ape")
        #...........
        
        
        #load the packages
        library(vegan) #used for conducting diversity analysis, multivariate tests, ordination and distance based dissimilarity
        library(ape) #used for analyzing, manipulating and simulating phylogenetic trees 
        library(ggplot2) #a general tool for R based data visualization 
        library(CCA) #used for canonical correspondence analysis
        library(tidyverse)  #used for data manipulation
        library(cluster)    #consists clustering algorithms like k-means with other packages in this workflow depend on
        library(factoextra) #used for determining adequate number of clusters and visualization of clusters using trees
        library(dendextend) #used for visualizing and comparing two dendrograms (the tangelegram function belongs to this function)  
        
        #set the number of list R prints to max
        options(max.print=999999) 
        
        ##Conduct Mantel test for the host genotypic and microbiome data sets
        ##First import the data sets and do the necessary formatting to match the mantel test function requirements  
        
        
        ###DArT Data Import, Processing, Analysis and Visualizations
        
        #Read the dart data matrix
        dart.matrix.data <- read.csv("data/snp_dis.csv", 
                                     header = T, 
                                     row.names = 1) #in case you have the raw SNP data frame, first you have to process that and generate a matrix 
        
        
        
        #explore the dart data
        head(dart.matrix.data) 
        rownames(dart.matrix.data)
        dim(dart.matrix.data)
        class(dart.matrix.data)
        str(dart.matrix.data)
        #View(dart.data)
        
        
        #make sure the dart data is in a matrix format
        dart.matrix.data <- as.matrix(dart.matrix.data)
        
        ##explore it 
        head(dart.matrix.data)
        dim(dart.matrix.data)
        class(str(dart.matrix.data))
        str(str(dart.matrix.data))
        
        
        ##Basic K-Mean and Hierarchical Clustering of the DArT data set
        
        
        ###K-mean clustering
        
        #Normalize the dart data matrix first
        
        ##omit na values if any
        dart.matx.norm <- na.omit(dart.matrix.data)
        head(dart.matx.norm)
        
        #Normalize
        dart.matx.scaled <- scale(dart.matx.norm)
        head(dart.matx.scaled)
        
        
        #compute the euclidean distance for the dart data
        dart.matx.eucl <- dist(dart.matx.scaled, 
                               method = "euclidean")
        
        #print it to see the distance table
        dart.matx.eucl
        
        
        #check the distance based correlation among the genotypes
        fviz_dist(dart.matx.eucl)
        
        
        #do k-means clustering
        
        ##first find the optimal culster number using Gap statistic
        dart.matx.gp.stat <- fviz_nbclust(dart.matx.scaled,
                                          kmeans,  
                                          method = "gap_stat", 
                                          nboot = 999)+
          geom_vline(xintercept = 5, linetype = 2) #xintercept and linetype are personal preferences based on the plot generated
        
        #plot the result:Gap statistic (k)数值最大时对应的k即为最佳
        dart.matx.gp.stat
        
        #based on the gapstat result compute the k-mean clustering
        ###the optimal cluster number was k = 5 in my case
        
        dart.matx.km.res <- kmeans(dart.matx.scaled,
                                   8, 
                                   nstart = 25)
        
        #Print the results
        print(dart.matx.km.res)
        
        #Hierarchical clustering
        ##use the dart euclidean distance matrix to do clustering
        dart.matx.hcl <- hclust(d=dart.matx.eucl, 
                                method = "ward.D2")
        
        
        #explore teh output
        dart.matx.hcl
        
        #Make a pliminary plot
        fviz_dend(dart.matx.hcl, 
                  cex = 0.5)
        
        
        #Compute cophentic distance
        cophentic.distance <- cophenetic(dart.matx.hcl)
        
        #explore the results
        cophentic.distance
        
        
        #Make a correlation between cophenetic distance and the original distance
        cor(dart.matx.eucl, cophentic.distance)
        
        ##examine the resulst and >80% (>0.8) correlation is a good one
        #[1] 0.8554528
        
        
        #Cut tree into 5 groups (based on gap stat result above)
        ##then examine the identity of genotypes and number of genotypes in each groups 
        groups <- cutree(dart.matx.hcl, k = 5)
        
        ##explore the results
        groups
        
        #Or print the number of members in each cluster in tabular fromat
        table(groups)
        
        
        ##Plot the tree (used a circular plot in my case)
        fviz_dend(dart.matx.hcl, k = 5, # Cut in four groups
                  cex = 0.8, # label size,
                  lwd = 1,
                  hang = 1,
                  repel = T,
                  k_colors = "lancet",
                  color_labels_by_k = TRUE, # color labels by groups
                  rect = TRUE, # Add rectangle around groups
                  rect_fill = TRUE,
                  rect_border = "jco", 
                  type = "circular")
        
        #Import the microbial abundance data (otu table) 
        
        ##Read the core files with ASVs converted to mean abundance
        # whole ASVs
        e19.mean.core <- t(read.csv("data/ASV_table_bacteria.tsv",
                                    header = T,
                                    sep = "\t",
                                    row.names = 1,check.names=F)) # 行为样本名
        
        # or core genera						  
        e19.mean.core <- read.csv("data/core bacterial genera.txt",
                                  header = T,
                                  sep = "\t",
                                  row.names = 1,check.names=F) # 行为样本名 
        
        # env factors
        e19.mean.env <- read.csv("data/env.txt",
                                 header = T,
                                 sep = "\t",
                                 row.names = 1,check.names=F) # 行为样本名
        e19.mean.env <- scale(e19.mean.env, center = TRUE, scale = TRUE)						  
        
        #explore it 
        head(e19.mean.core, n=2)
        dim(e19.mean.core)
        str(e19.mean.core)
        
        #change the dataframe to matrix
        e19.mean.core <- as.matrix(e19.mean.core)
        e19.mean.core
        
        
        #Conduct the Mantel Test
        
        ##First generate distance matrices for both data sets and check their dimensions match
        
        ##generate a bray-curtis distance matrix for the microbial data set
        pc1.e19.mean.core.bray <- as.matrix(vegdist(e19.mean.core, 
                                                    method = "bray"))
        
        #explore
        pc1.e19.mean.core.bray
        dim(pc1.e19.mean.core.bray) #check the dimensions
        
        
        
        ##generate a euclidean distance matrix for the dart data set
        pc2.dart.mat.euc <- as.matrix(vegdist(dart.matrix.data, 
                                              method = "euclidean"))
        
        #explore
        pc2.dart.mat.euc
        dim(pc2.dart.mat.euc)
        
        new_rows <- union(row.names(pc1.e19.mean.core.bray), row.names(pc2.dart.mat.euc))
        new_cols <- union(colnames(pc1.e19.mean.core.bray), colnames(pc2.dart.mat.euc))
        
        pc1.e19.mean.core.bray <- pc1.e19.mean.core.bray[new_rows, new_cols]
        pc2.dart.mat.euc <- pc2.dart.mat.euc[new_rows, new_cols]
        
        #Run mantel test
        ##both data sets have been confiremed to have 12x12 dimenstion
        ##so I am good to go here.
        
        #the mantel test formulat has the following format
        #mantel.test <- mantel(fist.distance.matrix, second.distance.matrix, 
        #                                  method="pearson", #method of correlation
        #                                  permutations = 999, #number of permutation
        #                                  strata = NULL, na.rm = FALSE)
        
        ##??matel.test for more
        
        mantel.test.otu.vs.dart <- mantel(pc1.e19.mean.core.bray, pc2.dart.mat.euc, 
                                          method="pearson", 
                                          permutations = 999, 
                                          strata = NULL, na.rm = FALSE)
        #explore the matel test result
        mantel.test.otu.vs.dart
#####################################################################################################
## C-D ## calculate bray_curtis distance and perform PCoA analysis
library(vegan)
library(ape)
library(ggplot2)

## 
otu <- read.delim('data/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- as.data.frame (t(otu))

group=read.delim("data/metadata.tsv",row.names = 1,check.names = F)
otu=otu[match(rownames(group),rownames(otu)),]
dis<-vegdist(otu,method="bray")

bray.pc<-pcoa(dis,correction="cailliez")
summary(bray.pc)
head(bray.pc$values)
bray.pc=bray.pc$vectors[,1:4] 
bray.pc=bray.pc[match(rownames(group),rownames(bray.pc)),]
bray.data <- cbind(bray.pc, group)

colors7=c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2","#96C37D")
ggplot(bray.data, aes(x = Axis.1, y = Axis.2, color=group2)) +
  geom_vline(xintercept=0,alpha = 0.3)+
  geom_hline(yintercept=0,alpha = 0.3)+ 
  geom_point(size =3, alpha = 0.9) + 
  scale_color_manual(values = colors7) +
  theme_bw() +
  labs(x = "PCo 1 (14.98%)", y = "PCo 2 (9.48%)",title="Bacteria") +
  theme(text = element_text(size = 15), legend.key = element_blank(),legend.position="right")+
  guides(color=guide_legend(title="Sites"))+
  theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  theme(axis.text.x=element_text(angle =0,size=10,vjust=0,hjust=0.5,color="black",face = "bold"))+
  theme(axis.text.y=element_text(size=10,color="black",face="bold"),legend.position = "right")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  theme(legend.text=element_text(size=10))

anosim(dis, group$group2,permutations = 999, distance = "bray")
adonis2(dis ~ group$group2) 
       
#############################################################################################################################
## E ## 

# this script refers to Thiergart et al.,(2020): Root microbiota assembly and adaptive differentiation among European Arabidopsis populations

# Collinearity analysis of environmental factors
library(pheatmap)
library(psych)

env=read.delim("data/env.txt",row.names=1,check.name=F,sep="\t")
res <- corr.test(env,env,method = "spearman",adjust = "fdr")

# Setting the Threshold
r_threshold <- 0.7
p_value_threshold <- 0.01

# filter
filtered_data <- res$r
filtered_data[abs(res$r) <= r_threshold | res$p.adj > p_value_threshold] <- NA

pheatmap(filtered_data, fontsize_number=10,fontsize = 10,cluster_rows = F,
cluster_cols = F,main="",cellheight=10,cellwidth = 10)


# the explanation of environmental factors on community composition
env_explain<-function(otu,group,env){
  library(vegan)
  library(tidyverse)
  otu <- as.data.frame (t(otu))#Rows are sample names
  otu=otu[match(rownames(group),rownames(otu)),]#
  mat.dist<-as.matrix(vegdist(otu,method="bray")) #bray-curtis distance
 env=env[match(rownames(group),rownames(env)),]
  colnames(env)=c("pH","NO3","NH4","AK","AP","TK","TP","TN","SOC","CN","Altitude","Longitude","Latitude",'MAT','AT','RH','MAP','AS','RA')
  env=select(env,-CN,-TN,-SOC,-Altitude,-Longitude,-Latitude,-RH,-AS,-AT) # to remove collinear environmental factors
  
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  env2=as.data.frame(lapply(env,normalize))
  rownames(env2)=rownames(env);colnames(env2)=colnames(env)
  
  metadata=cbind(group,env2)
  
  ### loop into the fractions to test separately each fraction ###
  result<-matrix(ncol=4,nrow=12)
  tab<-paste("t", "L", sep="")
  pat<-paste("L", sep = "-")
  aa_pat=grep(pat,rownames(mat.dist))
  aa_inter=intersect(rownames(mat.dist[aa_pat,]),rownames(metadata)) 
  res=adonis(formula= mat.dist[aa_inter,aa_inter] ~ as.vector(metadata[aa_inter,"pH"])+ as.vector(metadata[aa_inter,"NO3"])+ as.vector(metadata[aa_inter,"NH4"])+ as.vector(metadata[aa_inter,"AK"])+ as.vector(metadata[aa_inter,"AP"])+ as.vector(metadata[aa_inter,"TK"])+ as.vector(metadata[aa_inter,"TP"])+as.vector(metadata[aa_inter,"MAT"])+as.vector(metadata[aa_inter,"MAP"])+as.vector(metadata[aa_inter,"RA"]), data=as.data.frame(metadata[aa_inter,]))
  
  for (e in 1:12){
    
    result[e,1]<-res$aov.tab$Df[e]
    result[e,2]<-res$aov.tab$F.Model[e]
    result[e,3]<-res$aov.tab$`Pr(>F)`[e]
    result[e,4]<-res$aov.tab$R2[e]
  }
  
  colnames(result)=c("L_Df","L_F","L_Pr","L_R2")
  row.names(result)=c("pH","NO3","NH4","AK","AP","TK","TP",'MAT','MAP','RA',"Residuals", "Total")
  result
  result=as.data.frame(result)
  
  a=as.data.frame(t(select(result,contains("R2"))))
  a=select(a,-Residuals,-Total,)
  a_long <- as.data.frame(pivot_longer(a, cols = 1:ncol(a),names_to = "element", values_to = "content"))
  a_long$element=factor(a_long$element,levels = c("pH","NO3","NH4","AK","AP","TK","TP",'MAT','MAP','RA'))
  result2=result[match(a_long[,1],rownames(result)),]
  a_long$P=result2$L_Pr
  a_long[which(a_long$P>0.01),"sig"]<-"p>0.01";a_long[which(a_long$P< 0.01),"sig"]<-"p<0.01"
  a_long$sig=factor(a_long$sig,levels=c("p<0.01","p>0.01"))
  
  results=list(a=result,b=a_long)
  return(results)
}

##data of bacteria
otu <- read.delim('data/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
group=read.delim("data/metadata.tsv",row.names = 1,check.names = F)        
env=read.delim("data/env.txt",row.names=1,check.name=F,sep="\t")        
result1=env_explain(otu,group,env)        
df=result1$b
        
ggplot(df,aes(x = reorder(element,-content), y = content*100,fill=sig))+geom_bar(position="dodge",stat="identity")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
  labs(title="The whole community",x="",y="Explained variance (%) ")+
  theme(axis.text.x=element_text(size=10,angle=0,vjust=0,hjust=0.5,color = "black"),
        axis.text.y=element_text(size=10,color = "black"), legend.text=element_text(size=10),legend.position="right",legend.key.size=unit(0.3,"cm"))+
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  guides(fill=guide_legend(title=NULL))+theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  scale_fill_manual(values = c("#6DBB60","gray50"))+coord_flip()

###############################################################################################################################################################################
Fig. 3
###################################### the relative importance of species sorting and dispersal limitation
#Dist_Matrix: A dissimilarity matrix
#Env: Environmental table with samples as rows and environmental variables as columns.
#Geo: Geographic table with samples as rows and geographic/cartesian coordinates as columns.

source("VarPartDist.txt")

otu <- read.delim('data/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- as.data.frame (t(otu))
Dist_Matrix<-as.matrix(vegdist(otu,method="bray")) #bray-curtis


env1=read.delim("data/env.txt",row.names=1,check.name=F,sep="\t")
colnames(env1)=c("pH","NO3","NH4","AK","AP","TK","TP","TN","SOC","CN","Altitude","Longitude","Latitude",'MAT','AT','RH','MAP','AS','AR')
Env=select(env1,-CN,-Longitude,-Latitude)

Geo=select(env1,Longitude,Latitude)
VPD_bac=VarPartDist(Dist_Matrix,Env,Geo,Geo_Co=TRUE,Number_Permutations=999) # SDER=2.01


VPD_fun=VarPartDist(Dist_Matrix,Env,Geo,Geo_Co=TRUE,Number_Permutations=999) # SDER=0.40

df=data.frame(group=c("Bacteria","Fungi"),sder=c(2.01,0.40))

ggplot(df,aes(group,sder,fill=group))+
geom_bar(stat="identity")+
  scale_fill_manual(values=c("#BEB8DC", "#FFBE7A"))+
  theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  labs(x="",y="Sorting/dispersal effect ratio",title="")+
  theme(axis.text.x=element_text(angle =0,size=10,vjust=0,hjust=0.5,color="black",face = "bold"))+
  theme(axis.text.y=element_text(size=10,color="black",face="bold"),legend.position = "right")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  theme(legend.key=element_rect(fill="white"))+
  theme(legend.text=element_text(size=10))


########################## neutral model
source("codes/Neutral model.R")

otu <- read.delim('data/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- as.data.frame (t(otu))

nm_bac=Neutral.fit(otu) # m=0.28, Rsqr =0.92

# fungi: m=0.17, Rsqr =0.95

df=data.frame(group=c("Bacteria","Fungi"),sder=c(0.28,0.17))

ggplot(df,aes(group,sder,fill=group))+
geom_bar(stat="identity")+
  scale_fill_manual(values=c("#BEB8DC", "#FFBE7A"))+
  theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  labs(x="",y="The estimated migration rate (m)",title="")+
  theme(axis.text.x=element_text(angle =0,size=10,vjust=0,hjust=0.5,color="black",face = "bold"))+
  theme(axis.text.y=element_text(size=10,color="black",face="bold"),legend.position = "right")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  theme(legend.key=element_rect(fill="white"))+
  theme(legend.text=element_text(size=10))

########################## A community table with samples as rows and taxa as columns.
library(spaa)
library(plyr)
library(tidyverse)
library(reshape2)
library(ggpubr)

source("codes/habitat niche breadth analysis.R")

otu <- read.delim('data/ASV_table_CSS.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)
niche_bac=as.data.frame(Com.niche(otu));names(niche_bac)="Bacteria"

## fungi is the same to bacteria
niche_fun=as.data.frame(Com.niche(otu));names(niche_fun)="Fungi"

df=cbind(niche_bac,niche_fun); 
df_new=melt(df)

source("codes/summarySE.R")

gsd<-summarySE(df_new, measurevar="value",groupvars=c("variable"))

ggplot(df_new,aes(variable,value,fill=variable))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_boxplot(outlier.size=0.18,colour="black",size=0.18)+
  scale_fill_manual(values=c("#BEB8DC", "#FFBE7A"))+
  theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  labs(x="",y="Habitat niche breadth (Bcom)",title="")+
  theme(axis.text.x=element_text(angle =0,size=10,vjust=0,hjust=0.5,color="black",face = "bold"))+
  theme(axis.text.y=element_text(size=10,color="black",face="bold"),legend.position = "right")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  theme(legend.key=element_rect(fill="white"))+
  theme(legend.text=element_text(size=10))
 
############################import data
load("data/16s.iCAMP.Confidence.detail.rda")

tax<- read.delim("data/taxonomy.tsv",row.names=1,check.names = F,header=T, stringsAsFactors = FALSE)
sp.bin=res$detail$taxabin$sp.bin

tax %>% separate(Taxon,c("kingdom","phylum","class","order","family","genus","species"),sep="[;]")->tax
tax=tax[rownames(sp.bin),];sp.bin=sp.bin[match(rownames(tax),rownames(tax)),]
sp.bin.tax=cbind(sp.bin,tax)
sp.bin.tax$bin.id.new=factor(sp.bin.tax$bin.id.new)
write.table(sp.bin.tax,"data/sp.bin.tax.tsv",sep="\t")

## contribution of each bin to community assembly
b=read.delim("data/16s.BinContributeToProcess_EachGroup.csv",header=T, stringsAsFactors = FALSE,sep=",")
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
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+#
  guides(fill=guide_legend(title=NULL))+theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  scale_fill_manual(values = colors5) #+coord_polar(theta = "x")
p1_1

# top 5-bins
b_long$Bin_ID=gsub('[bin]','',b_long$Bin_ID)
spbin=read.delim("data/sp.bin.tax.tsv",sep="\t",check.name=F)
bin_top5<-subset(spbin,spbin$bin.id.new %in% a[1:5])
bin_top5=subset(bin_top5,!is.na(family))

#write.table(bin_top5,"data/sp.bin.tax.top5.tsv",sep="\t",row.names=F)

colors38<-c('dodgerblue4','darkseagreen','chartreuse4','darkorange','gold','burlywood2','darkolivegreen3','brown3','#984EA3','cyan3','#ed1299', '#09f9f5', '#cc8e12', '#d561dd', '#ddd53e','#4aef7b', '#e86502', '#8249aa', '#ff523f','#246b93','#FFCC00','#4A708B','#FFA54F','#40E0D0','#FF6347','#6DBB60','#00AED1','#8A64A7','#DF7469','#8ECFC9','#FFBE7A','#FA7F6F','#82B0D2','#BEB8DC','#E7DAD2','#96C37D', '#c93f00','grey50')

ggplot(bin_top5,aes(x = class,fill=family))+geom_bar(position="stack")+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
  labs(title="",x="At the class level",y="The number of ASVs")+
  theme(axis.text.x=element_text(size=10,angle=30,vjust=0.5,hjust=0,color = "black"),
        axis.text.y=element_text(size=10,color = "black"), legend.text=element_text(size=10),legend.position="right",legend.key.size=unit(0.3,"cm"))+
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  guides(fill=guide_legend(title=NULL))+theme(strip.text=element_text(size = rel(1),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  scale_fill_manual(values = colors38) #+coord_polar(theta = "x")


## The relative importance of each process to each bin
df=read.delim("data/16s.ProcessImportance_EachBin_EachGroup.csv",header=T,sep=",")
df=df[,-1:-3];rownames(df)=df[,1];df=df[,-1];df=as.data.frame(t(df));df$Bin_ID=rownames(df)
df_1=df[,c(1:5)];df_1=as.data.frame(lapply(df_1,as.numeric));
rownames(df_1)=rownames(df);colnames(df_1)=colnames(df[,1:5])
df_1$Bin_ID=rownames(df_1)
df_2=as.data.frame(df[,c(6,9)])
df_long<-melt(df_1,ID="Bin_ID")
df_long<-merge(df_long,df_2,by="Bin_ID",all=TRUE)

df_long$variable=factor(df_long$variable,levels=c("HeS","HoS","DL","HD","DR"))
df_long$DominantProcess=factor(df_long$DominantProcess,levels=c("HoS","DL","DR"))

###############################################################################################################################################################################
Fig. 4
###################################### core genera
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


###############################################################################################################################################################################
Fig. 5
######################################
################################################################################################################################
## A-B #
conda activate qiime2-2021.2

qiime tools import \
  --input-path isolates.fasta \
  --output-path seq.qza \
  --type 'FeatureData[Sequence]'

ime qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences seq.qza \
  --o-alignment aligned-seq.qza \
  --o-masked-alignment masked-aligned-seq.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza \
  --p-n-threads 70
  
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path 16S_isolates
  
## Phylogenetic tree was viewed in iTOL.
  
################################################################################################################################
## C ##
library(EasyAovWlxPlot)
library(ggplot2)

data <- read.delim("data/metadata2.txt", sep = "\t",  check.names = FALSE)

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
library(reshape2)

df=read.table("data/final_isolates_taxonomy.tsv",sep="\t",row.names=1,header=T)

df_new=melt(df)

colors3=c("#D95F02", "#8A64A7", "#00AED1")
df_new$Genus=factor(df_new$Genus,levels=c("Bacillus","Pantoea","Methylobacterium-Methylorubrum"))
df_new$Feature_ID=factor(df_new$Feature_ID,levels=c(df$Feature_ID))

ggplot(df_new,aes(Feature_ID,value,fill=Genus))+
    geom_bar(stat="identity")+
	facet_wrap(.~variable,1)+
	scale_fill_manual(values = colors3)+
    theme(strip.text=element_text(size = rel(0.7),face="bold"),strip.background=element_rect(fill="white",colour="transparent"))+
    theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent',color='black'))+
    labs(title="Bacterial isolates",x="",y="Percentage inhibition (100%)",fill="")+
    theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black',))+
    theme(axis.text.x=element_text(size=8,angle=0,vjust=0,hjust=0.5,color = "black",face = "bold"),
          axis.text.y=element_text(size=8,color = "black",face = "bold"),legend.position = "right")+
    #theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
	coord_flip()
















