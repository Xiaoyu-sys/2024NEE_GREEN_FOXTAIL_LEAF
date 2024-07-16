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
