#args = commandArgs(T)

library(vegan)
library(tidyverse)
library(philentropy)
#library(ggplot2)

meta <- read.table("./nasal_meta.tsv",sep = "\t", header = T) ##inout the metadata

##Explore the data distribution and Filter outliers
ggplot(meta) + geom_density(aes(x=Age)) +
  theme(panel.background = element_rect(fill='white',colour = 'black'))


#get_Standance <- function(meta,feature,n){
#  feature <- select(meta,feature) %>%  unlist %>% as.numeric()
#  m <- round(mean(feature, na.rm = T), 2)
#  sd <- round(sd(feature, na.rm = T), 2)
#  m1 <- m + n*sd
#  m2 <- m - n*sd
#  feature[feature >= m1 | feature <= m2] <- NA
#  return(feature)
#}

#meta <- data.frame(Simple_id = meta$Sample_id, 
#                   Age = get_Standance(meta = meta, feature = "Age", n=3),
#                   BMI = get_Standance(meta = meta, feature = "BMI", n=3),
#                   Sex = meta$Sex)


##mulit Permanova for Age, BMI, Sex
get_muliPerm <- function(df, meta, med){
  if(med=="JSD"){
    df<-as.matrix(df)
    dist <- JSD(df)}else{
      dist <- vegdist(df,method = med)}
  res_p <- adonis2(dist ~ Sex + BMI + Age, meta, permutations = 4999,
                   binary=F, na.action=na.omit, by="margin") %>% as.data.frame()
  return(res_p)
}

getRes <- function(meta, med){
  res <- list()
  n <- c("bacteria","fungi","total")
  for(i in 1:3){
    
    ##input the taxa profile
    df <- read.table(paste0(n[i],"_occ0.tsv"), sep = "\t", header = T)
    
    ##filter Occurrence rate > 10%
    df <- df[,which(colMeans(df > 0) > 0.1)]
    
    ##comfirm the taxa file sample names same as the metadata sample name
    if(identical(rownames(df),meta$Sample_id)){
      result <- get_muliPerm(df=df, meta = meta, med = med)
      res[[i]] <- result
    }
  }
  return(res)
}

#P_values FDR correction for permanova result
getFinal <- function(df1,df2,df3,df4,n){
  res <- list()
  n <- c("bacteria","fungi","total")
  for (i in 1:3) {
    res1 <- rbind(df1[[i]],df2[[i]],df3[[i]],df4[[i]])
    res1$FDR <- p.adjust(res1$`Pr(>F)`, method = "BH")
    res[[i]] <- res1
  }
  return(res)
}

meta <- data.frame(Sample_id = meta$Sample_id,
                   Age = meta$Age,
                   BMI = meta$BMI,
                   Sex = meta$Sex)


set.seed(1)
JSD <- getRes(med = 'JSD')
Bray <- getRes(med = "bray")
Jacca <- getRes(med ="jaccard" )
Eucli <- getRes(med ="euclidean")

result <- getFinal(df1 = Bray, df2 = JSD, df3 = Jacca, df4 = Eucli)

write.table(ressult,"mulit_Permenova_result.tsv")
