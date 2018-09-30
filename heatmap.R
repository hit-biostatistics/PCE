library(ggplot2)
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#library(limma)
library(Linnorm)
library(Rtsne)
library(DESeq2)
library(samr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clues)
library(pheatmap)

#得到聚类数 需要先跑PCA，princomp的结果传给x
get_clusternum<-function(x,percent=0.959){
  sumsqr<-sum(x^2)
  psum=0
  i=1
  while(psum<=percent){
    psum=psum+(x[[i]]^2)/(sumsqr)
    i=i+1
  }
  return(i-1)
}

#归一化方法
normalize<-function(x,y){
  if (y==1)
  {
    ES1<-log2(as.matrix(x)+1)
    return (ES1)
  }
  if (y==2)
  {
    ES1<-scale(t(x))
    
    return (t(ES1))
  }
  if (y==3)
  {
    ES1<-Linnorm(x)
    return (ES1)
  }
}

#降维方法
dim_reduce<-function(x,y){
  if (y==1)
  {
    TS2<-Rtsne(t(x),dim=2,perplexity = 15)
    return (as.matrix(TS2$Y))
  }
  if (y==2)
  {
    res <- prcomp(t(x),scale. = T)
    res<-predict(res)
    return(res[,1:2])
  }
  if (y==3)
  {
    kk2 <- as.data.frame(x)
    kk2 <- dimRedData(data=kk2)
    res <-  LLE()
    kk2 <- res@fun(kk2,res@stdpars)
    return(as.matrix(kk2@data@data))
  }
  if (y==4)
  {
    kk2 <- as.data.frame(x)
    kk2 <- dimRedData(data=kk2)
    res <-  LaplacianEigenmaps()
    kk2 <- res@fun(kk2,res@stdpars)
    return(as.matrix(kk2@data@data))
  }
  
}

#聚类方法
cluster<-function(x,y,z){
  if (y==1)
  {
    ires<-kmeans(x,z)
    return(ires$cluster)
  }
  if (y==2)
  {
    idis<-dist(x,method = "euclidean")
    clust2_1norm<-hclust(idis,method="average")
    ires<-cutree(clust2_1norm,k=z)
    return (ires)
  }
  
  if (y==3)
  {
    resc<-ecm(x=x,c=z)
    return(resc$y.pl)
  }
  if (y==4)
  {
    y0 <- sample(1:z,replace=TRUE,nrow(x))
    res<-EkNNclus(x=x,K=z,y0=y0)
    return (res$y.pl)
  }
  if (y==5)
  {
    res <- adpclust(x,centroids="auto")
    tans <- array(dim=1)
    limit=nrow(x)
    for (i in 1:limit)
    {
      tans[i]<-res$clusters[i]
    }
    return (tans)
  }
  if (y==6)
  {
    res <- dbscan(x=x,eps=10)
    tans <- array(dim=1)
    limit=nrow(x)
    for (i in 1:limit)
    {
      tans[i] <- res$cluster[i]
    }
    return (tans)
  }
  if (y==7)
  {
    ap_clust <- apcluster(negDistMat(r=2), x)
    limit=length(ap_clust@clusters)
    temprow=nrow(x)
    res <- array(dim=1)
    for (i in 1:limit)
    {
      limitj=length(ap_clust[[i]])
      for (j in 1:limitj)
      {
        temp <- as.numeric(ap_clust[[i]][j])
        res[temp]=i
      }
    }
    return(res)
  }
}

#t-test
T_matrix<- function(df,df_label){     #label 从1开始的整数 
  out<-data.frame(row.names = rownames(df))
  for( i in min(df_label):max(df_label) ){
    for( j in 1:dim(df)[1]){
      a<-as.vector(t(df[j,df_label==i])) #cluster i 的第j 个基因表达量
      b<-as.vector(t(df[j,df_label!=i]))
      if(sd(a)==0 && sd(b)==0){
      }else{
        #t检验
        out[j,i]<-t.test(a,b, alternative="greater", paired = FALSE, var.equal = FALSE)$p.value
      }
    }
    print(i)
  }
  colnames(out)<-paste("cluster",1:max(df_label))
  return(out)
} #t 检验

#求非参数检验得到的差异表达基因P值
W_matrix<- function(df,df_label){     #label 从1开始的整数 
  out<-data.frame(row.names = rownames(df))
  for( i in min(df_label):max(df_label) ){
    for( j in 1:dim(df)[1]){
      a<-as.vector(t(df[j,df_label==i])) #cluster i 的第j 个基因表达量
      b<-as.vector(t(df[j,df_label!=i]))
      if(sd(a)==0 && sd(b)==0){
      }else{
        #t检验
        out[j,i]<-wilcox.test(a,b, alternative="greater")$p.value
      }
    }
  }
  colnames(out)<-paste("cluster",1:max(df_label))
  return(out)
} # W检验

DE_matrix<-function(df, label,id){
  
  y<-c()
  for(i in 1:length(label)){
    if(label[i]==id) y[i]="treat"
    else y[i] = "control"
  }
  
  group_list <- factor(y,levels = c("control","treat"))
  gene_names = rownames(df)
  
  dds <- DESeqDataSetFromMatrix(countData = df,
                                colData = DataFrame(group_list,
                                                    row.names = colnames(df)),
                                design = ~ group_list)
  #dds <- estimateSizeFactors(dds, type="poscounts")
  dds2<-DESeq(dds,parallel=TRUE)
  res <- results(dds2)
  DE<-data.frame(gene_names,res$padj,res$log2FoldChange)
  #DE <-DE[DE$res.log2FoldChange>0,]
  #DE <- DE[order(DE$res.padj),]
  #DE.diff_gene <- as.vector(DE$gene_names[1:topn])
  return(DE)
}


#提取以上p值矩阵对应的基因名，k=提取 cluster k的差异表达基因
d_expr<-function(out_matrix,k=1,p=0.01,gene_names){   # k为输出第K类的差异表达，p是输出小于P的P值
  pvalue<-out_matrix[,k]
  times<-0
  gene_list<-c()
  for( i in 1:dim(out_matrix)[1]){
    if(!is.na(pvalue[i])&pvalue[i]<p){
      times<-times+1
      gene_list[times]<-gene_names[i]
    }
  }
  return(gene_list)
}

#yan
if(F){
yan<-read.csv("D:/single cell data/6golden/yan2013/yan.csv")
gene_names<-as.vector(yan$Gene_ID)
yan<-read.csv("D:/single cell data/6golden/yan2013/nsmb.2660-S2.csv")
yan_label<-rep(1:7,c(3,3,6,12,20,16,30))

exprs_matrix<-yan
rm(yan)

rownames(exprs_matrix)<-gene_names
author_label<-yan_label
rm(yan_label)
}

#goolam
if(F){
goolam<-read.table("D:/single cell data/6golden/gloom 2016/Goolam_et_al_2015_count_table.tsv")
goolam.label<-rep(c(1,2,3,4,5,3,1,2),c(6,40,16,6,6,16,10,24))
exprs_matrix<-goolam
author_label<-goolam.label
rm(goolam)
rm(goolam.label)
}

#biase
if(T){
  biase<-read.table("D:/single cell data/6golden/biase2014/GSE57249_fpkm.txt",
                    header = T)
  rownames(biase)<-as.vector(biase$ID)
  biase<-biase[-1]
  exprs_matrix<-biase[-(50:56)]
  rm(biase)
  author_label<-c(rep(1,9),rep(2,20),rep(3,20))
}

#Zeisel
if(F){
  Zeisel<-read.table("D:/single cell data/6silver/#Zeisel 2015/Zeisel.txt",
                     header = T)
  gene_names<-as.vector(Zeisel$cell_id)
  exprs_matrix<-Zeisel[,-1]
  #exprs_matrix<-Zeisel[rowSums(Zeisel[,-1]>2)>5,-1]
  rm(Zeisel)
  Zeisel_marker<-c('Aif1','Aldoc','Acta2','Cldn5','Thy1','Spink8','Mbp','Gad1','Tbr1')
  
  
  a<-dim_reduce(t(exprs_matrix),1)
  Z_plotdf<-data.frame(x1=a[,1],
                          x2=after_reduce[,2],
                          row.names = paste("Cell",1:dim(exprs_matrix)[2]))
  
  reds<-quantile(as.matrix(exprs_matrix),0.95)
  blues<-quantile(as.matrix(exprs_matrix),0.05)
  Z_color<-list()
  for (i in 1:length(Zeisel_marker)){
    expr_vec<-as.vector(exprs_matrix[gene_names==Zeisel_marker[i],])
    q<-c()
    for (j in expr_vec){
      if(j>=reds) q<-c(q,3)
      else if(j<=blues) q<-c(q,0)
      else q<-c(q,1)
    }
    Z_color[[i]]<-q
  }
  p<-ggplot(PEMC_plotdf,aes(x1,x2)) + 
    geom_point(aes(colour=as.character(Z_color[[4]])),size=1)+
    labs(y="tSNE2",x = "tSNE1")
  cbbPalette <- c("#0000FF","#EAADEA","#FF0000")
  p+scale_colour_manual(values=cbbPalette)+ theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
}

#Mariella
if(F){
  Mariella<-read.table("D:/single cell data/Mariella G/GSE102130_K27Mproject.RSEM.vh20170621.txt",
                       header = T)
  gene_names<-Mariella$Gene
  exprs_matrix<-Mariella[,-1]
  rm(Mariella)
  exprs_matrix<-as.matrix(exprs_matrix)
}

#Romanov
if(F){
  exprs_matrix<-read.csv("D:/single cell data/Romanov RA/hypoth_moldata_classification08-Mar-2017.csv")
  level1<-as.vector(t(exprs_matrix[1,-1]))
  level2<-as.vector(t(exprs_matrix[2,-1]))
  level2_cluster<-as.vector(t(exprs_matrix[3,-1]))
  gene_names<-as.vector(exprs_matrix[12:24356,1])
  exprs_matrix<-exprs_matrix[12:24356,-1]
  exprs_matrix<-as.matrix(exprs_matrix)
  author_label<-c()
  for (i in level1){
    if (i == "astrocytes") author_label<-c(author_label,1)
    else if (i == "endothelial") author_label<-c(author_label,2)
    else if (i == "ependymal") author_label<-c(author_label,3)
    else if (i == "microglia")author_label<-c(author_label,4)
    else if ( i == "neurons") author_label<-c(author_label,5)
    else if (i == "oligos") author_label<-c(author_label,6)
    else author_label<-c(author_label,7)
  }
}
#n=3 Linnorm,   r=1 Tsne,  c=2  hclust
n=3;r=1;c=2;k=3
#k=get_clusternum(pca1,0.96)

after_reduce<-dim_reduce(normalize(exprs_matrix,n),r)
PEMC_label<-cluster(after_reduce,c,k)
PEMC_hc<-hclust(dist(after_reduce, method = "euclidean"),method = "average")
adjustedRand(author_label,PEMC_label)


PEMC_plotdf<-data.frame(x1=after_reduce[,1],
            x2=after_reduce[,2],
            row.names = paste("Cell",1:dim(exprs_matrix)[2]))
clusters<-factor(PEMC_label)

p<-ggplot(PEMC_plotdf,aes(x1,x2)) + 
  geom_point(aes(colour=clusters),size=2)+
  labs(y="tSNE2",x = "tSNE1")
cbbPalette <- c("#000000","#FF0000","#0000FF","#33FF33","#999999","#FFFF00",
                "#00FFFF", "#FF00FF", "#993300", "#006699", "#FF9966", "#006600",
                "#AAAAAA")
p+scale_colour_manual(values=cbbPalette)+ theme_bw()+
  #scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+
  theme(title=element_text(size=15),axis.title=element_text(size=15),
        legend.text = element_text(size=15),
    panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

after_linnorm<-Linnorm(exprs_matrix)
T_diff_gene_pvalue<-T_matrix(after_linnorm, PEMC_label)
W_diff_gene_pvalue<-W_matrix(after_linnorm, PEMC_label)

get_DESeq_gene_names<-function(topn = 10, cut_positaion = 1){
  diff_gene_list <- list()
  for (i in 1:length(PEMC_label)){
    diff_gene_list[[i]] <-DE_matrix(PEMC_label,exprs_matrix,i)
  }
  
  pgene<-c()
  for (i in 1:k){
    df <- exprs_matrix
    DE<-na.omit(diff_gene_list[[i]])
    DE <-DE[DE$res.log2FoldChange>0,]
    DE <- DE[DE$res.padj<1e-2,]
    df<-df[DE$gene_names,]
    cur<-df[,PEMC_label==i]
    other<- df[, PEMC_label!=i]
    other_median<-apply(other, 1, median)
    #  other_max<-apply(other,1,max)
    #  cur_sd<-apply(cur,1,sd)
    #  cur_mean<-apply(cur,1,mean)
    #  cur_median<-apply(cur,1,median)
    #  cur_min<-apply(cur,1,min)
    cur<-cur[other_median<cut_positaion,]
    pgene<-c(pgene, rownames(cur[order(rowSums(cur),decreasing = TRUE),])[1:topn])
  }
  return(pgene)
}


get_diff_exprs_gene_names<-function(diff_pvalue_matrix,topn=10,p_value = 1e-2){
  pgene<-c()
  for (i in 1:k){
    df<-after_linnorm
    igene<-d_expr(diff_pvalue_matrix,i,p_value,gene_names = rownames(exprs_matrix))
    df<-df[igene,]
    cur<-df[,PEMC_label==i]
    other<-df[,PEMC_label!=i]
    other_median<-apply(other, 1, median)
    other_mean<-apply(other,1,mean)
    cur_median<-apply(cur,1,median)
    cur<-cur[other_median<0.5 & other_mean<0.5 &cur_median>3,]
    pgene<-c(pgene, rownames(cur[order(rowSums(cur),decreasing = T),])[1:topn])
  }
  return(pgene)
}

numsDiffGene<-10
T_pgene<-get_diff_exprs_gene_names(T_diff_gene_pvalue,topn =numsDiffGene,p_value = 1e-5)
W_pgene<-get_diff_exprs_gene_names(W_diff_gene_pvalue,topn = numsDiffGene, p_value = 1e-5)
#DE_pgene<-get_DESeq_gene_names(DE_diff_gene_pvalue)

yanclusterOrder<-c(8,7,3,2,1,6,5,4)
goolamClusterOrder<-c(3,1,2)
get_reOrderGene<-function(gene,clusterOrder){
  reOrderGene<-c()
  for (i in clusterOrder){
    reOrderGene<-c(reOrderGene,gene[((i-1)*numsDiffGene+1):(i*numsDiffGene)])
  }
  return(reOrderGene)
}
reOrderGene<-get_reOrderGene(T_gene,goolamClusterOrder)
reOrderGene

Zeisel_marker<-c('Aif1','Aldoc','Acta2','Cldn5','Thy1','Spink8','Mbp','Gad1','Tbr1')

ht_df<-data.frame(row.names = Zeisel_marker)
for (i in 1:length(Zeisel_marker)){
  expr_vec<-after_linnorm[gene_names==Zeisel_marker[i],]
  for(j in 1:length(expr_vec)){
    ht_df[i,j]<-expr_vec[j]
  }
}

annotation_col<-data.frame(Cluster = factor(as.character(PEMC_label)))
rownames(annotation_col)<-colnames(ht_df)
ann_colors = list(Cluster = c("1" = "#999999", "2"="#FF00FF","3"="#FF0000",
                              "4"="#993300","5" = "#000000","6"="#FFFF00",
                              "7" = "#33FF33","8"= "#00FFFF","9"="#0000FF",
                              "10"="#006699","11" =  "#FF9966","12"="#006600")[1:k])
pheatmap(ht_df,cluster_rows = FALSE,annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         cluster_cols =  PEMC_hc, 
         gaps_row = 1:7, 
         cutree_cols = k,
         show_colnames =  F)





ht_df<-after_linnorm[reOrderGene,]
annotation_col<-data.frame(Cluster = factor(as.character(PEMC_label)))
rownames(annotation_col)<-colnames(ht_df)
ann_colors = list(Cluster = c("1" = "#999999", "2"="#FF00FF","3"="#FF0000",
                              "4"="#993300","5" = "#000000","6"="#FFFF00",
                              "7" = "#33FF33","8"= "#00FFFF","9"="#0000FF",
                              "10"="#006699","11" =  "#FF9966","12"="#006600")[1:k])
pheatmap(ht_df,cluster_rows = FALSE,annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         cluster_cols =  PEMC_hc,
         gaps_row = seq(1,k-1)*numsDiffGene, 
         cutree_cols = k,
         show_colnames =  F)



