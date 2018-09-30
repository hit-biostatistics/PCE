#library(stats4)
#library(BiocGenerics)
#library(parallel) 
#library(cluster)
#library(clues)
#library(apcluster)
#library(evclust)
#library(ADPclust)
#library(kernlab)
#library(dbscan)
#library(fclust)
#library(lle)
#library(dimRed)
library(ggplot2)
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#library(limma)
library(Linnorm)
library(Rtsne)
library(DESeq2)
library(samr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clues)

get_clusternum<-function(x,percent=0.959){
  sumsqr<-sum(x^2)
  psum=0
  i=1
  while(psum<=percent){
    psum=psum+(x[[i]]^2)/(sumsqr)
    i=i+1
  }
  return(i=i-1)
}
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
  }
  colnames(out)<-paste("cluster",1:max(df_label))
  return(out)
}
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
}
d_expr<-function(out_matrix,k=1,p=0.01,gene_names){   # k为输出第K类的差异表达，p是输出小于P的P值
  pvalue<-out_matrix[,k]
  times<-0
  gene_list<-c()
  for( i in 1:dim(out_matrix)[1]){
    if(!is.na(pvalue[i])&&pvalue[i]<p){
      times<-times+1
      gene_list[times]<-gene_names[i]
    }
  }
  return(gene_list)
}
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
h_expr<-function(df,df_label,cluster_id,top = 200){
  clu_i<-df[,which(df_label==cluster_id)]
  clu_i$mean<-apply(df,1,mean)
  out<-as.vector(t(rownames(clu_i[order(clu_i$mean),])))[1:top]
  return(out)
}

#Zeisel  silver 聚类结果13类，对不上 人类脑细胞
if(F){
  Zeisel<-read.table("D:/Asingle cell data/6silver/Zeisel/Zeisel.txt",
                     header = T)
  Zeisel<-Zeisel[!duplicated(Zeisel[,1]),]
  gene<-as.vector(Zeisel[rowSums(Zeisel[,-1]>2)>5,1])
  #get_clusternum(Zeisel)
  Z_data<-Zeisel[rowSums(Zeisel[,-1]>2)>5,-1]
  rownames(Z_data)<-gene
}
#Usoskin   silver 
if(F){
  Uso<-read.csv("D:/Asingle cell data/6silver/Usoskin 2015/Usoskin.csv")
  gene<-as.vector(Uso[,1])[c(-1:-10)]
  type<-c('NF1','NF2','NF3','NF4','NF5','NP1','NP2','NP3','PEP1','PEP2','TH')
  RPM<-Uso[,t(Uso[8,])%in%type]
  author_label<-as.vector(t(RPM[8,]))
  for(i in 1:622){
    for(j in 1:11){
      if(author_label[i]==type[j]){
        author_label[i]=j
        break
      }
    }
  }
  RPM<-RPM[c(-1:-10),]
  RPM<-as.numeric(as.matrix(RPM))
  RPM<-matrix(data = RPM,ncol = 622 )
  data<-RPM[rowSums(RPM>2)>5,]
}
#Klein
if(F){
  Klein<-read.table("")
}
#Patel  按原作者聚类数ARI很高 数据已经做过预处理 基因只剩5000多 PCA聚类数200-300
if(F){
  Patel<-read.table("D:/Asingle cell data/6silver/Patel 2014/Patel.txt")
  Patel<-Patel[,1:430]
  author_label<-c(rep(1,118),rep(2,94),rep(3,75),rep(4,73),rep(5,70))
  pca1<-princomp(2^Patel-1)$sdev
  sumsqr<-sum(pca1^2)
  psum=0
  i=1
  while(psum<=0.959){
    psum=psum+(pca1[[i]]^2)/(sumsqr)
    i=i+1
  }
  i
  gene<-Patel
}
#Ting  基因只有9000多 基因名对不上
if(F){
  Ting<-read.table("D:/Asingle cell data/6silver/Ting2014/TING.txt",sep = '\t')
  gene<-as.vector(Ting$V4[-1])
  Ting<-Ting[2:9901,7:193]
}
#Treutlein
if(F){
  Tre<-read.table("D:/Asingle cell data/6silver/Treutlein 2014/nature13173-s4.txt")
  gene<-as.vector(t(Tre[1,5:23275]))
  FPKM<-as.numeric(as.matrix(Tre[2:81,5:23275]))
  data<-t(matrix(FPKM,nrow = 80))
  gene<-gene[rowSums(data>0)>2]
  data<-data[rowSums(data>0)>2,]
  author_label<-as.vector(Tre$V4[2:81])
  type<-c('AT1','AT2',"BP","ciliated","Clara")
  for(i in 1:length(author_label)){
    for(j in 1:length(type)){
      if(author_label[i]==type[j]) author_label[i]<-j
    } 
  }
  author_label<-as.integer(author_label)
  pca1<-princomp(FPKM)$sdev
  get_clusternum(pca1,percent = 0.96)
  
}
#Darmanis 
if(F){
  path <- "D:/Asingle cell data/Darmanis/GSE67835_RAW" ##文件目录
  fileNames <- dir(path)  ##获取该路径下的文件名
  filePath <- sapply(fileNames, function(x){ 
    paste(path,x,sep='/')})   ##生成读取文件路径
  data <- lapply(filePath, function(x){
    read.csv(x, header=T)})  ##读取数据，结果为list
}
#T_S  
if(F){
  T_S<-read.table("D:/Asingle cell data/Tirosh S/T_S.txt")
  gene<-as.vector(T_S[,1])[c(-1:-4)]
  type<-c(1:6)
  TPM<-T_S[,t(T_S[4,])%in%type]
  author_label<-as.vector(as.numeric(t(TPM[4,])))
  TPM<-TPM[c(-1:-4),]
  author_label<-scan("D:/Asingle cell data/Tirosh S/label.txt",what=character())
  TPM<-read.table("D:/Asingle cell data/Tirosh S/Ts.txt")
  gene<-scan("D:/Asingle cell data/Tirosh S/gene.txt",what = character(),sep = '\n')
  TPM1<-2^TPM-1
  get_clusternum(TPM)
  get_clusternum(TPM1)
}
#T_N 
if(F){
  T_N<-read.table("D:/Asingle cell data/Tirosh N/T_N.txt")
  TPM<-2^T_N-1
  TPM[is.nan(TPM)]<-0
  pca1<-scan(as.vector(pca1),what = numeric(),file = "D:/Asingle cell data/Tirosh N/TNpca.txt",sep = '\n')
  sumsqr<-sum(pca1^2)
  psum=0
  i=1
  while(psum<=0.959){
    psum=psum+(pca1[i]^2)/(sumsqr)
    i=i+1
  }
  
}
#Marques  UMI counts 23556*5069 分2类
if(F){
  Mar<-read.table("D:/Asingle cell data/Marques/Marques.tab")
  pca1<-scan("D:/Asingle cell data/Marques/Mar_pc5a.txt",what=numeric())
  sumsqr<-sum(pca1^2)
  psum<-0
  i<-1
  while(psum<=0.965){
    psum=psum+(pca1[i]^2)/(sumsqr)
    i=i+1
  }
  i<-i-1
}


#data<-TPM[rowSums(TPM>2)>5,]
#clustering and plot
data<-Patel
n=3;r=1;c=2;k=5 #n=3 Linnorm,   r=1 Tsne,  c=2  hclust
after_norm<-normalize(data,n)
after_reduce<-dim_reduce(data,r)
mylabel<-cluster(after_reduce,c,k)
adjustedRand(author_label,mylabel)
df<-data.frame(x1=after_reduce[,1],x2=after_reduce[,2],row.names = paste("Cell",1:dim(data)[2]))
ggplot(df,aes(x1,x2)) + 
  geom_point(colour=mylabel,size=2)


#t-test or wilcoxon
label<-mylabel;p_value=1e-2;topn=1000;n=1
diffw<-W_matrix(data,label)
difft<-T_matrix(data,label)

w<-d_expr(diffw,k=n,p=p_value,gene_names = gene)
t<-d_expr(difft,k=n,p=p_value,gene_names = gene)
#high1<-h_expr(data,label,cluster_id =n,top = topn)
#out<-intersect(high1,cluster1)
#DESeq2
if(T){
  y<-c()
  for(i in 1:length(label)){
    if(label[i]==n) y[i]=1
    else y[i] = 2
  }
  group_list=factor(y)
  test<-as.integer(round(3^data-1))
  test<-matrix(test,ncol = 80)
  test[is.na(test)]<-0
  dds <- DESeqDataSetFromMatrix(countData = test,
                                colData = DataFrame(group_list),
                                design = ~ group_list)
  dds <- estimateSizeFactors(dds, type="poscounts")
  dds2<-DESeq(dds)
  res <- results(dds2)
  DE<-data.frame(gene,res$pvalue)
  DE<-as.vector(t(DE[DE$res.pvalue<p_value,]$gene))
}
# SAMseq y must take values 1,2
if(T){
  y<-c()
  for(i in 1:length(label)){
    if(label[i]==n) y[i]=1
    else y[i] = 2
  }
  samfit <- SAMseq(test, y, resp.type = "Two class unpaired", nperms = 10, random.seed = NULL, nresamp = 20, fdr.output = 0.2) 
  SA<-gene[as.numeric(samfit$siggenes.table$genes.up[,2])]
}
out1<-intersect(t,w)
out2<-intersect(DE,SA)
out<-intersect(out1,out2)

tran<-select(org.Mm.eg.db, keys = out1 ,columns ="ENSEMBL",keytype ="SYMBOL")
out<-as.vector(na.omit(tran)[,2])
write(out,file = "D:/c++ code/cell type/test.in",sep = "\n")

