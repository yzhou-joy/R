#01.vector
vector1<-c(1,4,5)
vector1[2]=5  #change the second element
vector1<-c(vector1,7)  #Add one element
for (i in vector1){
print (i)
}  # for loop

#01.matrix
matrix1<-matrix(seq(1:8),nrow=2,ncol=4,byrow=TRUE)
matrix1[1,2]=9
matrix2<-rbind(matrix1,c(1,1,1,1))
for (m in 1:nrow(matrix2)){
  for (n in 1:ncol(matrix2)){
    cat(matrix2[m,n])
  }
  cat("\n")
}   #for loop to print matrix

#01.dataframe
names<-c("xx","yy","zz")
a<-c(24,19,30)
d<-data.frame(name=names,age=a)
d[1,2]=26
d[["weight"]]<-c(50,60,55)
for (m in 1:nrow(d)){
  print(d[m, ])
}   #for loop to print dataframe

#01.list
l=list(id='aa',height=1.7)
l[["height"]]=1.8
l[["dob"]]=c(1960,12,1)
for (m in 1:length(l)){
  print(l[m])
}   #for loop to print matrix

#02
setwd("E:\\研究僧\\研究生课程\\用R语言分析新一代测序数据\\第一次上机") #这句可以用中文
genes=read.table("hg19_gene_table.txt", comment="", header=TRUE)      #该句路径中不能包括中文

#03
nrow(genes)   #导入时第一行作为标题，统计行数时不算标题行

#04
table(genes$chrom)  #统计chrom列每个染色体上基因数量

#05
table(genes$strand)  #根据strand列统计正负链基因数量

#06
geneLen=genes$txEnd-genes$txStart   #调用txEnd列和txStart列求每个基因的长度

#07
hist(geneLen,50, main="Gene lengths", xlab="base pairs")  #基因长度的画图

#08
hist(log(geneLen),50, main="log(Gene lengths)", xlab="base pairs") #基因长度取log的直方图

#09
boxplot(geneLen~genes$strand) ## this doesn’t show very well

#10
boxplot(log(geneLen)~genes$strand) ## no difference

#11  Hypothesis testing
t.test(geneLen~genes$strand)

#12
boxplot(geneLen ~ genes$chrom,par(las="2"),cex.axis=1)

#13
idx=c(grep("random", genes$chrom), grep("hap", genes$chrom))  #带random、hap的染色体的基因信息行
par(las=3)
boxplot(log(geneLen[-idx]) ~ as.character(genes$chrom[-idx]),cex.axis=0.5) #cex.axis为坐标字的大小；as.character转化为字符
boxplot(log(geneLen[-idx]) ~ (genes$chrom[-idx]),cex.axis=0.5)   #这样是凌乱的

#14 最长的gene
idx=which.max(geneLen)
genes[idx, ]

#15 外显子个数
nExon=genes$exonCount
hist(nExon,40, main="Number of exons", xlab="Exon count")

#16  拥有最多外显子的基因
idx=which.max(nExon)
genes[idx,1:9]   #显示idx基因行的1-9列，10、11两列为外显子起始、终止位置，拥有最多外显子的这两行显示呈现多页

#17  基因长度和外显子数之间的关系
plot(geneLen, nExon, pch=".", xlab="Gene length", ylab="Number of exons")
cor(nExon, geneLen) ## weak correlation, need to run a test (lm or glm).
glm(nExon~geneLen)  #family等参数还不清楚含义

#18 染色体长度与染色体号的对应定义
chrlen=c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954)
names(chrlen)=paste("chr",c(1:22,"X","Y"),sep="")

#19  每个染色体上每个兆碱基对的基因数
chr_gene<-data.frame(chrlen, genes_num=rep(0,24))   #创建数据框，基因数目为0
for (i in rownames(chr_gene)){
  chr_gene[rownames(chr_gene) %in% i,2]<-nrow(genes[genes$chrom %in% i,])    #每个染色体上的基因数目 %in% 相当于 ==
}
gene_density_per_M<-round(chr_gene$genes_num/chr_gene$chrlen*10^6,2)  #round四舍五入
genes_density<-cbind(chr_gene, gene_density_per_M)
genes_density
genes_density[order(-genes_density$gene_density_per_M),]  #按照gene_density_per_M排序 -为从高到低


