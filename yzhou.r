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
setwd("E:\\�о�ɮ\\�о����γ�\\��R���Է�����һ����������\\��һ���ϻ�") #������������
genes=read.table("hg19_gene_table.txt", comment="", header=TRUE)      #�þ�·���в��ܰ�������

#03
nrow(genes)   #����ʱ��һ����Ϊ���⣬ͳ������ʱ���������

#04
table(genes$chrom)  #ͳ��chrom��ÿ��Ⱦɫ���ϻ�������

#05
table(genes$strand)  #����strand��ͳ����������������

#06
geneLen=genes$txEnd-genes$txStart   #����txEnd�к�txStart����ÿ������ĳ���

#07
hist(geneLen,50, main="Gene lengths", xlab="base pairs")  #���򳤶ȵĻ�ͼ

#08
hist(log(geneLen),50, main="log(Gene lengths)", xlab="base pairs") #���򳤶�ȡlog��ֱ��ͼ

#09
boxplot(geneLen~genes$strand) ## this doesn��t show very well

#10
boxplot(log(geneLen)~genes$strand) ## no difference

#11  Hypothesis testing
t.test(geneLen~genes$strand)

#12
boxplot(geneLen ~ genes$chrom,par(las="2"),cex.axis=1)

#13
idx=c(grep("random", genes$chrom), grep("hap", genes$chrom))  #��random��hap��Ⱦɫ��Ļ�����Ϣ��
par(las=3)
boxplot(log(geneLen[-idx]) ~ as.character(genes$chrom[-idx]),cex.axis=0.5) #cex.axisΪ�����ֵĴ�С��as.characterת��Ϊ�ַ�
boxplot(log(geneLen[-idx]) ~ (genes$chrom[-idx]),cex.axis=0.5)   #���������ҵ�

#14 ���gene
idx=which.max(geneLen)
genes[idx, ]

#15 �����Ӹ���
nExon=genes$exonCount
hist(nExon,40, main="Number of exons", xlab="Exon count")

#16  ӵ����������ӵĻ���
idx=which.max(nExon)
genes[idx,1:9]   #��ʾidx�����е�1-9�У�10��11����Ϊ��������ʼ����ֹλ�ã�ӵ����������ӵ���������ʾ���ֶ�ҳ

#17  ���򳤶Ⱥ���������֮��Ĺ�ϵ
plot(geneLen, nExon, pch=".", xlab="Gene length", ylab="Number of exons")
cor(nExon, geneLen) ## weak correlation, need to run a test (lm or glm).
glm(nExon~geneLen)  #family�Ȳ��������������

#18 Ⱦɫ�峤����Ⱦɫ��ŵĶ�Ӧ����
chrlen=c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954)
names(chrlen)=paste("chr",c(1:22,"X","Y"),sep="")

#19  ÿ��Ⱦɫ����ÿ���׼���ԵĻ�����
chr_gene<-data.frame(chrlen, genes_num=rep(0,24))   #�������ݿ򣬻�����ĿΪ0
for (i in rownames(chr_gene)){
  chr_gene[rownames(chr_gene) %in% i,2]<-nrow(genes[genes$chrom %in% i,])    #ÿ��Ⱦɫ���ϵĻ�����Ŀ %in% �൱�� ==
}
gene_density_per_M<-round(chr_gene$genes_num/chr_gene$chrlen*10^6,2)  #round��������
genes_density<-cbind(chr_gene, gene_density_per_M)
genes_density
genes_density[order(-genes_density$gene_density_per_M),]  #����gene_density_per_M���� -Ϊ�Ӹߵ���


