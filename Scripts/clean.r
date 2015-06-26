ref.raw = read.table("Data/new_ref.txt",header=TRUE,sep=",")
non.raw = read.table("Data/new_non.txt",header=TRUE,sep=",")

pos = as.matrix(ref.raw[,c(1:4)])

ref = as.matrix(ref.raw[,-c(1:4)])
non = as.matrix(non.raw[,-c(1:4)])

total = ref+non
total[(ref+non)<20] = NA;

cs = colSums(is.na(total))

pdf('Figures/supplementary_1.pdf')
plot(sort(cs),pch=20)
abline(12000,0,col="blue",lwd=2,lty=2);

bad.sams = which(cs > 12000)

ref = ref[,-bad.sams]
non = non[,-bad.sams]

total = ref+non;
total[(ref+non)<20] = NA

rs = rowSums(is.na(total))

bad.snps = which(rs > 0)

plot(sort(rs),pch=20,xlab = "Sample (sorted by missingness)",ylab="Number of SNPs missing (counts<20)")
dev.off()
ref = ref[-bad.snps,];
non = non[-bad.snps,];
pos = pos[-bad.snps,];

variable =which((rowSums(ref)>20)&(rowSums(non)>20))

ref = ref[variable,];
non = non[variable,];
pos = pos[variable,];

allele.freq = rowSums(non)/rowSums(ref+non)

low.freq = which((allele.freq<0.015)|(allele.freq>0.985))

ref = ref[-low.freq,]
non = non[-low.freq,]
pos = pos[-low.freq,];

frac = non/(ref+non)
num.samples = dim(frac)[2];

pw.mat = array(0,dim=c(num.samples,num.samples))
for(i in 1:(num.samples-1))
{
	#print(i)
	for(j in (i+1):num.samples)
	{
		pw.mat[i,j] = mean((frac[,i]-frac[,j])^2)
		pw.mat[j,i] = pw.mat[i,j]
	}
}

par(mfrow=c(1,3))
pca = princomp(pw.mat)

plot(pca$loadings[,1],pca$loadings[,2],col="blue",pch=20)
plot(pca$loadings[,1],pca$loadings[,3],col="blue",pch=20)
plot(pca$loadings[,2],pca$loadings[,3],col="blue",pch=20)

bad = which(pca$loadings[,1]>-0.07)
## noticeably less mixed than average 

ref = ref[,-bad]
non = non[,-bad];

pw.mat = pw.mat[-bad,-bad];
library(ape)
pdf('Figures/supplementary_2.pdf')
par(mfrow=c(2,2))
pca = princomp(pw.mat)

plot(pca$loadings[,1],pca$loadings[,2],col="blue",pch=20,xlab="PCA 1",ylab="PCA 2")
plot(pca$loadings[,1],pca$loadings[,3],col="blue",pch=20,xlab = "PCA 1", ylab = "PCA 3")
plot(pca$loadings[,2],pca$loadings[,3],col="blue",pch=20, xlab="PCA 2",ylab = "PCA 3")
plot(nj(pw.mat),"unrooted",show.tip.label=FALSE)

dev.off();

allele.freq = rowSums(non)/rowSums(ref+non)

low.freq = which((allele.freq<0.015)|(allele.freq>0.985))

ref = ref[-low.freq,]
non = non[-low.freq,]
pos = pos[-low.freq,]

frac = non/(ref+non)

write.table(x= ref,file = "Data/ref.clean.csv",col.names=TRUE,row.names=FALSE,sep=",")
write.table(x= non,file = "Data/non.clean.csv",col.names=TRUE,row.names=FALSE,sep=",")
write.table(x=colnames(ref),file="Data/samples.list",col.names=FALSE,row.names=FALSE,sep=",")	
