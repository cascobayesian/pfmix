ref.raw = read.table("Data/ref_pseudo.tab",header=TRUE,sep="\t")
non.raw = read.table("Data/non_pseudo.tab",header=TRUE,sep="\t")

ref = as.matrix(ref.raw);
non = as.matrix(non.raw);

frac = non/(ref+non);
allele.freq   = rowSums(non)/(rowSums(non+ref))
allele.freq.w = apply(frac,1,mean)

good = which((allele.freq<0.995)&(allele.freq>0.005))
ref = ref[good,];
non = non[good,];
frac = non/(ref+non);
allele.freq   = rowSums(non)/(rowSums(non+ref))

