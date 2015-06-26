require(VGAM)
require(gtools)
source("Scripts/init.r")
source("Scripts/funcs.r")

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2)
{
	cat('\n\n')
	cat('\nCommand requires two arguments: sample number and number of components!\n')
	cat('\n\nThank you! Come again!\n\n')
	q();
}
sample   = as.numeric(args[1]); 
num.test = as.numeric(args[2]);
num.iter = 4000;
thin     = 100;



data.set = cbind(non[,sample],non[,sample]+ref[,sample]);
sample.name  = colnames(ref)[sample];

model.out = run.mcmc(num.test,allele.freq,data.set,num.iter,thin)

cat(paste("Model plotted to Figures/sample",sample.name,".pdf\n",sep=" "))
plot.figure(model.out,sample.name,data.set,allele.freq)

out.file = paste("Output/mcmc",sample.name,"RData",sep=".")
cat(paste("Full model writted to", out.file,"\n",sep=" "))
save(file = out.file,list=ls())

cat("\n\nComplete!\n")


