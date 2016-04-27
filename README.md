# pfmix
R code for visualizing and inferring the structure of Plasmodium falciparum mixture

This README covers two distinct but related functions that you can find in the pfmix package: how to infer the structure of strain mixture using a Bayesian mixture model, and how to infer the inbreeding coefficient for mixed samples using one of the four estimators. 

Both of these functionalities applies to the same sort of data: whole-genome sequence read counts collected from clinical samples. These data are going to be best if there is no PCR performed (direct sequencing) but should work tolerably even if PCR was used. 

To get started with either calculation, first download the package and install it. There's some helpful advice elsewhere on the interwebs for this: [here](http://stackoverflow.com/questions/1474081/how-do-i-install-an-r-package-from-source) and [here](https://cran.r-project.org/doc/manuals/r-devel/R-admin.html). Then load the library using `library(pfmix)`.

Both calculations rely on data sets that assume a biallelic structure, so that data is naturally formulated as two matrices: one with all the reference read counts, and one with all the non-reference read counts. Once you've loaded the library, you can load some example data by typing `data(pf_data)`. This will load two matrices: `ref` and `non`. Both include the first 199 SNPs from 344 Ghanaian samples (rows are SNPs, columns are samples). These data derive from the genuinely excellent [PF3K resource](https://www.malariagen.net/data/pf3k-3).

In both use cases, it's going be useful to calculate the allele frequency for the SNPs. You can do that with this code:

`library(abind)`

`allele.freq <- calc.mle.p(abind(ref,non,along=3)`

### Mixture model

This model works sample-by-sample, so the first thing you'll need is to pull out the data for a sample of interest, say `s`. You are also going to want to choose a maximum number of strains to check (Eight is a reasonable value, but smaller also works).

You will also need to set the number of iterations and the amount of thinning in the chain. I'll use 4000 and 10, respectively. 

Here's the code for all of this:

`num.iter <- 4000; thin <- 10`

`sample.data <- data.set = cbind(non[,s],non[,s]+ref[,s])`

'out <- run.mcmc(num.test = 8,allele.freq,sample.data,num.iter,thin)`

If you'd like to plot the resulting model over the data to check for any suspicious fit or to jazz up a paper:

`plot.figure(model.out,sample.name,data.set,allele.freq)`

where `sample.name' is whatever name you want for the file. (It will try to write it to a Figures folder, so if there is no such folder below your working directory, you will get an error!) 

### Inbreeding coefficients

You can also use a similar framework to estimate the inbreeding coefficients for samples. This requires that you have an ambient population that the sample is drawn from, which you most likely already have if you've calculated `allele.freq` above. 

Other than that the framework is pretty similar: calculations are on sample-by-sample basis so it's necessary to pull out the data for each sample first. In this case, the functions require the within-sample allele frequency rather than the full data set, so we first transform them:

`wsaf <- non[,s]/(non[,s]+ref[,s])`

Then, any of the frequentist estimators are easy. Here's the direct estimator:

`f <- calc.f.fst(wasf,allele.freq)`

The other two are the initial (`calc.f.ini`) and the regressed (`calc.f.reg`), which take the same inputs. There is also a fully Bayesian approach that will simultaneously estimate all of the allele frequencies and the inbreeding coefficients. To do this, you'll have to refactor the matrices a bit into an array by:

`data.set = abind(ref,non,along=3)`

Then, just run the MCMC: 

`mcmc.run <- mcmc.fstat(data.set,num.iter,thin)`








