# pfmix
### Inferring the structure of Plasmodium falciparum mixture from WGS data

This forms a basic README for the pfmix code and workflows. Hopefully, it will provide sufficient information 
(or, more precisely pointer to information) so that any individual with basic R knowledge could retool these scripts
to implement the model on their own machine.

In general, once the data is in the right format and has been cleaned for a few different attributes, the model
should run easily and (relatively) quickly. At least that's the theory. 

The key scripts are:
  * init.r
  * clean.r
  * run.model.r
  * funcs.r 

In order for the algorithm to run the data files need to be placed in the Data folder. The files give the number of reference and non-reference read counts at each site at each position in two separate files. They need to be in tab-delimited format with no header, the columns giving the samples and the rows the variants (SNPs). The algorithm operates on the assumption that the set of samples is sufficient to infer a local population level allele frequency. If that is not the case, strange results will result. 

The scripts can be run most easily by going to the terminal and running the command:

Rscript Scripts/run.r i j 

where i is the sample number you want and j is the total number of components to test. WARNING: you probably don't want j to be more than 8. 

Soon, all of this will be turned into an R package.

If you have any questions, please feel free to contact me: jobrien@bowdoin.edu.
  
  
  
