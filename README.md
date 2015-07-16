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

In order for the algorithm to run the appropriate data files need to be placed in the Data folder. The files give the number of reference and non-reference read counts at each site at each position in two separate files. They need to be in tab-delimited format with no header, the columns giving the samples and the rows the variants (SNPs). They also need to be the same files names as appear at the top of the init.r file. The algorithm operates on the assumption that the set of samples is sufficient to infer a non-zero local population level allele frequency. If that is not the case, strange results will occur, mostly likely in the form of an error. You can ensure that this is the case by undertaking the steps in clear.r

The scripts can be run most easily by going to the terminal and typing the command:

Rscript Scripts/run.r i j 

where i is the sample number you want and j is the total number of components to test. A reasonable value for j is 5 and, for time reasons, you probably don't want j to be more than 8. 

Setting aside the cleaning issue, these scripts should produce the results required. The output will be saved in a *.RData file with the name of the sample in the Output folder. A figure file will also be produced with the sample name in the Figures folder, showing the raw data and the inferred model on a single plot. For each sample, you should check that the output seems reasonable since in a small number of cases the procedure produces errant inferences. These can usually be resolved by fixing the number of strains inferred. 

If you have any questions, please feel free to contact me: jobrien@bowdoin.edu.
  
  
  
