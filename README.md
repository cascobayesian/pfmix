# pfmix
### Inferring the structure of Plasmodium falciparum mixture

This forms a basic README for the pfmix code and workflows. Hopefully, it will provide sufficient information 
(or, more precisely pointer to information) so that any individual with basic R knowledge could retool these scripts
to implement the model on their own machine.

In general, once the data is in the right format and has been cleaned for a few different attributes, the model
should run easily and (relatively) quickly. At least that's the theory. 

The key document for getting the data into the right format, doing initial qc checks, and running the model is:
  * manual.pdf
  
The key scripts are:
  * init.r
  * clean.r
  * run.model.r
  * funcs.r 

Soon, all of this will be turned into an R package.

If you have any questions, please feel free to contact me: jobrien@bowdoin.edu.
  
  
  
