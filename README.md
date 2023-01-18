# Paper-InformationSharingAcrossGenes

In toxicological concentration-response studies, a frequent goal is the determination of an `alert concentration’, i.e. the lowest concentration where a notable change in the response in comparison to the control is observed. 
In high-throughput gene expression experiments, e.g. based on microarray or RNA-seq technology, concentration-response profiles can be measured for thousands of genes simultaneously.
One approach for determining the alert concentration is given by fitting a parametric model to the data which allows interpolation between the tested concentrations.

It is well known that the quality of a model fit improves with the number of measured data points.
However, adding new replicates for existing concentrations or even several replicates for new concentrations is time-consuming and expensive.
Here, we propose an empirical Bayes approach to information sharing across genes, where in essence a weighted mean of the individual estimate for one specific parameter of a fitted model and the mean of all estimates of the entire set of genes is calculated as a result.
Results of a controlled plasmode simulation study show that for many genes a notable improvement in terms of the mean squared error between estimate and true underlying value of the parameter can be observed.
However, for some genes, the MSE increases, and this cannot be prevented by using a more sophisticated prior distribution in the Bayesian approach.

In this GitHub repository, the code needed to reproduce the simulation studies as well as the resulting datasets is stored.