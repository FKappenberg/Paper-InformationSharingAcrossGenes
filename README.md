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


## VPA-Datasets

The underlying data, originally published in Krug et al. (2013, https://link.springer.com/article/10.1007/s00204-012-0967-3), are stored here in the way that they are used in the publication.

- VPA.RData: Appropriately pre-processed expression values for all 54675 probe sets and all 27 samples of the gene expression data.
- VPA-Parameter-StatSign-BiolRel.RData: parameters of 4pLL models (using the parameterization $\tilde{e} = \log(e)$) fitted to those 7191 probe sets that fulfilled the criteria for statistical significance and biological relevance

## Simulation and Preparation

- EmpiricalBayes-Simulation-MixedPrior-CompleteDatasets.R: contains the code for the simulation procedure. The simulation results are stored in the files in the folder "SimResults"
- Preparation-Analysis-MixedPrior.R: Summarizing the simulation results to datasets containing e.g. all estimates for the Bayes, ML approach to be analyzed later

## Main directory
 
The following results files are created by the code in "Preparaion-Analysis-MixedPrior.R" and contain the following:

- Results-BaPerGeneML.RData: Estimates for parameter $\tilde{e}$ as obtained by the Bayes procedure with the ML prior for all 7191 genes and all 1000 simulation runs
- Results-BaPerGeneMix.RData: Estimates for parameter $\tilde{e}$ as obtained by the Bayes procedure with the mixed prior for all 7191 genes and all 1000 simulation runs
- Results-BaPerGeneRob.RData: Estimates for parameter $\tilde{e}$ as obtained by the Bayes procedure with the robust prior for all 7191 genes and all 1000 simulation runs
- Results-EstPerGene.RData: Estimates for parameter $\tilde{e}$ as obtained by directly estimating it
- Results-PostVarMLPerGene.RData: Posterior variance of parameter $\tilde{e}$ as obtained by the Bayes procedure with the ML prior for all 7191 genes and all 1000 simulation runs

The final analysis of the resulting MSEs is performed in "AnalysisInformationSharing.R"
