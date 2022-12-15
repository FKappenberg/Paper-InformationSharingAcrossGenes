library(ggplot2)
library(scales)
#load("Results-EmpiricalBayes-MixedPrior-Complete-5Priors.RData")
load("Results-EstPerGene.RData")
load("Results-BaPerGeneML.RData")
load("Results-BaPerGeneRob.RData")
load("Results-BaPerGeneMix.RData")
load("Results-PostvarMLPerGene.RData")

load("VPA datasets/VPA-Parameter-StatSign-BiolRel.RData")
par.mat <- pars.considered


inds.missing <- which(apply(ba.per.gene.rob, 2, function(x) sum(is.na(x))) > 200 | 
                        apply(ba.per.gene.ml, 2, function(x) sum(is.na(x))) > 200 | 
                        apply(ba.per.gene.mix, 2, function(x) sum(is.na(x))) > 200)


##
# Calculation of MSE for the direct estimation
mse.est <- sapply(((1:7191)[-inds.missing]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.ml[,gene]))
  mse <- mean((est.per.gene[ind.considered,gene] - e.truth)^2)
  mse
})

# Calculation of MSE for the Bayes approach with ML prior
mse.bayes.ml <- sapply(((1:7191)[-inds.missing]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.ml[,gene]))
  mse <- mean((ba.per.gene.ml[ind.considered,gene] - e.truth)^2)
  mse
})

# Calculation of MSE for the Bayes approach with robust prior
mse.bayes.rob <- sapply(((1:7191)[-inds.missing]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.rob[,gene]))
  mse <- mean((ba.per.gene.rob[ind.considered,gene] - e.truth)^2, na.rm=TRUE)
  mse
})

# Calculation of MSE for the Bayes approach with mixed prior
mse.bayes.mix <- sapply(((1:7191)[-inds.missing]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.mix[,gene]))
  mse <- mean((ba.per.gene.mix[ind.considered,gene] - e.truth)^2, na.rm=TRUE)
  mse
})





#Second part: Consider only those gnees for which the originally underlying parameter e is smaller than 6.91
inds.missing.new1 <- which(apply(ba.per.gene.rob, 2, function(x) sum(is.na(x))) > 200 | 
                             apply(ba.per.gene.ml, 2, function(x) sum(is.na(x))) > 200 | 
                             apply(ba.per.gene.mix, 2, function(x) sum(is.na(x))) > 200 |
                             par.mat[,4] > 6.91)

##
# Calculation of MSE for the direct estimation
mse.est.new1 <- sapply(((1:7191)[-inds.missing.new1]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.ml[,gene]))
  mse <- mean((est.per.gene[ind.considered,gene] - e.truth)^2)
  mse
})

# Calculation of MSE for the Bayes approach with ML prior
mse.bayes.ml.new1 <- sapply(((1:7191)[-inds.missing.new1]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.ml[,gene]))
  mse <- mean((ba.per.gene.ml[ind.considered,gene] - e.truth)^2)
  mse
})

# Calculation of MSE for the Bayes approach with robust prior
mse.bayes.rob.new1 <- sapply(((1:7191)[-inds.missing.new1]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.rob[,gene]))
  mse <- mean((ba.per.gene.rob[ind.considered,gene] - e.truth)^2, na.rm=TRUE)
  mse
})

# Calculation of MSE for the Bayes approach with mixed prior
mse.bayes.mix.new1 <- sapply(((1:7191)[-inds.missing.new1]), function(gene){
  e.truth <- par.mat[gene, 4]
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.mix[,gene]))
  mse <- mean((ba.per.gene.mix[ind.considered,gene] - e.truth)^2, na.rm=TRUE)
  mse
})


# Coverage probabilities for the confidence intervals based on the direct estimation
cp.est <- sapply(((1:7191)[-inds.missing]), function(gene){
  e.truth <- par.mat[gene, 4]
  
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.ml[,gene]))
  
  quant <- qt(p=0.975, df=(27-4))
  lower.boundary <- pmax(0, est.per.gene[ind.considered, gene] - quant* sqrt(Var.par.e.per.gene[ind.considered, gene]))
  upper.boundary <- est.per.gene[ind.considered, gene] + quant* sqrt(Var.par.e.per.gene[ind.considered, gene])
  covered <- (e.truth > lower.boundary & e.truth < upper.boundary)
  
  cp <- mean(covered)
  cp
})

# Coverage probabilities for the credible intervals based on the direct estimation
cp.bayes.ml <- sapply(((1:7191)[-inds.missing]), function(gene){
  e.truth <- par.mat[gene, 4]
  
  ind.considered <- which(!is.na(est.per.gene[,gene]) & !is.na(ba.per.gene.ml[,gene]))
  
  lower.boundary <- qnorm(0.025, mean=ba.per.gene.ml[ind.considered, gene], sd = sqrt(Post.var.ML.per.gene[ind.considered, gene]))
  upper.boundary <- qnorm(0.975, mean=ba.per.gene.ml[ind.considered, gene], sd = sqrt(Post.var.ML.per.gene[ind.considered, gene]))
  covered <- (e.truth > lower.boundary & e.truth < upper.boundary)
  
  cp <- mean(covered)
  cp
})
