######################
## Empirical Bayes  ##
######################

### parallelizing is performed here as explained in https://www.statistik.tu-dortmund.de/rechnerdoku/tutorial/Computecluster.html
p <- as.integer(Sys.getenv("PBS_ARRAYID"))
print(p)

# initializing the random number generator
library("rlecuyer")
.lec.SetPackageSeed(c(2250, 3376, 5270, 3920, 3760, 3186)) # rlecuyer equivalent for set.seed()
nstreams <- 50 # number of streams of random numbers
names <- paste("myrngstream",1:nstreams,sep="")
.lec.CreateStream(names) # generate the streams of random variables
.lec.CurrentStream(names[p]) # choose the p-th stream


library(drc)
library(LearnBayes)
library(mixtools)

load("../VPA dataset/VPA-Parameter-StatSign-BiolRel.RData")
par.b <- pars.considered[,1]
par.c <- pars.considered[,2]
par.d <- pars.considered[,3]
par.e <- pars.considered[,4]

##log-logistische function where e equals \tilde{e}
my.loglog <- function(x, b, c, d, e){
  return(c+(d-c)/(1+exp(b*(log(x)-e))))
}
conc <- c(rep(0, 6), rep(c(25, 150, 350, 450, 550, 800, 1000), each=3))

##Posterior in the Normal-Normal-Modell:
posterior.fct <- function(x, sigma2, e0, tau2){
  mean.post <- (tau2*x + sigma2*e0)/(tau2+sigma2)
  var.post <- (tau2*sigma2)/(tau2+sigma2)
  out <- c(mean.post, var.post)
  names(out) <- c("Post.Mean", "Post.Variance")
  return(out)
}


nr.datasets <- 20
nr.genes <- length(par.e)

emp.results <- lapply(1:nr.datasets, function(nn){
  print(nn)
  
  ##1. Step: Simulation of the entire dataset
  dataset <- t(sapply(1:nr.genes, function(gg){
    resp <- my.loglog(conc, par.b[gg], par.c[gg], par.d[gg], par.e[gg]) + rnorm(27, 0, 0.1)
    resp
  }))
  colnames(dataset) <- conc
  
  ##2. Step: Fit a 4pLL model to each gene
  pars.var <- sapply(1:nr.genes, function(gg){
    resp <- dataset[gg,]
    obj <- try(drm(resp ~ conc, fct=LL2.4()), silent=TRUE)
    
    if(inherits(obj, "try-error")) { out <- rep(NA, 8) }
    
    if(!inherits(obj, "try-error")){
      out <- c(coef(summary(obj))[, 1], (coef(summary(obj))[, 2])^2)
    }
    names(out) <- c("Par.b", "Par.c", "Par.d", "Par.e", "Var.b", "Var.c", "Var.d", "Var.e")
    out
  })
  
  ##3. Step: Estimation of prior mean and prior variance 
  ## a) based on empirical mean and empirical variance
  mean.prior <- mean(pars.var[4,], na.rm=TRUE)
  var.prior <- var(pars.var[4,], na.rm=TRUE)
  ## b) based on the empirical median and the empirical MAD
  median.prior <- median(pars.var[4,], na.rm=TRUE)
  mad.prior <- (stats::mad(pars.var[4,], na.rm=TRUE))^2
  ## c) based on the EM algorithm, as a mixture of 5 normal distributions
  gm <- try(normalmixEM(pars.var[4,], k=5, lambda=c(1/5, 1/5, 1/5, 1/5, 1/5), 
                        mu=c(6.2, 5, 8, 12, 8), sigma=c(0.3,0.5, 0.5, 1, 0.3), maxit=2000), silent=TRUE)
  
  if(inherits(gm, "try-error")) { 
    lambda.mix.prior <- NA
    mean.mix.prior <- NA
    var.mix.prior <- NA
    mix.posterior <- rep(NA, 5)
}
  
  if(!inherits(gm, "try-error")){
    lambda.mix.prior <- gm$lambda
    mean.mix.prior <- gm$mu
    var.mix.prior <- gm$sigma^2
    mix.posterior <- gm$posterior
  }

  ##4. Step: Apply "Bayes" for each gene individually for each version of the prior
  bayes.e <- sapply(1:nr.genes, function(gg){
    xx <- pars.var[4,gg]
    var.modell <- (pars.var[8,gg])
    
    out.bayes.ml <- posterior.fct(x = xx, sigma2 = var.modell, e0 = mean.prior, tau2 = var.prior)
    out.bayes.rob <- posterior.fct(x = xx, sigma2 = var.modell, e0 = median.prior, tau2 = mad.prior)
    
    ##Posterior for the mixed prior
    post.mix <- normal.normal.mix(probs=lambda.mix.prior, normalpar = rbind(c(mean.mix.prior[1], var.mix.prior[1]),
                                                                            c(mean.mix.prior[2], var.mix.prior[2]),
                                                                            c(mean.mix.prior[3], var.mix.prior[3]),
                                                                            c(mean.mix.prior[4], var.mix.prior[4]),
                                                                            c(mean.mix.prior[5], var.mix.prior[5])),
                                  data = c(xx, var.modell))
    post.mix.mean <- sum(post.mix$probs*post.mix$normalpar[,1])
    post.mix.var <- sum(post.mix$probs*(post.mix$normalpar[,2] + post.mix$normalpar[,1]^2 - post.mix.mean^2))

    if(inherits(gm, "try-error")) mix.posterior.out <- mix.posterior
    if(! inherits(gm, "try-error")) mix.posterior.out <- mix.posterior[gg,]
    
    out.bayes <- c(out.bayes.ml, out.bayes.rob, post.mix.mean, post.mix.var, 
                   post.mix$probs, post.mix$normalpar[,1], post.mix$normalpar[,2],
                   mix.posterior.out)
    names(out.bayes) <- c("Post.Mean.ML", "Post.Var.ML", "Post.Mean.Robust", "Post.Var.Robust",
                          "Post.Mean.Mix", "Post.Var.Mix", 
                          "Post.Mix.Lambda1", "Post.Mix.Lambda2", "Post.Mix.Lambda3", "Post.Mix.Lambda4", "Post.Mix.Lambda5", 
                          "Post.Mean.1", "Post.Mean.2", "Post.Mean.3", "Post.Mean.4", "Post.Mean.5", 
                          "Post.Var.1", "Post.Var.2", "Post.Var.3", "Post.Var.4", "Post.Var.5",
                          "Latent.Prior.1", "Latent.Prior.2", "Latent.Prior.3", "Latent.Prior.4", "Latent.Prior.5")
    out.bayes
  })
  
  ##5. Step: Output of the Table with Estimates, Standard Errors,
      ## Posterior Menas, Posterior Variances
  out.complete <- rbind(pars.var, bayes.e)
  colnames(out.complete) <- names(par.e)[1:nr.genes]
  out.complete
  
})

save(emp.results, file=paste0("SimResults/Results-Original-Nr", p, ".RData"))


