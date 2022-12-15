#################################################
##    Vorbereitungen für die                   ##
## Auswertung der empirischen Bayes-Methode    ##
##        für den originalen Datensatz         ##
#################################################

setwd("../02 Cluster Results/")

##rownames of one "set" in emp.results:
# [1] "Par.b"            "Par.c"            "Par.d"            "Par.e"            "Var.b"            "Var.c"            "Var.d"           
# [8] "Var.e"            "Post.Mean.ML"     "Post.Var.ML"      "Post.Mean.Robust" "Post.Var.Robust"  "Post.Mean.Mix"    "Post.Var.Mix"    
# [15] "Post.Mix.Lambda1" "Post.Mix.Lambda2" "Post.Mix.Lambda3" "Post.Mix.Lambda4" "Post.Mix.Lambda5" "Post.Mean.1"      "Post.Mean.2"     
# [22] "Post.Mean.3"      "Post.Mean.4"      "Post.Mean.5"      "Post.Var.1"       "Post.Var.2"       "Post.Var.3"       "Post.Var.4"      
# [29] "Post.Var.5"       "Latent.Prior.1"   "Latent.Prior.2"   "Latent.Prior.3"   "Latent.Prior.4"   "Latent.Prior.5"     

##goal: extract the estimates and write them in a 1000x7191 matrix
estimate.per.gene <- lapply(1:50, function(i){
  ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
  load(list.files()[i])
  ##12x7191 matrix
  est <- t(sapply(emp.results, function(set){
    set["Par.e",]
  }))
  est
})
##bind all matrices together to the final large matrix
est.per.gene <- do.call(rbind, estimate.per.gene)


##Bayes estimate via ML
bayes.per.gene <- lapply(1:50, function(i){
  load(list.files()[i])
  est <- t(sapply(emp.results, function(set){
    set["Post.Mean.ML",]
  }))
  est
})
ba.per.gene.ml <- do.call(rbind, bayes.per.gene)


##Bayes estimate via robust estimation
bayes.per.gene <- lapply(1:50, function(i){
  load(list.files()[i])
  est <- t(sapply(emp.results, function(set){
    set["Post.Mean.Robust",]
  }))
  est
})
ba.per.gene.rob <- do.call(rbind, bayes.per.gene)


##Bayes estimate via mixed prior
bayes.per.gene <- lapply(1:50, function(i){
  load(list.files()[i])
  est <- t(sapply(emp.results, function(set){
    set["Post.Mean.Mix",]
  }))
  est
})
ba.per.gene.mix <- do.call(rbind, bayes.per.gene)

# 
# ##Squared standard error of parameter e
# Var.Pare.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set[8,]
#   }))
#   est
# })
# Var.par.e.per.gene <- do.call(rbind, Var.Pare.per.gene)
# 
# 


##Posterior variance via ML estimation
PostVar.ML.per.gene <- lapply(1:50, function(i){
  load(list.files()[i])
  est <- t(sapply(emp.results, function(set){
    set["Post.Var.ML",]
  }))
  est
})
Post.var.ML.per.gene <- do.call(rbind, PostVar.ML.per.gene)
# 
# 

# ##First posterior probability
# Lambda1.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Post.Mix.Lambda1",]
#   }))
#   est
# })
# Lambda1.per.gene <- do.call(rbind, Lambda1.per.gene)
# 
# ##Second posterior probability
# Lambda2.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Post.Mix.Lambda2",]
#   }))
#   est
# })
# Lambda2.per.gene <- do.call(rbind, Lambda2.per.gene)
# 
# ##Third posterior probability
# Lambda3.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Post.Mix.Lambda3",]
#   }))
#   est
# })
# Lambda3.per.gene <- do.call(rbind, Lambda3.per.gene)
# 
# ##Fourth posterior probability
# Lambda4.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Post.Mix.Lambda4",]
#   }))
#   est
# })
# Lambda4.per.gene <- do.call(rbind, Lambda4.per.gene)
# 
# 
# ##Fifth posterior probability
# Lambda5.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Post.Mix.Lambda5",]
#   }))
#   est
# })
# Lambda5.per.gene <- do.call(rbind, Lambda5.per.gene)

# 
# ##First Latent Prior Probability
# Prior1.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Latent.Prior.1",]
#   }))
#   est
# })
# Prior1.per.gene <- do.call(rbind, Prior1.per.gene)
# 
# 
# 
# ##Second Latent Prior Probability
# Prior2.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Latent.Prior.2",]
#   }))
#   est
# })
# Prior2.per.gene <- do.call(rbind, Prior2.per.gene)
# 
# 
# 
# ##Third Latent Prior Probability
# Prior3.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Latent.Prior.3",]
#   }))
#   est
# })
# Prior3.per.gene <- do.call(rbind, Prior3.per.gene)
# 
# 
# 
# ##Fourth Latent Prior Probability
# Prior4.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Latent.Prior.4",]
#   }))
#   est
# })
# Prior4.per.gene <- do.call(rbind, Prior4.per.gene)
# 
# 
# 
# ##Fifth Latent Prior Probability
# Prior5.per.gene <- lapply(1:50, function(i){
#   ##load the i-th result file: list with 20 elements, each is a 4x7191 matrix
#   load(list.files()[i])
#   ##20x7191 matrix
#   est <- t(sapply(emp.results, function(set){
#     set["Latent.Prior.5",]
#   }))
#   est
# })
# Prior5.per.gene <- do.call(rbind, Prior5.per.gene)


#save(est.per.gene, ba.per.gene.ml, ba.per.gene.rob, ba.per.gene.mix,
#     #Var.par.e.per.gene
#     Post.var.ML.per.gene,
#     #Lambda1.per.gene, Lambda2.per.gene, Lambda3.per.gene, Lambda4.per.gene, Lambda5.per.gene,
#     #Prior1.per.gene, Prior2.per.gene, Prior3.per.gene, Prior4.per.gene, Prior5.per.gene,
#     file="../Results-EmpiricalBayes-MixedPrior-Complete-5Priors.RData")

save(est.per.gene, file="../Results-EstPerGene.RData")
save(ba.per.gene.ml, file="../Results-BaPerGeneML.RData")
save(ba.per.gene.rob, file="../Results-BaPerGeneRob.RData")
save(ba.per.gene.mix, file="../Results-BaPerGeneMix.RData")
save(Post.var.ML.per.gene, file="../Results-PostvarMLPerGene.RData")











