
# #################
#  MAR example
# #################

rm(list = ls(all=TRUE))
cat("\f")

library(ggplot2)
library(reshape2)
library(MASS)
library(gam)
library(dplyr)

set.seed(0)

setwd("~/seq_mar/")

#################################################
# Functions
#################################################

source("dgp_mar.R")  # dgp()

expit <- function(x) 1/(1+exp(-x))


#################################################
# Initialization and DGP
#################################################

# args <- commandArgs(TRUE)
# iter = as.numeric(args[[1]]) # iterations
# min_N = as.numeric(args[[2]]) # sample size from min_N to max_N
# max_N = as.numeric(args[[3]]) #
# increments = as.numeric(args[[4]]) # increments in sample size
# type_var = args[[5]] # "binary" vs "continuous"
# type_missing = args[[6]] # "SeqMAR" vs "Permutation"


iter = 100
min_N = 500
max_N = 15000
increments = 500

type_var= "binary"
type_missing = "SeqMAR"

N = seq(min_N, max_N, by=increments)

# number of variables
K = 4

file_name = paste0(type_missing, "_", type_var, ".csv")

#################################################
# Analysis 
#################################################

results = as.data.frame(matrix(0, ncol = 3, nrow=length(N)))
colnames(results) = c("n", "accept_rate", "complete_rate") 

for (m in seq(1, length(N))){

n = N[m]
cat("sample size: ", n, "\n")

# pvals = as.data.frame(matrix(0, ncol = K-1, nrow=iter))
# colnames(pvals) = c("pval_R1", "pval_R2", "pval_R3")
count = iter
complete_rate = 0
for (i in 1:iter){
  
output = dgp(n, type_var, type_missing)
fmla_null = output$fmla_null
fmla_alt = output$fmla_alt
dat = output$dat
attach(dat, warn.conflicts=FALSE)

complete_rate = complete_rate + dim(dat[complete.cases(dat), ])[1]/n

R = dat[, 1:K] 
# X = dat[, (K+1):(2*K)] 
# Xs = dat[, (2*K+1):(3*K)]

model_base = glm(fmla_null[[K]], dat, family = binomial) 
W = predict(model_base, dat, type = "response") # p(RK = 1 | R_{\prec K}, Xs_{\prec K})

for (k in seq(K-1, 1)){

  # ++++++++++++++++++
  # null 
  # ++++++++++++++++++
  model_null = glm(fmla_null[[k]], dat, family = binomial) 
  Wk_null = predict(model_null, dat, type = "response") # p(Rk = 1 | R_{\prec k}, XR_{\prec k}) 

  # ++++++++++++++++++
  # alternative 
  # ++++++++++++++++++
  idx = seq(1, n) # R_{\succ k} == 1
  for (j in seq(k+1, K)){  
    idx_t = which(R[,j] == 1)
    idx = intersect(idx, idx_t)
  }
  
  model_alt = gam(fmla_alt[[k]], dat[idx, ], family = quasibinomial, weights = W[idx])
  Wk_alt = predict(model_alt, dat, type = "response") # p(Rk = 1 | R_{\prec k}, XR_{\prec k}, X_{\succ k})

  # ++++++++++++++++++
  # likelihood ratio test 
  # ++++++++++++++++++
  idx_Rk = which(R[, k] == 1) # Rk == 1

  Wk_null_natural = Wk_null  
  Wk_null_natural[-idx_Rk] = 1 - Wk_null_natural[-idx_Rk] 
  
  Wk_alt_natural = Wk_alt  
  Wk_alt_natural[-idx_Rk] = 1 - Wk_alt_natural[-idx_Rk] 
  
  rho = -2*sum((1/W[idx])*log(Wk_null_natural[idx]/Wk_alt_natural[idx]))
  pval = pchisq(rho, df = K-k, lower.tail=F)
  
  # pvals[i, k] = pval
  
  if (pval < 0.05){
    count = count - 1
    break
  }

  W = W*Wk_null 
  
  # ++++++++++++++++++
  # sanity check 
  # ++++++++++++++++++
  
  # -2*sum(log(Wk_null_natural))
  # model_null$deviance
  # 
  # model_alt$deviance
  # -2*sum(log(Wk_alt_natural[idx]))
  # -2*sum((1/W)*log(Wk_alt_natural))
  # -2*sum((1/W[idx])*log(Wk_alt_natural[idx]))
  # 
  # rho_2 = model_null$deviance - model_alt$deviance
  # pchisq(rho_2, df = K-k, lower.tail=F)
  
} # tests in a single model

} # end iterations

results[m, ] = c(n, count/iter, complete_rate/iter)

}# end sample size

#################################################
# Resutls 
#################################################

write.csv(results, file_name, row.names = FALSE)

plot(N, results$accept_rate, xlab = "Sample size (n)", ylab = "Acceptance rate", ylim=c(0, 1), type="b", main = type_missing)


