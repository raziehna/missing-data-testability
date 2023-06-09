
# #################
#  Block-par example
# #################

rm(list = ls(all=TRUE))
cat("\f")

library(ggplot2)
library(reshape2)
library(MASS)
library(gam)
library(dplyr)

set.seed(0)

setwd("~/block_parallel/")

#################################################
# Functions
#################################################

source("dgp_bp.R")  # dgp()

expit <- function(x) 1/(1+exp(-x))

compute_odds <- function(dat, fmla, n, k, q){
  
  R = dat[, 1:K]
  
  numerator = rep(1, nrow(dat))
  denominator = rep(1, nrow(dat)) 
  
  # go over each Ri in 1, ..., K
  for (i in seq(1, K)) {
    
    # if Ri is not Rk or Rq
    if (i != k & i != q) {
      
      numerator = numerator * R[, i]
      denominator = denominator * R[, i]
    }
    
    # otherwise if Ri is Rk or Rq
    else {
      # drop any NAs _except_ for NAs in Xs_i
      # this line assumes that all Xs come after the R and X (hence 2*K)
      datRi = dat %>% filter(complete.cases(.[,-(2*K+i)]))
      model_Ri = gam(fmla[[i]], datRi, family = quasibinomial)
      proba_Ri = predict(model_Ri, dat, type = "response")
      numerator = numerator * (1 - R[, i])
      denominator = denominator * R[, i] * (1 - proba_Ri) / proba_Ri
    }
  }
  
  rho = mean(numerator)/mean(denominator)
  
  return(rho)
}

#################################################
# Initialization and DGP
#################################################

# args <- commandArgs(TRUE)
# iter = as.numeric(args[[1]]) # iterations
# min_N = as.numeric(args[[2]]) # sample size from min_N to max_N
# max_N = as.numeric(args[[3]]) #
# increments = as.numeric(args[[4]]) # increments in sample size 
# type_var = args[[5]] # "binary" vs "continuous"
# type_missing = args[[6]] # "BP" vs "NoSelf"


iter = 100
min_N = 1000
max_N = 10000
increments = 2000

type_var= "binary"
type_missing = "BP"

N = seq(min_N, max_N, by=increments)

# number of variables
K = 4

file_name = paste0(type_missing, "_", type_var, ".csv")

#################################################
# Analysis
#################################################

results = as.data.frame(matrix(0, ncol = 2, nrow=length(N)))
colnames(results) = c("n", "complete_rate")
odds_dat = as.data.frame(matrix(0, nrow = iter, ncol=length(N)))

for (m in seq(1, length(N))){

  n = N[m]
  cat("sample size: ", n, "\n")
  
  complete_rate = 0
  for (i in 1:iter){
  
    output = dgp(n, type_var, type_missing)
    fmla = output$fmla
    dat = output$dat
    attach(dat, warn.conflicts=FALSE)
    
    complete_rate = complete_rate + dim(dat[complete.cases(dat), ])[1]/n
    
    R = dat[, 1:K]
    
    k = 1
    q = 2
    
    rho = compute_odds(dat, fmla, n, k, q)
    
    odds_dat[i, m] = rho
    
} # end iterations

results[m, ] = c(n, complete_rate/iter)

}# end sample size

#################################################
# Resutls
#################################################


for (i in seq(1, length(N))){
  box = boxplot(odds_dat[, i], main=type_var)
  idx_out = which(odds_dat[, i] %in% box$out)
  odds_dat[idx_out, i] = median(odds_dat[, i])
}

write.csv(odds_dat, file_name, row.names = FALSE)


boxplot(odds_dat, main=type_var)
abline(h=1)

