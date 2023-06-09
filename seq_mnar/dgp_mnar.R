  
# ####################
#  DGP for Seq MNAR model
# ####################

dgp <- function(n, type_var, type_missing){
  
  if (type_var == "binary"){
    x_coeff = runif(9, -1, 1) 
    X1 = rbinom(n, 1, expit(x_coeff[1])) 
    X2 = rbinom(n, 1, expit(x_coeff[2] + x_coeff[3]*X1))
    X3 = rbinom(n, 1, expit(x_coeff[4] + x_coeff[5]*X1 + x_coeff[6]*X2))
    X4 = rbinom(n, 1, expit(x_coeff[7] + x_coeff[8]*X1 + x_coeff[9]*X2))
  }
  
  if (type_var == "continuous"){
    mu = rep(0, 4)
    sigma = matrix(c(1, 0.75, 0.5, 0.25, 0.75, 1, 0.75, 0.5, 0.5, 0.75, 1, 0.75, 0.25, 0.5, 0.75, 1), nrow = 4)
    X = mvrnorm(n, mu, sigma)
    X1 = X[, 1]
    X2 = X[, 2]
    X3 = X[, 3]
    X4 = X[, 4]
  }
  
  a_r1 = runif(4, -1, 1) # about 6% complete cases
  a_r2 = runif(5, -1, 1)
  a_r3 = runif(6, -1, 1)
  a_r4 = runif(7, -1, 1)
  
  # a_r1 = runif(4, -0.5, 1.5) # about 30% complete cases
  # a_r2 = runif(5, -0.5, 1.5)
  # a_r3 = runif(6, -0.5, 1.5)
  # a_r4 = runif(7, -0.5, 1.5)
  
  # a_r1 = runif(4, 0, 2) # about 48% complete cases
  # a_r2 = runif(5, 0, 2)
  # a_r3 = runif(6, 0, 2)
  # a_r4 = runif(7, 0, 2)

  R1 = rbinom(n, 1, expit(a_r1[1] + a_r1[2]*X2 + a_r1[3]*X3 + a_r1[4]*X4))
  
  if (type_missing=="SeqMNAR"){
    R2 = rbinom(n, 1, expit(a_r2[1] + a_r2[2]*R1 + a_r2[3]*X3 + a_r2[4]*X4))
    R3 = rbinom(n, 1, expit(a_r3[1] + a_r3[2]*R1 + a_r3[3]*R2 + a_r3[4]*X4))
    R4 = rbinom(n, 1, expit(a_r4[1] + a_r4[2]*R1 + a_r4[3]*R2 + a_r4[4]*R3))
  }
  
  if (type_missing=="Permutation"){
    R2 = rbinom(n, 1, expit(a_r2[1] + a_r2[2]*R1 + a_r2[3]*X3 + a_r2[4]*X4 + a_r2[5]*X1*R1))
    R3 = rbinom(n, 1, expit(a_r3[1] + a_r3[2]*R1 + a_r3[3]*R2 + a_r3[4]*X4 + a_r3[5]*X1*R1 + a_r3[6]*X2*R2))
    R4 = rbinom(n, 1, expit(a_r4[1] + a_r4[2]*R1 + a_r4[3]*R2 + a_r4[4]*R3 + a_r4[5]*X1*R1 + a_r4[6]*X2*R2 + a_r4[7]*X3*R3))  
  }
  
  X1s = X1
  X1s[R1==0] = NA
  
  X2s = X2
  X2s[R2==0] = NA
  
  X3s = X3
  X3s[R3==0] = NA
  
  X4s = X4
  X4s[R4==0] = NA
  
  dat = data.frame(R1, R2, R3, R4, 
                   X1, X2, X3, X4, 
                   X1s, X2s, X3s, X4s, 
                   XR1 = X1*R1, XR2 = X2*R2, XR3 = X3*R3)  
  
  # formulas
  fmla_null = list(as.formula(R1 ~ X2 + X3 + X4), 
                   as.formula(R2 ~ R1 + X3 + X4), 
                   as.formula(R3 ~ R1 + R2 + X4), 
                   as.formula(R4 ~ R1 + R2 + R3))
  
  fmla_alt = list(as.formula(R1 ~ X2 + X3 + X4), 
                  as.formula(R2 ~ R1 + X3 + X4 + XR1), 
                  as.formula(R3 ~ R1 + R2 + X4 + XR1 + XR2), 
                  as.formula(R4 ~ R1 + R2 + R3 + XR1 + XR2 + XR3))
  
  return(list(dat=dat, 
              fmla_null=fmla_null, 
              fmla_alt=fmla_alt))
}

