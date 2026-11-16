### test if numerics suggest globally attracting biologically relevant single species equilibrium
### Author: Joe Brennan

### perennial competition

library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(1245)

germ.data = read.table('germMat.txt', header = T, sep = '')

comptogerm = read.table('perTrans.txt', header = F, sep='')
colnames(comptogerm) = c('beta_S', 'beta_A', 'beta_M', 'v')

#### SP RELATED PARAMETERS

SP.data = read.table('SP.txt', header = F, sep = '')
colnames(SP.data) = c('r', 'p', 'alpha_{SP,SP}', 'alpha_{SP,EG}', 'alpha_{SP,BH}', 'alpha_{SP,BD}', 'alpha_{SP,AB}', 'leaf damage impact')

SP.lambda = SP.data$r * (1 - SP.data$p) / SP.data$p
alpha_SP = SP.data[,3:7]
g_SP = germ.data$SP

omega_SP = read.table('SPAdultSurv.txt', header = F)[,1]

#### EG RELATED PARAMETERS

EG.data = read.table('EG.txt', header = F, sep = '')
colnames(EG.data) = c('r', 'p', 'alpha_{EG,SP}', 'alpha_{EG,EG}', 'alpha_{EG,BH}', 'alpha_{EG,BD}', 'alpha_{EG,AB}', 'leaf damage impact')

EG.lambda = EG.data$r * (1 - EG.data$p) / EG.data$p
alpha_EG = EG.data[,3:7]
g_EG = germ.data$EG

omega_EG = read.table('EGAdultSurv.txt', header = F)[,1]

### perennial competition simulator: SP vs EG

simulate.perennial.comp = function(times = 2000, k = 1, n.init = c(0,0), p.init = c(0,0)){
  
  #### species names
  spp = c("SP", "EG")
  
  #### demographic parameters
  g = list(
    SP = g_SP,
    EG = g_EG
  )
  
  v = list(
    SP = comptogerm$v[k],
    EG = comptogerm$v[k]
  )
  
  lambda = list(
    SP = SP.lambda[k],
    EG = EG.lambda[k]
  )
  
  omega = list(
    SP = omega_SP[k],
    EG = omega_EG[k]
  )
  
  #### alpha competition matrix
  alpha =
    matrix(
      c(alpha_SPSP=alpha_SP[k, "alpha_{SP,SP}"],
        alpha_SPEG=alpha_SP[k, "alpha_{SP,EG}"],
        alpha_EGSP=alpha_EG[k, "alpha_{EG,SP}"],
        alpha_EGEG=alpha_EG[k, "alpha_{EG,EG}"]),
      nrow=2,
      byrow=T
    )
  colnames(alpha)<- rownames(alpha) <- c('SP','EG')
  
  #### beta parameters
  beta = list(
    s = comptogerm$beta_S[k],
    m = comptogerm$beta_M[k]
  )
  
  ### storage vectors
  
  n = matrix(NA, nrow = times + 1, ncol = 2)
  p = matrix(NA, nrow = times + 1, ncol = 2)
  IGR = matrix(NA, nrow = times + 1, ncol = 2)
  
  colnames(n) = spp
  colnames(p) = spp
  colnames(IGR) = spp
  
  ### initial conditions
  
  n[1, ] = n.init
  p[1, ] = p.init
  
  ### simulation loop
  
  for(t in 1:times){
    
    for (i in 1:2){
      
      A = alpha[i,1] * (g$SP * (1-p[t,1]) + p[t,1]) * n[t,1] + alpha[i,2] * (g$EG * (1-p[t,2]) + p[t,2]) * n[t,2]
      B = beta$s * sum(as.numeric(g) * (1-p[t,]) * n[t,]) + beta$m * sum(p[t,] * n[t,])
      
      seedlings.contribution.pop = (1-p[t,i]) * (1 - g[[i]] + g[[i]] * v[[i]] / (1 + B))
      adult.contribution.pop = p[t,i] * (omega[[i]] + lambda[[i]]/(1 + A))
      
      n[t+1,i] = n[t,i] * (seedlings.contribution.pop + adult.contribution.pop)
      
      p.num = omega[[i]] * p[t,i] + (1-p[t,i]) * g[[i]] * v[[i]] / (1+B)
      p.denom = p[t,i] * (omega[[i]] + lambda[[i]] / (1+A)) + (1-p[t,i]) * (1 - g[[i]] + g[[i]] * v[[i]] / (1+ B))
      p[t+1,i] = p.num / p.denom
      
      IGR[t+1,i] = log(seedlings.contribution.pop + adult.contribution.pop)
      
    }
    
  }
  
  return(list(SP.density = n[,1], SP.adultprop = p[,1], EG.density = n[,2], EG.adultprop = p[,2], IGR.SP = IGR[,1], IGR.EG = IGR[,2]))
}

#####

vars <- c("SP.density", "EG.density", "SP.adultprop", "EG.adultprop")
n.runs <- 2000 * 2 * 100

res <- matrix(
  NA_real_, nrow = n.runs, ncol = 8 + 2 * length(vars),
  dimnames = list(NULL, c("iter", "k", "j", "i",
                          "SP.n.init", "EG.n.init", "SP.p.init", "EG.p.init",
                          paste0(vars, ".mean"), paste0(vars, ".var")))
)
row <- 0L
# 
# for (k in 1:2000) {
#   t0 <- Sys.time()
#   for (j in 1:2) {
#     for (i in 1:100) {
#       n.init <- c(0, 0)
#       n.init[j] <- runif(1, 0.0001, 20000)
#       p.init <- runif(2, 0, 1)
#       output <- simulate.perennial.comp(k = k, n.init = n.init, p.init = p.init)
#       
#       last <- lapply(output[vars], tail, 500)
#       row <- row + 1L
#       res[row, ] <- c(row, k, j, i, n.init, p.init,
#                       vapply(last, mean, numeric(1)),
#                       vapply(last, var, numeric(1)))
#     }
#   }
#   cat(sprintf("k = %d: %.2f s\n", k, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
# }

############# parallel

library(parallel)
library(pbapply)

cols <- c("iter", "k", "j", "i", "SP.n.init", "EG.n.init", "SP.p.init", "EG.p.init",
          paste0(vars, ".mean"), paste0(vars, ".var"))

run.k <- function(k) {
  set.seed(1245 + k)
  out <- matrix(NA_real_, nrow = 200, ncol = length(cols), dimnames = list(NULL, cols))
  row <- 0L
  for (j in 1:2) {
    for (i in 1:100) {
      n.init <- c(0, 0)
      n.init[j] <- runif(1, 0.0001, 20000)
      p.init <- runif(2, 0, 1)
      output <- simulate.perennial.comp(k = k, n.init = n.init, p.init = p.init)
      last <- lapply(output[vars], tail, 500)
      row <- row + 1L
      out[row, ] <- c(NA, k, j, i, n.init, p.init,
                      vapply(last, mean, numeric(1)),
                      vapply(last, var, numeric(1)))
    }
  }
  out
}

n.cores <- detectCores(logical = FALSE) - 1
cl <- makeCluster(n.cores)

clusterExport(cl, c("simulate.perennial.comp", "run.k", "cols", "vars",
                    "g_SP", "g_EG", "comptogerm", "SP.lambda", "EG.lambda",
                    "alpha_SP", "alpha_EG", "omega_SP", "omega_EG"))

pboptions(type = "timer")
res <- pblapply(1:2000, run.k, cl = cl) |> do.call(what = rbind)
stopCluster(cl)

res[, "iter"] <- seq_len(nrow(res))
res <- as.data.frame(res)

tol <- 1e-6
means <- paste0(vars, ".mean")
vrs   <- paste0(vars, ".var")

## spread of the long-run means across the 100 i's, for each (k, j)
spread <- aggregate(res[means], res[c("k", "j")], \(x) diff(range(x)) / max(abs(x), 1))

## max within-run variance over the last 500 steps (non-zero means cycling, not an equilibrium)
osc <- aggregate(res[vrs], res[c("k", "j")], max)

chk <- merge(spread, osc, by = c("k", "j"))
chk$same    <- rowSums(chk[means] > tol) == 0
chk$settled <- rowSums(chk[vrs] > tol) == 0

table(same = chk$same, settled = chk$settled)
subset(chk, !same | !settled) |> head()
