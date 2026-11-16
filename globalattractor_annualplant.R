### global attractor for single species when perenniality removed?

library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(303)

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

simulate.perennial.comp.noEGperenniality = function(times = 5000, k = 1, n.init = c(0,0), p.init = c(0,0)){
  
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
  # seedlings -> beta_S
  # adults -> beta_M
  
  beta = list(
    s = comptogerm$beta_S[k],
    m = comptogerm$beta_M[k],
    a = comptogerm$beta_A[k]
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
    
    #### Competition terms
    
    A1 <- alpha["SP","SP"] *
      (g$SP * (1 - p[t,1]) + p[t,1]) * n[t,1] +
      alpha["SP","EG"] *
      g$EG * n[t,2]
    
    B <- beta$s * (g$SP * (1 - p[t,1]) * n[t,1]) +
      beta$m * (p[t,1] * n[t,1]) + 
      beta$a * g$EG * n[t,2]
    
    #### Species 1 (perennial)
    
    seedlings.contribution.pop <-
      (1 - p[t,1]) *
      (1 - g$SP + g$SP * v$SP / (1 + B))
    
    adult.contribution.pop <-
      p[t,1] *
      (omega$SP + lambda$SP / (1 + A1))
    
    growth1 <- seedlings.contribution.pop +
      adult.contribution.pop
    
    n[t+1,1] <- n[t,1] * growth1
    
    p.num <-
      omega$SP * p[t,1] +
      (1 - p[t,1]) * g$SP * v$SP / (1 + B)
    
    p.denom <-
      p[t,1] *
      (omega$SP + lambda$SP / (1 + A1)) +
      (1 - p[t,1]) *
      (1 - g$SP + g$SP * v$SP / (1 + B))
    
    p[t+1,1] <- p.num / p.denom
    
    IGR[t+1,1] <- log(growth1)
    
    #### Species 2
    
    A2 <- alpha["EG","SP"] *
      (g$SP * (1 - p[t,1]) + p[t,1]) * n[t,1] +
      alpha["EG","EG"] *
      g$EG * n[t,2]
    
    growth2 <-
      (1 - g$EG) +
      (lambda$EG * g$EG) / (1 + A2)
    
    n[t+1,2] <- n[t,2] * growth2
    
    IGR[t+1,2] <- log(growth2)
    
    ## species 2 has no p-dynamics
    p[t+1,2] <- 0
  }
  
  return(list(SP.density = n[,1], SP.adultprop = p[,1], EG.density = n[,2], EG.adultprop = p[,2], IGR.SP = IGR[,1], IGR.EG = IGR[,2]))
}

simulate.perennial.comp.noSPperenniality <- function(times = 2000,
                                                     k = 1,
                                                     n.init = c(0,0),
                                                     p.init = c(0,0)){
  
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
      c(alpha_SPSP = alpha_SP[k, "alpha_{SP,SP}"],
        alpha_SPEG = alpha_SP[k, "alpha_{SP,EG}"],
        alpha_EGSP = alpha_EG[k, "alpha_{EG,SP}"],
        alpha_EGEG = alpha_EG[k, "alpha_{EG,EG}"]),
      nrow = 2,
      byrow=T
    )
  
  colnames(alpha) <- rownames(alpha) <- c("SP", "EG")
  
  #### beta parameters
  beta = list(
    s = comptogerm$beta_S[k],
    m = comptogerm$beta_M[k],
    a = comptogerm$beta_A[k]
  )
  
  #### storage vectors
  
  n = matrix(NA, nrow = times + 1, ncol = 2)
  p = matrix(NA, nrow = times + 1, ncol = 2)
  IGR = matrix(NA, nrow = times + 1, ncol = 2)
  
  colnames(n) = spp
  colnames(p) = spp
  colnames(IGR) = spp
  
  #### initial conditions
  
  n[1, ] = n.init
  p[1, ] = p.init
  
  #### simulation loop
  
  for(t in 1:times){
    
    #### Species 1 (SP) -- annual / no perenniality

    A1 <- alpha["SP","SP"] *
      g$SP * n[t,1] +
      alpha["SP","EG"] *
      (g$EG * (1 - p[t,2]) + p[t,2]) * n[t,2]
    
    growth1 <-
      (1 - g$SP) +
      (lambda$SP * g$SP) / (1 + A1)
    
    n[t+1,1] <- n[t,1] * growth1
    
    IGR[t+1,1] <- log(growth1)
    
    ## no stage structure
    p[t+1,1] <- 0
    
    
    #### Species 2 (EG) -- perennial

    A2 <- alpha["EG","SP"] *
      g$SP * n[t,1] +
      alpha["EG","EG"] *
      (g$EG * (1 - p[t,2]) + p[t,2]) * n[t,2]
    
    B <- beta$a *
      g$SP * n[t,1] +
      beta$s *
      (g$EG * (1 - p[t,2]) * n[t,2]) +
      beta$m *
      (p[t,2] * n[t,2])
    
    seedlings.contribution.pop <-
      (1 - p[t,2]) *
      (1 - g$EG + g$EG * v$EG / (1 + B))
    
    adult.contribution.pop <-
      p[t,2] *
      (omega$EG + lambda$EG / (1 + A2))
    
    growth2 <-
      seedlings.contribution.pop +
      adult.contribution.pop
    
    n[t+1,2] <- n[t,2] * growth2
    
    p.num <-
      omega$EG * p[t,2] +
      (1 - p[t,2]) *
      g$EG * v$EG / (1 + B)
    
    p.denom <-
      p[t,2] *
      (omega$EG + lambda$EG / (1 + A2)) +
      (1 - p[t,2]) *
      (1 - g$EG + g$EG * v$EG / (1 + B))
    
    p[t+1,2] <- p.num / p.denom
    
    IGR[t+1,2] <- log(growth2)
  }
  
  return(
    list(
      SP.density = n[,1],
      SP.adultprop = p[,1],
      EG.density = n[,2],
      EG.adultprop = p[,2],
      IGR.SP = IGR[,1],
      IGR.EG = IGR[,2]
    )
  )
}

vars <- c("SP.density", "EG.density", "SP.adultprop", "EG.adultprop")
n.runs <- 2000 * 2 * 100

res <- matrix(
  NA_real_, nrow = n.runs, ncol = 8 + 2 * length(vars),
  dimnames = list(NULL, c("iter", "k", "j", "i",
                          "SP.n.init", "EG.n.init", "SP.p.init", "EG.p.init",
                          paste0(vars, ".mean"), paste0(vars, ".var")))
)
row <- 0L

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
      if (j==1) output <- simulate.perennial.comp.noSPperenniality(times = 3000, k = k, n.init = n.init, p.init = p.init)
      if (j==2) output <- simulate.perennial.comp.noEGperenniality(times = 3000, k = k, n.init = n.init, p.init = p.init)
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

clusterExport(cl, c("simulate.perennial.comp.noSPperenniality", "simulate.perennial.comp.noEGperenniality", "run.k", "cols", "vars",
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

############ CHECK THOSE THAT DID NOT CONVERGE. Do converge when given longer, therefore are globally stable 

## (k, j) pairs to re-check
targets <- data.frame(k = c(1456, 1552, 1584, 1745, 1871), j = 2)
ks <- unique(targets$k)

run.k <- function(k) {
  set.seed(1245 + k)
  js <- targets$j[targets$k == k]
  out <- matrix(NA_real_, nrow = 100 * length(js), ncol = length(cols),
                dimnames = list(NULL, cols))
  row <- 0L
  for (j in 1:2) {
    for (i in 1:100) {
      n.init <- c(0, 0)
      n.init[j] <- runif(1, 0.0001, 20000)
      p.init <- runif(2, 0, 1)
      if (!j %in% js) next  # keep RNG stream identical to the full run
      if (j == 1) output <- simulate.perennial.comp.noSPperenniality(times = 10000, k = k, n.init = n.init, p.init = p.init)
      if (j == 2) output <- simulate.perennial.comp.noEGperenniality(times = 10000, k = k, n.init = n.init, p.init = p.init)
      last <- lapply(output[vars], tail, 500)
      row <- row + 1L
      out[row, ] <- c(NA, k, j, i, n.init, p.init,
                      vapply(last, mean, numeric(1)),
                      vapply(last, var, numeric(1)))
    }
  }
  out
}

n.cores <- min(length(ks), detectCores(logical = FALSE) - 1)
cl <- makeCluster(n.cores)

clusterExport(cl, c("simulate.perennial.comp.noSPperenniality", "simulate.perennial.comp.noEGperenniality",
                    "run.k", "cols", "vars", "targets",
                    "g_SP", "g_EG", "comptogerm", "SP.lambda", "EG.lambda",
                    "alpha_SP", "alpha_EG", "omega_SP", "omega_EG"))

pboptions(type = "timer")
res <- pblapply(ks, run.k, cl = cl) |> do.call(what = rbind)
stopCluster(cl)

res[, "iter"] <- seq_len(nrow(res))
res <- as.data.frame(res)
