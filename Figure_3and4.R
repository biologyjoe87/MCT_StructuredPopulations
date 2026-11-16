### Figure 3 and 4
### perennial competition

library(ggplot2)
library(dplyr)
library(tidyr)

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

######################## calculate mutual invasibility for all k

IGRs = matrix(NA, 2000, 2, dimnames = list(NULL, c("SP.IGR", "EG.IGR")))

for (k.temp in 1:2000) {
  
  # EG resident, SP invades
  EGpresent = simulate.perennial.comp(k = k.temp, n.init = c(0, 1), p.init = c(0.1, 0.1))
  burnin = (length(EGpresent$IGR.SP) * 1/4 + 1):length(EGpresent$IGR.SP)
  IGRs[k.temp, 1] = mean(EGpresent$IGR.SP[burnin])
  
  # SP resident, EG invades
  SPpresent = simulate.perennial.comp(k = k.temp, n.init = c(1, 0), p.init = c(0.1, 0.1))
  burnin = (length(SPpresent$IGR.EG) * 1/4 + 1):length(SPpresent$IGR.EG)
  IGRs[k.temp, 2] = mean(SPpresent$IGR.EG[burnin])
  print(k.temp)
}


library(ggplot2)

IGRs_df <- data.frame(IGRs)

IGRs_df$quadrant <- case_when(
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR > 0 ~ "Coexistence",
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR < 0 ~ "SP wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR > 0 ~ "EG wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR < 0 ~ "Priority effects"
)

tab1=prop.table(table(IGRs_df$quadrant))

mutinv1=ggplot(IGRs_df, aes(x = SP.IGR, y = EG.IGR, color = quadrant)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Coexistence"    = "green3",
                                "SP wins"        = "orange2",
                                "EG wins"        = "yellow4",
                                "Priority effects" = "cornflowerblue")) +
  labs(title = 'Both with perenniality',
       x = "SP invasion growth rate when EG at equilibrium",
       y = "EG invasion growth rate when SP at equilibrium",
       color = NULL) +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")

############# EXAMPLE OF EACH CASE

# COEXISTENCE

output=simulate.perennial.comp(k=303, n.init = c(10,10), p.init = c(0.1,0.1))

out_df <- data.frame(
  time = 0:(length(output$SP.density) - 1),
  
  SP_Adult = log(output$SP.adultprop * output$SP.density+1),
  EG_Adult = log(output$EG.adultprop * output$EG.density + 1),
  
  SP_Seed = log((1-output$SP.adultprop) * output$SP.density + 1),
  EG_Seed = log((1-output$EG.adultprop) * output$EG.density + 1)
)

plot_df <- out_df |>
  tidyr::pivot_longer(-time, names_to = "series", values_to = "log_density") |>
  tidyr::separate(series, into = c("species", "stage"), sep = "_")

ggplot(plot_df, aes(time, log_density, color = species, linetype = stage)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 16) +
  labs(x = "Time", y = "Log (Density + 1)")

# EG WINS

output=simulate.perennial.comp(k=331, n.init = c(10,10), p.init = c(0.1,0.1))

out_df <- data.frame(
  time = 0:(length(output$SP.density) - 1),
  
  SP_Adult = log(output$SP.adultprop * output$SP.density+1),
  EG_Adult = log(output$EG.adultprop * output$EG.density + 1),
  
  SP_Seed = log((1-output$SP.adultprop) * output$SP.density + 1),
  EG_Seed = log((1-output$EG.adultprop) * output$EG.density + 1)
)

plot_df <- out_df |>
  tidyr::pivot_longer(-time, names_to = "series", values_to = "log_density") |>
  tidyr::separate(series, into = c("species", "stage"), sep = "_")

ggplot(plot_df, aes(time, log_density, color = species, linetype = stage)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 14) +
  labs(x = "Time", y = "Log (Density + 1)")

# SP WINS

output=simulate.perennial.comp(k=305, n.init = c(10,10), p.init = c(0.1,0.1))

out_df <- data.frame(
  time = 0:(length(output$SP.density) - 1),
  
  SP_Adult = log(output$SP.adultprop * output$SP.density+1),
  EG_Adult = log(output$EG.adultprop * output$EG.density + 1),
  
  SP_Seed = log((1-output$SP.adultprop) * output$SP.density + 1),
  EG_Seed = log((1-output$EG.adultprop) * output$EG.density + 1)
)

plot_df <- out_df |>
  tidyr::pivot_longer(-time, names_to = "series", values_to = "log_density") |>
  tidyr::separate(series, into = c("species", "stage"), sep = "_")

ggplot(plot_df, aes(time, log_density, color = species, linetype = stage)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 14) +
  labs(x = "Time", y = "Log (Density + 1)")

# PRIORITY EFFECT

output=simulate.perennial.comp(k=307, n.init = c(10,10), p.init = c(0.5,0.1), times = 1500)

out_df <- data.frame(
  time = 0:(length(output$SP.density) - 1),
  
  SP_Adult = log(output$SP.adultprop * output$SP.density+1),
  EG_Adult = log(output$EG.adultprop * output$EG.density + 1),
  
  SP_Seed = log((1-output$SP.adultprop) * output$SP.density + 1),
  EG_Seed = log((1-output$EG.adultprop) * output$EG.density + 1)
)

plot_df <- out_df |>
  tidyr::pivot_longer(-time, names_to = "series", values_to = "log_density") |>
  tidyr::separate(series, into = c("species", "stage"), sep = "_")

p1 = ggplot(plot_df, aes(time, log_density, color = species, linetype = stage)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = "Time", y = "Log (Density + 1)")

output=simulate.perennial.comp(k=307, n.init = c(10,10), p.init = c(0.1,0.5), times = 1500)

out_df <- data.frame(
  time = 0:(length(output$SP.density) - 1),
  
  SP_Adult = log(output$SP.adultprop * output$SP.density+1),
  EG_Adult = log(output$EG.adultprop * output$EG.density + 1),
  
  SP_Seed = log((1-output$SP.adultprop) * output$SP.density + 1),
  EG_Seed = log((1-output$EG.adultprop) * output$EG.density + 1)
)

plot_df <- out_df |>
  tidyr::pivot_longer(-time, names_to = "series", values_to = "log_density") |>
  tidyr::separate(series, into = c("species", "stage"), sep = "_")

p2 = ggplot(plot_df, aes(time, log_density, color = species, linetype = stage)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_cartesian(ylim = c(NA, 10)) +
  labs(x = "Time", y = "")

###########################################################

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

## calculate mutual invasibility for all k

IGRs = matrix(NA, nrow = 2000, ncol = 2)
colnames(IGRs) = c('SP.IGR', 'EG.IGR') 

for (k.temp in 1:2000){
  
  EGpresent = simulate.perennial.comp.noEGperenniality(k=k.temp, n.init = c(0,1), p.init = c(0.1, 0.1))
  IGRs[k.temp,1] = mean(EGpresent$IGR.SP[(length(EGpresent$IGR.SP) * 1/4 + 1):length(EGpresent$IGR.SP)])
  
  SPpresent = simulate.perennial.comp.noEGperenniality(k=k.temp, n.init = c(1,0), p.init = c(0.1, 0.1))
  IGRs[k.temp,2] = mean(SPpresent$IGR.EG[(length(SPpresent$IGR.EG) * 1/4 + 1):length(SPpresent$IGR.EG)])
  
  print(k.temp)
}

IGRs_df <- data.frame(IGRs)

IGRs_df$quadrant <- case_when(
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR > 0 ~ "Coexistence",
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR < 0 ~ "SP wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR > 0 ~ "EG wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR < 0 ~ "Priority effects"
)

tab2=prop.table(table(IGRs_df$quadrant))

mutinv2=ggplot(IGRs_df, aes(x = SP.IGR, y = EG.IGR, color = quadrant)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Coexistence"    = "green3",
                                "SP wins"        = "orange2",
                                "EG wins"        = "yellow4",
                                "Priority effects" = "cornflowerblue")) +
  labs(title = 'No perenniality for EG',
       x = "SP invasion growth rate when EG at equilibrium",
       y = "EG invasion growth rate when SP at equilibrium",
       color = NULL) +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")

######################  NO SP PERENNIALITY

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

## calculate mutual invasibility for all k

IGRs = matrix(NA, nrow = 2000, ncol = 2)
colnames(IGRs) = c('SP.IGR', 'EG.IGR') 

for (k.temp in 1:2000){
  
  EGpresent = simulate.perennial.comp.noSPperenniality(k=k.temp, n.init = c(0,1), p.init = c(0.1, 0.1))
  IGRs[k.temp,1] = mean(EGpresent$IGR.SP[(length(EGpresent$IGR.SP) * 1/4 + 1):length(EGpresent$IGR.SP)])
  
  SPpresent = simulate.perennial.comp.noSPperenniality(k=k.temp, n.init = c(1,0), p.init = c(0.1, 0.1))
  IGRs[k.temp,2] = mean(SPpresent$IGR.EG[(length(SPpresent$IGR.EG) * 1/4 + 1):length(SPpresent$IGR.EG)])
  
  print(k.temp)
}

IGRs_df <- data.frame(IGRs)

IGRs_df$quadrant <- case_when(
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR > 0 ~ "Coexistence",
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR < 0 ~ "SP wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR > 0 ~ "EG wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR < 0 ~ "Priority effects"
)

tab3=prop.table(table(IGRs_df$quadrant))

mutinv3=ggplot(IGRs_df, aes(x = SP.IGR, y = EG.IGR, color = quadrant)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Coexistence"    = "green3",
                                "SP wins"        = "orange2",
                                "EG wins"        = "yellow4",
                                "Priority effects" = "cornflowerblue")) +
  labs(title = 'No perenniality for SP',
       x = "SP invasion growth rate when EG at equilibrium",
       y = "EG invasion growth rate when SP at equilibrium",
       color = NULL) +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")

###################### NEITHER WITH PERENNIALITY

simulate.comp.annual <- function(times = 2000,
                                 k = 1,
                                 n.init = c(0,0)){
  
  #### species names
  spp <- c("SP", "EG")
  
  #### demographic parameters
  g <- c(
    SP = g_SP,
    EG = g_EG
  )
  
  lambda <- c(
    SP = SP.lambda[k],
    EG = EG.lambda[k]
  )
  
  #### alpha competition matrix
  alpha <- matrix(
    c(alpha_SP[k,"alpha_{SP,SP}"],
      alpha_SP[k,"alpha_{SP,EG}"],
      alpha_EG[k,"alpha_{EG,SP}"],
      alpha_EG[k,"alpha_{EG,EG}"]),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(alpha) <- colnames(alpha) <- spp
  
  #### storage
  
  n <- matrix(NA, nrow = times + 1, ncol = 2)
  IGR <- matrix(NA, nrow = times + 1, ncol = 2)
  
  colnames(n) <- spp
  colnames(IGR) <- spp
  
  #### initial conditions
  
  n[1,] <- n.init
  
  #### simulation
  
  for(t in 1:times){
    
    #### competition terms
    
    A_SP <-
      alpha["SP","SP"] * g["SP"] * n[t,1] +
      alpha["SP","EG"] * g["EG"] * n[t,2]
    
    A_EG <-
      alpha["EG","SP"] * g["SP"] * n[t,1] +
      alpha["EG","EG"] * g["EG"] * n[t,2]
    
    #### SP
    
    growth_SP <-
      (1 - g["SP"]) +
      (lambda["SP"] * g["SP"]) / (1 + A_SP)
    
    n[t+1,1] <- n[t,1] * growth_SP
    
    IGR[t+1,1] <- log(growth_SP)
    
    #### EG
    
    growth_EG <-
      (1 - g["EG"]) +
      (lambda["EG"] * g["EG"]) / (1 + A_EG)
    
    n[t+1,2] <- n[t,2] * growth_EG
    
    IGR[t+1,2] <- log(growth_EG)
  }
  
  list(
    SP.density = n[,1],
    EG.density = n[,2],
    IGR.SP = IGR[,1],
    IGR.EG = IGR[,2]
  )
}

IGRs <- matrix(NA, nrow = 2000, ncol = 2)
colnames(IGRs) <- c("SP.IGR","EG.IGR")

for(k.temp in 1:2000){
  
  #### SP invading EG resident
  
  EGpresent <- simulate.comp.annual(
    k = k.temp,
    n.init = c(0,1)
  )
  
  IGRs[k.temp,1] <-
    mean(
      tail(EGpresent$IGR.SP,
           floor(length(EGpresent$IGR.SP)*3/4))
    )
  
  #### EG invading SP resident
  
  SPpresent <- simulate.comp.annual(
    k = k.temp,
    n.init = c(1,0)
  )
  
  IGRs[k.temp,2] <-
    mean(
      tail(SPpresent$IGR.EG,
           floor(length(SPpresent$IGR.EG)*3/4))
    )
  
  print(k.temp)
}

IGRs_df <- data.frame(IGRs)

IGRs_df$quadrant <- case_when(
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR > 0 ~ "Coexistence",
  IGRs_df$SP.IGR > 0 & IGRs_df$EG.IGR < 0 ~ "SP wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR > 0 ~ "EG wins",
  IGRs_df$SP.IGR < 0 & IGRs_df$EG.IGR < 0 ~ "Priority effects"
)

tab4=prop.table(table(IGRs_df$quadrant))

mutinv4=ggplot(IGRs_df, aes(x = SP.IGR, y = EG.IGR, color = quadrant)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Coexistence"    = "green3",
                                "SP wins"        = "orange2",
                                "EG wins"        = "yellow4",
                                "Priority effects" = "cornflowerblue")) +
  labs(title = 'No perenniality for either',
       x = "SP invasion growth rate when EG at equilibrium",
       y = "EG invasion growth rate when SP at equilibrium",
       color = NULL) +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")

########### PLOT

library(patchwork)

(mutinv1 + mutinv2) / (mutinv3+mutinv4)
