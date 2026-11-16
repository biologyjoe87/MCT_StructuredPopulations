### Figure 2
### Wang et al model from Park et al model Science
### Author: Joe Brennan

library(deSolve)
library(ggplot2)
library(reshape2)

FD_NO_func = function(parms){

  with(as.list(parms), {

    FD <- ((mu + delta1) * (mu + gamma2) * (beta1 - mu - gamma1)) /
      ((mu + delta2) * (mu + gamma1) * (beta2 - mu - gamma2)) *
      sqrt(
        (sigma21 * (mu + gamma1) * (delta2 + sigma22 * (mu + gamma2))) /
          (sigma12 * (mu + gamma2) * (delta1 + sigma11 * (mu + gamma1)))
      )

    NO <- sqrt(
      (sigma12 * sigma21 * (mu + gamma1) * (mu + gamma2)) /
        (
          (delta1 + sigma11 * (mu + gamma1)) *
            (delta2 + sigma22 * (mu + gamma2))
        )
    )

    return(list(FD = FD, NO = NO))

  })

}


params_H1N1_H3N2 <- c(
  mu     = 1/(52*80),
  delta1 = 1/(3.12*52),
  delta2 = 1/(3.08*52),
  beta1  = 1.44 * (1/(52*80) + 7/(2.64)),
  beta2  = 1.60 * (1/(52*80) + 7/3.03),
  gamma1 = 7/2.64,
  gamma2 = 7/3.03,
  sigma11 = 1,
  sigma12 = 0.17,
  sigma21 = 0.4,
  sigma22 = 1
)

sigma12.seq = seq(0.15,1,by=0.01)
FD.array = c()
NO.array = c()

for(i in sigma12.seq){

  params_H1N1_H3N2['sigma12'] = i
  FDNO = FD_NO_func(params_H1N1_H3N2)
  FD.array = c(FD.array, FDNO$FD)
  NO.array = c(NO.array, FDNO$NO)

}

niche_overlap <- seq(from = -1, to = 11, by = 0.01)
fitness_ratio_min <- niche_overlap
fitness_ratio_max <- 1 / niche_overlap
boundary_df <- data.frame(
  niche_overlap = niche_overlap,
  min_fitness_ratio = fitness_ratio_min,
  max_fitness_ratio = fitness_ratio_max
)

df <- data.frame(
  ND = NO.array,
  FD = FD.array,
  sigma12 = sigma12.seq
)

p1=ggplot() +
  theme_classic(base_size=18) +
  geom_ribbon(data = boundary_df,
              aes(x = niche_overlap, ymin = min_fitness_ratio, ymax = max_fitness_ratio),
              fill = 'grey80') +
  geom_line(data = boundary_df, aes(x = niche_overlap, y = max_fitness_ratio)) +
  geom_line(data = boundary_df, aes(x = niche_overlap, y = min_fitness_ratio)) +
  geom_path(data = df, aes(x = ND, y = FD, color = sigma12), size = 1.2) +
  geom_point(data = df, aes(x = ND, y = FD, color = sigma12), size = 2) +
  scale_color_viridis_c(name = expression(sigma[12])) +
  annotate("text", x = 1.5, y = 1, label = "Priority \nEffect", size = 5, fontface = "bold") +
  annotate("text", x = 0.6, y = 1.1, label = "Coexistence", size = 5, fontface = "bold") +
  annotate("text", x = 1, y = 0.4, label = "H3N2 wins", size = 5, fontface = "bold") +
  annotate("text", x = 1.15, y = 1.8, label = "H1N1 wins", size = 5, fontface = "bold") +
  labs(
    x = expression("Niche Overlap (" * rho * ")"),
    y = expression("Fitness Differences (" * italic(F) * ")")
  ) +
  coord_cartesian(xlim=c(0.1,2),ylim = c(0.1, 2))

######### SIMULATION COEXISTENCE

two_strain_model <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    dS1 <- mu * (1 - S1) + delta1 * (1 - S1 - I1) -
      sigma11 * beta1 * I1 * S1 -
      sigma12 * beta2 * I2 * S1
    
    dS2 <- mu * (1 - S2) + delta2 * (1 - S2 - I2) -
      sigma21 * beta1 * I1 * S2 -
      sigma22 * beta2 * I2 * S2
    
    dI1 <- beta1 * S1 * I1 - (gamma1 + mu) * I1
    
    dI2 <- beta2 * S2 * I2 - (gamma2 + mu) * I2
    
    list(c(dS1, dS2, dI1, dI2))
  })
}

params_H1N1_H3N2 <- c(
  mu     = 1/(52*80),
  delta1 = 1/(3.12*52),
  delta2 = 1/(3.08*52),
  beta1  = 1.44 * (1/(52*80) + 7/(2.64)),
  beta2  = 1.60 * (1/(52*80) + 7/3.03),
  gamma1 = 7/2.64,   
  gamma2 = 7/3.03,
  sigma11 = 1.0,
  sigma12 = 0.17,
  sigma21 = 0.4, 
  sigma22 = 1.0 
)

init <- c(
  S1 = 0.999,
  S2 = 0.999,
  I1 = 0.001,
  I2 = 0.001
)

times <- seq(0, 5000, by = 0.1)

out <- ode(y = init, times = times, func = two_strain_model, parms = params_H1N1_H3N2)
out_df <- as.data.frame(out)

plot_df <- out_df |>
  pivot_longer(-time, names_to = "pop", values_to = "density") |>
  mutate(
    number = gsub("[A-Z]", "", pop),
    letter = gsub("[0-9]", "", pop)
  ) |>
  filter(letter == "I")

p2=ggplot(plot_df,
          aes(time, density,
              color = number)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 18) +
  coord_cartesian(ylim = c(0, 0.005),xlim=c(0,1500)) +
  scale_color_discrete(
    labels = c("H1N1", "H3N2")
  ) +
  labs(
    x = "Time",
    y = "Proportion Infected",
    color = "Strain",
    linetype = "Class"
  )

######### SIMULATION H3N2 WIN

two_strain_model <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    dS1 <- mu * (1 - S1) + delta1 * (1 - S1 - I1) -
      sigma11 * beta1 * I1 * S1 -
      sigma12 * beta2 * I2 * S1
    
    dS2 <- mu * (1 - S2) + delta2 * (1 - S2 - I2) -
      sigma21 * beta1 * I1 * S2 -
      sigma22 * beta2 * I2 * S2
    
    dI1 <- beta1 * S1 * I1 - (gamma1 + mu) * I1
    
    dI2 <- beta2 * S2 * I2 - (gamma2 + mu) * I2
    
    list(c(dS1, dS2, dI1, dI2))
  })
}

params_H1N1_H3N2 <- c(
  mu     = 1/(52*80), 
  delta1 = 1/(3.12*52),
  delta2 = 1/(3.08*52),
  beta1  = 1.44 * (1/(52*80) + 7/(2.64)),
  beta2  = 1.60 * (1/(52*80) + 7/3.03), 
  gamma1 = 7/2.64,
  gamma2 = 7/3.03,
  sigma11 = 1.0,
  sigma12 = 0.75,
  sigma21 = 0.6,
  sigma22 = 1.0 
)

init <- c(
  S1 = 0.995,
  S2 = 0.999,
  I1 = 0.005,
  I2 = 0.001
)

times <- seq(0, 5000, by = 0.1)

out <- ode(y = init, times = times, func = two_strain_model, parms = params_H1N1_H3N2)
out_df <- as.data.frame(out)

plot_df <- out_df |>
  pivot_longer(-time, names_to = "pop", values_to = "density") |>
  mutate(
    number = gsub("[A-Z]", "", pop),
    letter = gsub("[0-9]", "", pop)
  ) |>
  filter(letter == "I")

p3=ggplot(plot_df,
          aes(time, density,
              color = number)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 18) +
  coord_cartesian(ylim = c(0, 0.005),xlim=c(0,2000)) +
  scale_color_discrete(
    labels = c("H1N1", "H3N2")
  ) +
  labs(
    x = "Time",
    y = "Proportion Infected",
    color = "Strain",
    linetype = "Class"
  )

library(patchwork)

p1 + (p2/p3)
