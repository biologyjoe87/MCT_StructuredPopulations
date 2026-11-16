### Figure 5
### Author: Joe Brennan

library(dplyr)
library(ggplot2)
library(deSolve)

parms <- c(
  r = 3.0,        
  a_IP = 1.0,
  a_FP = 4.0,     
  h = 2.0,       
  tau_FI = 6.0,   
  tau_IF = 1.0,
  c = 0.1,
  mu = 0.01,
  invading.r = 2.64,
  invading.a_IP = 2, 
  invading.a_FP = 5,
  invading.tau_FI = 3,
  invading.tau_IF = 2
)

# model
foraging_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    dxdt <- x * ( r*y*(1 - x) - (P * (a_IP*(1 - y) + a_FP*y)) / (1 + h * x * (a_IP*(1 - y) + a_FP*y)) )
    
    dydt <- r * y * (1 - y) * (1 - x) + (1 / tau_FI) * (1 - y) - (1 / tau_IF) * y -
      (P * y * (1 - y) * (a_FP - a_IP)) / (1 + h * x * (a_IP*(1 - y) + a_FP*y))
    
    dPdt <- P * ( c * x * (a_IP*(1 - y) + a_FP*y) / (1 + h * x * (a_IP*(1 - y) + a_FP*y)) - mu )
    
    dy2dt = invading.r * y2 * (1 - y2) + (1 / invading.tau_FI) * (1 - y2) - (1 / invading.tau_IF) * y2 -
      (P * y2 * (1 - y2) * (invading.a_FP - invading.a_IP)) / (1 + h * x * (a_IP*(1 - y) + a_FP*y)) #assume x2=0
    
    list(c(dxdt, dydt, dPdt, dy2dt))
  })
}

# initial conditions
state <- c(x = 0.5, y = 0.5, P = 0.1, y2 = 0.1)
times <- seq(0, 1000, by = 0.1)

# simulate
out <- ode(y = state, times = times, func = foraging_model, parms = parms)
out_df <- as.data.frame(out)
plot_out_df = out_df[7000:10000,]
plot_out_df$time = plot_out_df$time-699
plot_out_df$x2 = rep(0,length(plot_out_df$x))

# reshape for ggplot
out_long <- pivot_longer(plot_out_df, cols = c(x, y, P, x2, y2), names_to = "species", values_to = "abundance")
out_long$behavior <- ifelse(out_long$species %in% c("y", "y2"),
                            "Proportion foraging",
                            "Density")
out_long$speciesidentity <- factor(
  case_when(
    out_long$species == "x"  ~ "Resident Prey",
    out_long$species == "x2" ~ "Invading Prey",
    out_long$species == "P"  ~ "Predator",
    out_long$species == "y"  ~ "Resident Prey",
    out_long$species == "y2" ~ "Invading Prey",
    TRUE ~ out_long$species
  ),
  levels = c("Resident Prey", "Invading Prey", "Predator")
)

out_long$abundance_rescaled <- ifelse(out_long$species %in% c("y", "y2"),
                                      out_long$abundance,
                                      out_long$abundance^(1/5))
out_long$behavior <- ifelse(out_long$species %in% c("y", "y2"),
                            "Proportion~foraging",
                            "Density^{1/5}")

p1 <- ggplot(out_long,
             aes(time, abundance_rescaled, color = speciesidentity)) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 18) +
  facet_wrap(~ behavior, scales = "free_y", labeller = label_parsed) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time", y = "Value", color = "Species")

################ SPECIES 2 CANNOT INVADE

twospp_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    denom = (1 + h * (x1 * (a_IP*(1 - y1) + a_FP*y1) + x2 * (invading.a_IP * (1-y2) + invading.a_FP * y2)))
    
    dx1dt <- x1 * ( r*y1*(1 - x1) - (P * (a_IP*(1 - y1) + a_FP*y1)) / denom )
    
    dy1dt <- r * y1 * (1 - y1) * (1 - x1) + (1 / tau_FI) * (1 - y1) - (1 / tau_IF) * y1 -
      (P * y1 * (1 - y1) * (a_FP - a_IP)) / denom
    
    dPdt <- P * ( c * (x1 * (a_IP*(1 - y1) + a_FP*y1) + x2 * (invading.a_IP*(1 - y2) + invading.a_FP*y2)) / denom - mu )
    
    dx2dt <- x2 * ( invading.r*y2*(1 - x2) - (P * (invading.a_IP*(1 - y2) + invading.a_FP*y2)) / denom )
    
    dy2dt = invading.r * y2 * (1 - y2) + (1 / invading.tau_FI) * (1 - y2) - (1 / invading.tau_IF) * y2 -
      (P * y2 * (1 - y2) * (invading.a_FP - invading.a_IP)) / denom 
    
    list(c(dx1dt, dy1dt, dPdt, dx2dt, dy2dt))
  })
}

# parameters
parms <- c(
  r = 3.0,         # faster prey reproduction -> stronger rebounds
  a_IP = 1.0,
  a_FP = 4.0,      # much higher attack on fearful individuals 4 worked
  h = 2.0,         # strong handling/saturation -> nonlinear functional response
  tau_FI = 6.0,    # slow F -> I transition (adds lag)
  tau_IF = 1.0,
  f = 0.0,         # strong suppression of transition by P (nonlinear) 8 worked
  c = 0.1,         
  mu = 0.01,       # low predator mortality => predators persist longer. was 0.01
  invading.r = 2.64, # was 2.6
  invading.a_IP = 2, # was 1
  invading.a_FP = 5, # was 4
  invading.tau_FI = 3,
  invading.tau_IF = 2,
  invading.f = 0
)

# initial conditions
state <- c(x1 = 0.5, y1 = 0.5, P = 0.1, x2 = 0.0, y2 = 0.0)

# simulate until time 1000 (no x2)
times1 <- seq(0, 1000, by = 0.1)
out1 <- ode(y = state, times = times1, func = twospp_model, parms = parms)

# take last state and add x2 at low density
state2 <- as.list(out1[nrow(out1), -1])  # drop time column
state2$x2 <- 1e-4  # introduce x2 at low density

# simulate from 1000 to, say, 2000
times2 <- seq(1000, 2000, by = 0.1)
out2 <- ode(y = unlist(state2), times = times2, func = twospp_model, parms = parms)

library(tidyr)

out_df <- rbind(out1, out2[-1, ]) |> as.data.frame()

out_long <- out_df |>
  pivot_longer(c(x1, y1, P, x2, y2),
               names_to = "species", values_to = "abundance") |>
  transform(
    is_y = species %in% c("y1", "y2"),
    speciesidentity = factor(
      c(x1 = "Resident Prey", x2 = "Invading Prey", P = "Predator",
        y1 = "Resident Prey", y2 = "Invading Prey")[species],
      levels = c("Resident Prey", "Invading Prey", "Predator")
    )
  ) |>
  transform(
    behavior = ifelse(is_y, "Proportion~foraging", "Density^{1/5}"),
    abundance_rescaled = ifelse(is_y, abundance, abundance^(1/5))
  )

blank <- data.frame(
  behavior = "Proportion~foraging",
  time = min(out_long$time),
  abundance_rescaled = c(0, 1),
  speciesidentity = "Resident Prey"
)

p2 <- ggplot(out_long, aes(time, abundance_rescaled, color = speciesidentity)) +
  geom_line(linewidth = 1) +
  geom_blank(data = blank) +
  theme_bw(base_size = 18) +
  facet_wrap(~ behavior, scales = "free_y", labeller = label_parsed) +
  labs(x = "Time", y = "Value", color = "Species")

########## DECOMPOSITION OF GROWTH RATE OF INVADING SPP

# initial conditions
state <- c(x = 0.5, y = 0.5, P = 0.1, y2 = 0.1)
times <- seq(0, 20000, by = 0.1)

# simulate
out <- ode(y = state, times = times, func = foraging_model, parms = parms)
out_df <- as.data.frame(out)

# evaluate impact of fluctuating y and P
IGR.resident = mean(parms['r']*out_df$y*(1 - out_df$x) - (out_df$P * (parms['a_IP']*(1 - out_df$y) + parms['a_FP']*out_df$y)) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - out_df$y) + parms['a_FP']*out_df$y)))
IGR.invader = mean(parms['invading.r']*out_df$y2 - (out_df$P * (parms['invading.a_IP']*(1 - out_df$y2) + parms['invading.a_FP']*out_df$y2)) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - out_df$y) + parms['a_FP']*out_df$y)))

mean.y = mean(out[,'y'])
mean.P = mean(out[,'P'])

baseline = mean(parms['invading.r']*out[,'y2'] - (mean.P * (parms['invading.a_IP']*(1 - out[,'y2']) + parms['invading.a_FP']*out[,'y2'])) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - mean.y) + parms['a_FP']*mean.y))) -
  mean(parms['r']*mean.y*(1 - out_df$x) - (mean.P * (parms['a_IP']*(1 - mean.y) + parms['a_FP']*mean.y)) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - mean.y) + parms['a_FP']*mean.y)))
fluc.P = mean(parms['invading.r']*out[,'y2'] - (out[,'P'] * (parms['invading.a_IP']*(1 - out[,'y2']) + parms['invading.a_FP']*out[,'y2'])) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - mean.y) + parms['a_FP']*mean.y))) -
  mean(parms['r']*mean.y*(1 - out_df$x) - (out[,'P'] * (parms['a_IP']*(1 - mean.y) + parms['a_FP']*mean.y)) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - mean.y) + parms['a_FP']*mean.y))) -
  baseline
fluc.y = mean(parms['invading.r']*out[,'y2'] - (mean.P * (parms['invading.a_IP']*(1 - out[,'y2']) + parms['invading.a_FP']*out[,'y2'])) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - out[,'y']) + parms['a_FP']*out[,'y']))) -
  mean(parms['r']*out[,'y']*(1 - out_df$x) - (mean.P * (parms['a_IP']*(1 - out[,'y']) + parms['a_FP']*out[,'y'])) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - out[,'y']) + parms['a_FP']*out[,'y']))) -
  baseline

shuffled.y = sample(out_df$y)

variance.int = IGR.invader-fluc.y-fluc.P-baseline
nocov.varianceint = (
  mean(parms['invading.r']*out[,'y2'] - (out[,'P'] * (parms['invading.a_IP']*(1 - out[,'y2']) + parms['invading.a_FP']*out[,'y2'])) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - shuffled.y) + parms['a_FP']*shuffled.y))) -
    mean(parms['r']*shuffled.y*(1 - out_df$x) - (out[,'P'] * (parms['a_IP']*(1 - shuffled.y) + parms['a_FP']*shuffled.y)) / (1 + parms['h'] * out_df$x * (parms['a_IP']*(1 - shuffled.y) + parms['a_FP']*shuffled.y))) -
    baseline-fluc.P-fluc.y
)

storage = variance.int - nocov.varianceint

bars <- data.frame(
  Component = c('Total',"Baseline", "Fluctuations in\npredator density", "Fluctuations in\nproportion of\nresident prey foraging", 'Variance interaction', "Storage effect"),
  Value = c(baseline + fluc.P + fluc.y + storage + nocov.varianceint,baseline, fluc.P, fluc.y, nocov.varianceint,storage)
)
bars$Component <- factor(bars$Component, levels = rev(unique(bars$Component)))

p3=ggplot(bars, aes(x = Component, y = Value, fill = Component)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none") +
  ylab("Contribution to IGR of Invading Prey") +
  xlab("") +
  scale_fill_manual(values=rep("grey",6))+
  ggtitle("Decomposition of Invasion\nGrowth Rate")+
  geom_vline(xintercept=5.5, linetype = "dotted", color = "black", linewidth =2) +
  geom_hline(yintercept=0) +
  coord_flip()
