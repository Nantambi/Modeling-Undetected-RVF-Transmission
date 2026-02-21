###############################################################################
# SEIR MODEL FOR RIFT VALLEY FEVER (RVF) TRANSMISSION AMONG ABATTOIR STAFF
# NAMISINDWA DISTRICT, UGANDA (POPULATION = 283)
#
# Detection occurs ONLY from the infectious compartment.
# Detection does NOT affect transmission, so epidemic dynamics are identical
# under both surveillance scenarios.
#
# Five plots are generated:
#   1. SEIR dynamics (all compartments, Status Quo) – original.
#   2. Infectious individuals over time (both scenarios, overlapping).
#   3. Monthly detected cases (bar plot, both scenarios).
#   4. Sensitivity analysis: effect of varying spillover coefficient (klh).
#   5. Panel: S (top), R (middle), and E+I together (bottom) – Status Quo.
###############################################################################

library(deSolve)
library(rootSolve)
library(tidyverse)
library(patchwork)   # for combining plots

# NDVI DATA FOR NAMISINDWA 
# I used my cattle RVF model to load the NDVI dataset; you can download it as well; `final_ndvi_dataset` is already loaded.
ndvi_nam <- final_ndvi_dataset %>%
  filter(site_name == "Namisindwa")

if (nrow(ndvi_nam) == 0) {
  stop("No NDVI data found for Namisindwa")
}

# Monthly mean NDVI
monthly_ndvi <- ndvi_nam %>%
  group_by(month) %>%
  summarise(ndvi = mean(mean_ndvi, na.rm = TRUE)) %>%
  arrange(month)

# Approximate day-of-year for mid-month

month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
cum_days <- c(0, cumsum(month_days))
month_mid_days <- (cum_days[1:12] + cum_days[2:13]) / 2

# Interpolate to daily values

daily_ndvi <- approx(x = month_mid_days,
                     y = monthly_ndvi$ndvi,
                     xout = 1:365,
                     rule = 2)$y
mean_ndvi_year <- mean(daily_ndvi)

# Function to get NDVI for any time t (periodic)

ndvi_for_time <- function(t) {
  doy <- ((t - 1) %% 365) + 1
  approx(x = 1:365, y = daily_ndvi, xout = doy, rule = 2)$y
}

#  POPULATION & DEMOGRAPHY
N <- 283
mu <- 1 / (65 * 365)      # Natural death / replacement rate (per day); changed the life span to 65 years

# DISEASE PARAMETERS (HUMANS)
eh <- 1 / 4                # Incubation rate (mean 4 days)
ch <- 1 / 10               # Recovery rate (mean 10 days infectious)
lh <- 0.01                 # Disease-induced mortality per day
omega <- 1 / (3 * 365)     # Waning immunity rate (immunity ~3 years)

# LIVESTOCK-TO-HUMAN SPILLOVER
Il <- 20                   # Infected livestock (constant)
L  <- 100                  # Total livestock
klh <- 0.03                # Spillover transmission coefficient
lambda_mean <- klh * (Il / L)   # Baseline force of infection at mean NDVI

# SURVEILLANCE SCENARIOS 
# Detection occurs ONLY from infectious compartment; does NOT affect transmission.

scenarios <- list(
  Status_Quo = list(
    testing_rate = 0.001,    # Very low testing (only very ill suspects)
    Se_I = 0.74),               # Moderate sensitivity
  Improved_Surveillance = list(
    testing_rate = 0.05,      # Active screening of high‑risk group
    Se_I = 0.99))              # High sensitivity

# INITIAL CONDITIONS
# Endemic equilibrium without seasonality
seir_equilibrium <- function(x) {
  S <- x[1]; E <- x[2]; I <- x[3]; R <- x[4]
  Ph <- mu * N + lh * I
  dS <- Ph - lambda_mean * S - mu * S + omega * R
  dE <- lambda_mean * S - (eh + mu) * E
  dI <- eh * E - (mu + lh + ch) * I
  dR <- ch * I - mu * R - omega * R
  c(dS, dE, dI, dR)
}

EE <- multiroot(f = seir_equilibrium, start = c(50, 5, 5, 223))$root
names(EE) <- c("S", "E", "I", "R")
# EE: S=41.84, E=1.004, I=2.281, R=237.87

initial_state <- c(EE, Detected_I = 0)   # cumulative detections (state variable)
names(initial_state) <- c("S", "E", "I", "R", "Detected_I")

#  MODEL FUNCTION 
seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # NDVI‑modulated force of infection
    ndvi <- ndvi_for_time(time)
    lambda_t <- lambda_mean * (ndvi / mean_ndvi_year)
    lambda_t <- max(lambda_t, 0)
    
    # Population replacement
    Ph <- mu * N + lh * I
    
    # SEIR dynamics (no detection feedback)
    dS <- Ph - lambda_t * S - mu * S + omega * R
    dE <- lambda_t * S - (eh + mu) * E
    dI <- eh * E - (mu + lh + ch) * I
    dR <- ch * I - mu * R - omega * R
    
    # Cumulative detection events from infectious compartment only; I just wanted to start small with only I compartments but we can extend this with IgM, IgG diagnostics and see which one to use in such a circulating case
    
    dDetected_I <- testing_rate * Se_I * I
    
    list(c(dS, dE, dI, dR, dDetected_I))
  })
}

# SIMULATION SETTINGS: Running for endemic equilibrium

times <- seq(0, 365 * 50, by = 1)   # 50 years
burnin <- 365 * 40                   # Use last 10 years for analysis

#  RUN SCENARIOS

results <- list()

for (sc in names(scenarios)) {
  params <- c(
    mu = mu, eh = eh, ch = ch, lh = lh,
    lambda_mean = lambda_mean, omega = omega,
    testing_rate = scenarios[[sc]]$testing_rate,
    Se_I = scenarios[[sc]]$Se_I)
  
  out <- ode(y = initial_state,
             times = times,
             func = seir_model,
             parms = params,
             method = "lsoda")
  
  df <- as.data.frame(out)
  df$Scenario <- sc
  results[[sc]] <- df
}

results_all <- bind_rows(results)

# SEIR DYNAMICS
# Since dynamics are identical for both scenarios, we plot Status Quo as example.
seir_long <- results_all %>%
  filter(Scenario == "Status_Quo") %>%
  select(time, S, E, I, R) %>%
  pivot_longer(cols = c(S, E, I, R), names_to = "Compartment", values_to = "Count")

p1 <- ggplot(seir_long, aes(x = time / 365, y = Count, color = Compartment)) +
  geom_line(linewidth = 0.8) +
  labs(title = "SEIR Dynamics (Status Quo Scenario)",
       x = "Time (years)", y = "Number of individuals") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1)

# INFECTIOUS OVER TIME (This is the core for our research question)

p2 <- ggplot(results_all, aes(x = time / 365, y = I, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  labs(title = "Infectious Individuals Over Time",
       x = "Time (years)", y = "Number of infectious individuals") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2)

# MONTHLY DETECTED CASES; I intentionally did it to see the risk of RVF among abattoir works over time
# Extract last 365 days (endemic year)
end_time <- max(times)
start_analysis <- end_time - 364
results_endemic <- results_all %>%
  filter(time >= start_analysis) %>%
  arrange(Scenario, time)

# Compute daily new detections (difference in cumulative Detected_I)

results_endemic <- results_endemic %>%
  group_by(Scenario) %>%
  mutate(new_detections = Detected_I - lag(Detected_I, default = first(Detected_I))) %>%
  ungroup()

# Assign month (1-12) for the endemic year

results_endemic <- results_endemic %>%
  mutate(month = ((time - start_analysis) %/% (365/12)) + 1,
         month = pmin(month, 12))

# Sum new detections by month and scenario

monthly_detected <- results_endemic %>%
  group_by(Scenario, month) %>%
  summarise(Detected_I = sum(new_detections), .groups = "drop")

# Annual totals
annual_detected <- monthly_detected %>%
  group_by(Scenario) %>%
  summarise(Annual_Detected = sum(Detected_I))
print(annual_detected)

p3 <- ggplot(monthly_detected, aes(x = month, y = Detected_I, fill = Scenario)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Detected RVF Cases (Infectious Compartment Only)",
       subtitle = "Endemic equilibrium – last 12 months of 50‑year simulation",
       x = "Month", y = "Number of detected cases") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p3)

# SENSITIVITY ANALYSIS 
# Vary spillover coefficient klh by ±20%
klh_values <- lambda_mean * c(0.8, 1.0, 1.2)
sensitivity_results <- list()

for (klh_i in klh_values) {
  params <- c(
    mu = mu, eh = eh, ch = ch, lh = lh,
    lambda_mean = klh_i, omega = omega,
    testing_rate = scenarios$Status_Quo$testing_rate,
    Se_I = scenarios$Status_Quo$Se_I)
  
  out <- ode(y = initial_state,
             times = times,
             func = seir_model,
             parms = params,
             method = "lsoda")
  df <- as.data.frame(out)
  df$klh <- klh_i
  sensitivity_results[[paste0("klh_", klh_i)]] <- df
}

sensitivity_all <- bind_rows(sensitivity_results)

p4 <- ggplot(sensitivity_all %>% filter(time <= 50*365),
             aes(x = time/365, y = I, color = as.factor(klh))) +
  geom_line(linewidth = 0.8) +
  labs(title = "Sensitivity to Spillover Coefficient (klh)",
       x = "Time (years)", y = "Infectious Individuals",
       color = expression(lambda[mean])) +
  theme_minimal()
print(p4)

# PANEL FOR S, R, AND E+I
# Use Status Quo data (identical dynamics)
data_sq <- results_all %>% filter(Scenario == "Status_Quo")

# Top plot: Susceptible
p_S <- ggplot(data_sq, aes(x = time / 365, y = S)) +
  geom_line(color = "blue", linewidth = 0.8) +
  labs(title = "Susceptible (S)", x = "Time (years)", y = "Count") +
  theme_minimal()

# Middle plot: Recovered
p_R <- ggplot(data_sq, aes(x = time / 365, y = R)) +
  geom_line(color = "green", linewidth = 0.8) +
  labs(title = "Recovered (R)", x = "Time (years)", y = "Count") +
  theme_minimal()

# Bottom plot: Exposed and Infectious together
p_EI <- ggplot(data_sq, aes(x = time / 365)) +
  geom_line(aes(y = E, color = "Exposed"), linewidth = 0.8) +
  geom_line(aes(y = I, color = "Infectious"), linewidth = 0.8) +
  scale_color_manual(values = c("Exposed" = "orange", "Infectious" = "red")) +
  labs(title = "Exposed (E) and Infectious (I)", x = "Time (years)", y = "Count",
       color = "Compartment") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine using patchwork
p5 <- (p_S / p_R / p_EI) + plot_annotation(title = "SEIR Compartments (Status Quo)")
print(p5)
