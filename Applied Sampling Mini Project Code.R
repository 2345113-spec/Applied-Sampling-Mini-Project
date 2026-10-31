# Applied Sampling: Employment confidence among final-year students
# Wits and UJ

library(survey)

# SECTION B: SAMPLE SIZE

N_wits <- 5000
N_uj <- 5000
N <- N_wits + N_uj

p_wits <- 0.75
p_uj <- 0.65
P0 <- 0.70
z <- 1.96
d <- 0.05 * P0

# Neyman allocation with equal sampling costs
S_wits <- sqrt(p_wits * (1 - p_wits))
S_uj <- sqrt(p_uj * (1 - p_uj))
W_wits <- N_wits / N
W_uj <- N_uj / N

A <- W_wits * S_wits + W_uj * S_uj
F <- W_wits^2 * S_wits^2 / N_wits +
  W_uj^2 * S_uj^2 / N_uj

n_required <- ceiling(z^2 * A^2 / (d^2 + z^2 * F))

# Allocation used for the study
n_wits <- 290
n_uj <- 320
n <- n_wits + n_uj

# Check the sample size and allocation
n_required
n
n_wits
n_uj


# Sample size sensitivity plots
# The same Neyman formula is used with different planning values.
sample_size <- function(p, margin) {
  p1 <- p + 0.05
  p2 <- p - 0.05
  s1 <- sqrt(p1 * (1 - p1))
  s2 <- sqrt(p2 * (1 - p2))
  a <- W_wits * s1 + W_uj * s2
  f <- W_wits^2 * s1^2 / N_wits + W_uj^2 * s2^2 / N_uj
  ceiling(z^2 * a^2 / (margin^2 + z^2 * f))
}

# How sample size changes with the precision requirement
margins <- seq(0.01, 0.15, by = 0.005)
needed <- sapply(margins, function(r) sample_size(P0, r * P0))
plot(100 * margins, needed, type = "l", lwd = 2,
     xlab = "Margin of error (% of baseline proportion)",
     ylab = "Required sample size",
     main = "Sample size and required precision")
abline(v = 5, lty = 2)

# Sensitivity to different assumed population proportions
fractions <- c(0.25, 0.50, 0.75, 1.00)
assumed_p <- fractions * P0
sizes <- sapply(assumed_p, function(p) {
  sapply(margins, function(r) sample_size(p, r * p))
})
matplot(100 * margins, sizes, type = "l", lty = 1, lwd = 2,
        xlab = "Margin of error (% of assumed proportion)",
        ylab = "Required sample size",
        main = "Sample size sensitivity")
legend("topright", legend = paste0(100 * fractions, "% of baseline"),
       col = 1:4, lty = 1, bty = "n")

# Effect of the assumed proportion at a fixed margin of 0.035
p_values <- seq(0.10, 0.90, by = 0.01)
n_values <- sapply(p_values, function(p) sample_size(p, d))
plot(p_values, n_values, type = "l", lwd = 2,
     xlab = "Assumed population proportion",
     ylab = "Required sample size",
     main = "Sample size across assumed proportions")
abline(v = 0.50, lty = 2)

# SECTION C: GENERATE THE POPULATION AND SELECT A SAMPLE

set.seed(17)

population <- data.frame(
  student_id = 1:N,
  university = c(rep("Wits", N_wits), rep("UJ", N_uj))
)

population$field_of_study <- sample(
  c("Commerce", "Science", "Engineering", "Humanities", "Education"),
  N, replace = TRUE
)
population$internship <- rbinom(N, 1, 0.35)
population$work_experience <- rbinom(N, 1, 0.40)
population$career_programme <- rbinom(N, 1, 0.45)
population$confidence_probability <- ifelse(
  population$university == "Wits", p_wits, p_uj
)
population$employment_confidence <- rbinom(
  N, 1, population$confidence_probability
)

# Select students without replacement from each university
set.seed(17)
wits <- subset(population, university == "Wits")
uj <- subset(population, university == "UJ")

sample_wits <- wits[sample(1:N_wits, n_wits), , drop = FALSE]
sample_uj <- uj[sample(1:N_uj, n_uj), , drop = FALSE]
sample_data <- rbind(sample_wits, sample_uj)

sample_data$inclusion_probability <- ifelse(
  sample_data$university == "Wits", n_wits / N_wits, n_uj / N_uj
)
sample_data$sampling_weight <- 1 / sample_data$inclusion_probability

#population and selected sample are saved
write.csv(population, "career_readiness_population.csv", row.names = FALSE)
write.csv(sample_data, "career_readiness_sample.csv", row.names = FALSE)

P_true <- mean(population$employment_confidence)
# Check the generated population and selected sample
P_true
table(sample_data$university)

# SECTION D: ANALYSE THE SELECTED SAMPLE

sample_data$fpc <- ifelse(sample_data$university == "Wits", N_wits, N_uj)

design <- svydesign(
  ids = ~1,
  strata = ~university,
  weights = ~sampling_weight,
  fpc = ~fpc,
  data = sample_data
)

result <- svymean(~employment_confidence, design)
P_hat <- as.numeric(coef(result))
variance <- as.numeric(vcov(result))
se <- as.numeric(SE(result))

estimation_results <- data.frame(
  Estimate = P_hat,
  Variance = variance,
  SE = se,
  Lower_CI = P_hat - z * se,
  Upper_CI = P_hat + z * se,
  CV = se / P_hat
)

weights_table <- data.frame(
  University = c("Wits", "UJ"),
  Population = c(N_wits, N_uj),
  Sample = c(n_wits, n_uj),
  Probability = c(n_wits / N_wits, n_uj / N_uj),
  Weight = c(N_wits / n_wits, N_uj / n_uj)
)

comparison <- data.frame(
  Method = c("True population", "Unweighted sample", "Stratified estimate"),
  Proportion = c(P_true, mean(sample_data$employment_confidence), P_hat)
)

# Compare with one SRSWOR sample of the same size
set.seed(17)
srs_sample <- population[sample(1:N, n), , drop = FALSE]
srs_sample$fpc <- N
srs_design <- svydesign(ids = ~1, fpc = ~fpc, data = srs_sample)
srs_result <- svymean(~employment_confidence, srs_design)

precision_table <- data.frame(
  Design = c("Stratified", "SRSWOR"),
  Estimate = c(P_hat, as.numeric(coef(srs_result))),
  SE = c(se, as.numeric(SE(srs_result)))
)

srs_variance <- (1 - n / N) * var(population$employment_confidence) / n
deff <- variance / srs_variance
required_margin <- d
achieved_margin <- z * se

# Review the sample estimates and design comparison
estimation_results
weights_table
comparison
precision_table
deff
required_margin
achieved_margin
achieved_margin <= required_margin

# SECTION E: REPEATED SAMPLING

set.seed(17)
B <- 1000

strat_results <- data.frame(
  Estimate = numeric(B), SE = numeric(B),
  Lower = numeric(B), Upper = numeric(B)
)
srs_results <- strat_results

for (b in 1:B) {
  # Draw a stratified sample
  sw <- wits[sample(1:N_wits, n_wits), , drop = FALSE]
  su <- uj[sample(1:N_uj, n_uj), , drop = FALSE]

  p_strat <- W_wits * mean(sw$employment_confidence) +
    W_uj * mean(su$employment_confidence)

  v_strat <- W_wits^2 * (1 - n_wits / N_wits) *
    var(sw$employment_confidence) / n_wits +
    W_uj^2 * (1 - n_uj / N_uj) *
    var(su$employment_confidence) / n_uj

  se_strat <- sqrt(v_strat)
  strat_results[b, ] <- c(
    p_strat, se_strat,
    p_strat - z * se_strat,
    p_strat + z * se_strat
  )

  # Draw an SRSWOR sample
  srs <- population[sample(1:N, n), , drop = FALSE]
  p_srs <- mean(srs$employment_confidence)
  se_srs <- sqrt((1 - n / N) * var(srs$employment_confidence) / n)

  srs_results[b, ] <- c(
    p_srs, se_srs,
    p_srs - z * se_srs,
    p_srs + z * se_srs
  )
}

# Calculate the performance measures for each design
evaluate <- function(results) {
  bias <- mean(results$Estimate) - P_true
  mse <- mean((results$Estimate - P_true)^2)
  coverage <- mean(results$Lower <= P_true & results$Upper >= P_true)

  data.frame(
    Bias = bias,
    Relative_Bias = 100 * bias / P_true,
    Empirical_SE = sd(results$Estimate),
    Average_SE = mean(results$SE),
    MSE = mse,
    RMSE = sqrt(mse),
    Coverage = coverage,
    Average_Width = mean(results$Upper - results$Lower),
    Coverage_MCSE = sqrt(coverage * (1 - coverage) / B)
  )
}

comparison_E <- rbind(evaluate(strat_results), evaluate(srs_results))
comparison_E$Design <- c("Stratified", "SRSWOR")
comparison_E <- comparison_E[, c("Design", names(comparison_E)[names(comparison_E) != "Design"])]

empirical_deff <- var(strat_results$Estimate) / var(srs_results$Estimate)
estimated_deff <- mean(strat_results$SE^2) / mean(srs_results$SE^2)

# Review the Monte Carlo results
comparison_E
empirical_deff
estimated_deff

# SECTION E: PLOTS
# Red lines show the true proportion in the simulated population.
hist(strat_results$Estimate, breaks = 25,
     main = "Stratified sampling estimates",
     xlab = "Estimated employment confidence",
     col = "lightblue")
abline(v = P_true, col = "red", lwd = 2)

hist(srs_results$Estimate, breaks = 25,
     main = "SRSWOR estimates",
     xlab = "Estimated employment confidence",
     col = "lightgreen")
abline(v = P_true, col = "red", lwd = 2)

