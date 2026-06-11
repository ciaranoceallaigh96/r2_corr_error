#SCRIPT TO COMPARE PROPERTIES OF CORRELATION-BASED R2 VS. ERROR-BASED R2 IN TERMS OF HERITABILITY BENCHMARKING
#R2_corr is useful for heritability benchmarking
#R2_error is useful for detecting calibration problems
#
#A PRS may contain the right genetic signal, but be badly calibrated.
#R2_corr is unchanged by intercept/scale errors.
#R2_error is penalized by intercept/scale errors.
#
#R2_corr is desirable for variance-capture benchmarking.
#R2_error is desirable for checking raw prediction calibration.



set.seed(123)

n <- 50000

#Target heritabilities
target_h2 <- 0.40
target_H2 <- 0.60

#define R2 metrics
r2_corr <- function(y, yhat) {
  cor(y, yhat)^2
}

r2_error <- function(y, yhat) {
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}




#Additive genetic component
A <- rnorm(n, mean = 0, sd = 1)

#Non-additive genetic component (independent of A)
NonAdd <- rnorm(n, mean = 0, sd = 1)

#Var(A) = target_h2
#Var(NonAdd) = target_H2 - target_h2

A <- scale(A)[, 1] #mean = 0, variance = 1

A <- A * sqrt(target_h2) #variance = target_h2

NonAdd <- scale(NonAdd)[, 1] #mean = 0, variance = 1

NonAdd <- NonAdd * sqrt(target_H2 - target_h2) #variance = non-additive component

#Total genetic component
G <- A + NonAdd

#Environmental component
E <- rnorm(n, mean = 0, sd = sqrt(1 - target_H2))

#Phenotype
Phen <- G + E

cat("\nVariance components:\n")
print(round(data.frame(
  Var_A     = var(A),
  Var_NonAdd = var(NonAdd),
  Var_E     = var(E),
  Var_Phen  = var(Phen),
  Sum_components = var(A) + var(NonAdd) + var(E)
), 4))



#heritabilities
h2 <- var(A) / var(Phen)
H2 <- var(G) / var(Phen)

cat("\nHeritability benchmarks:\n")
print(round(data.frame(
  h2 = h2,
  H2 = H2
), 4))


#Additive PRS
add_prs_good <- A

#Non-additive PRS:
#represents a non-linear model capturing the full genetic component
nonadd_prs_good <- G


#MIS-CALIBRATION
#scores contain the same signal as the good PRS but on wrong scale and wrong intercept
#a well-calibrated PRS should have intercept near 0 and slope near 1

#add intercept of 1 and slope/scaling factor of 2
add_prs_bad <- 1 + 2 * add_prs_good
nonadd_prs_bad <- 1 + 2 * nonadd_prs_good

#NOISE + ERROR for atteunated score
#noise standard deviation is set to equal to the square root of each signal variance for each score
#signal and noise now contribute equally to score variation and expected signal is reduced by approximately 50%

sd_add_noise <- sqrt(var(A))
sd_nonadd_noise <- sqrt(var(G))

add_prs_noisy <- A + rnorm(n, mean = 0, sd = sd_add_noise)
nonadd_prs_noisy <- G + rnorm(n, mean = 0, sd = sd_nonadd_noise)

results <- data.frame(
  scenario = c(
    "Additive PRS, well calibrated",
    "Additive PRS, badly calibrated",
    "Additive PRS, noisy",
    "Non-additive PRS, well calibrated",
    "Non-additive PRS, badly calibrated",
    "Non-additive PRS, noisy"
  ),

  benchmark = c(
    "h2",
    "h2",
    "h2",
    "H2",
    "H2",
    "H2"
  ),

  benchmark_value = c(
    h2,
    h2,
    h2,
    H2,
    H2,
    H2
  ),

  R2_corr = c(
    r2_corr(Phen, add_prs_good),
    r2_corr(Phen, add_prs_bad),
    r2_corr(Phen, add_prs_noisy),
    r2_corr(Phen, nonadd_prs_good),
    r2_corr(Phen, nonadd_prs_bad),
    r2_corr(Phen, nonadd_prs_noisy)
  ),

  R2_error = c(
    r2_error(Phen, add_prs_good),
    r2_error(Phen, add_prs_bad),
    r2_error(Phen, add_prs_noisy),
    r2_error(Phen, nonadd_prs_good),
    r2_error(Phen, nonadd_prs_bad),
    r2_error(Phen, nonadd_prs_noisy)
  )
)

results$difference <- results$R2_corr - results$R2_error

results_plot <- results

num_cols <- sapply(results, is.numeric)
results[num_cols] <- round(results[num_cols], 4)

cat("\nR2 comparison:\n")
print(results)





#Scatter plots of Phenotype vs PRS [all plots use non-additive scenario]
prs_list <- list(
  "Additive good"      = add_prs_good,
  "Additive bad"       = add_prs_bad,
  "Additive noisy"     = add_prs_noisy,
  "Non-additive good"  = nonadd_prs_good,
  "Non-additive bad"   = nonadd_prs_bad,
  "Non-additive noisy" = nonadd_prs_noisy
)

#Sample points for visualization purposes
set.seed(123)
idx <- sample(seq_len(n), 3000)


#A: well-calibrated prediction score 
#C: badly calibrated prediction score with free axes
#
#Dashed line = perfect calibration: y = x
#Red line = fitted auxiliary regression line: phenotype ~ prediction score

plot_cases <- list(
  "panel_a_well_calibrated_matched_axes.png" = list(
    x_full = nonadd_prs_good,
    x_plot = nonadd_prs_good[idx],
    y_plot = Phen[idx],
    matched_axes = TRUE
  ),
  "panel_c_miscalibrated_free_axes.png" = list(
    x_full = nonadd_prs_bad,
    x_plot = nonadd_prs_bad[idx],
    y_plot = Phen[idx],
    matched_axes = FALSE
  )
)

for (file_name in names(plot_cases)) {
  
  png(
    filename = file_name,
    width = 6,
    height = 6,
    units = "in",
    res = 600
  )
  
  par(mar = c(4.5, 4.5, 1, 1))
  
  x_full <- plot_cases[[file_name]]$x_full
  y_full <- Phen
  
  x <- plot_cases[[file_name]]$x_plot
  y <- plot_cases[[file_name]]$y_plot
  
  fit_full <- lm(y_full ~ x_full)
  fit_sample <- lm(y ~ x)
  
  r2c <- r2_corr(y_full, x_full)
  r2e <- r2_error(y_full, x_full)
  
  if (plot_cases[[file_name]]$matched_axes) {
    
    axis_lim <- range(c(x, y), na.rm = TRUE)
    
    plot(
      x, y,
      pch = 16,
      cex = 0.45,
      col = rgb(0, 0, 0, 0.2),
      xlab = "Prediction Score",
      ylab = "Observed Phenotype",
      main = "",
      xlim = axis_lim,
      ylim = axis_lim,
      asp = 1
    )
    
  } else {
    
    plot(
      x, y,
      pch = 16,
      cex = 0.45,
      col = rgb(0, 0, 0, 0.2),
      xlab = "Prediction Score",
      ylab = "Observed Phenotype",
      main = ""
    )
  }
  
  #Perfect calibration line
  abline(a = 0, b = 1, lty = 2, lwd = 2)
  
  #Fitted auxiliary regression line
  abline(fit_sample, col = "red", lwd = 2)
  
  legend(
    "topleft",
    legend = c(
      bquote(R[corr]^2 == .(format(round(r2c, 2), nsmall = 2))),
      bquote(R[error]^2 == .(format(round(r2e, 2), nsmall = 2)))
    ),
    bty = "n",
    cex = 0.8,
    y.intersp = 1.6
  )
  
  legend(
    "bottomright",
    legend = c("Perfect calibration", "Fitted auxiliary regression"),
    lty = c(2, 1),
    col = c("black", "red"),
    lwd = 2,
    bty = "n",
    cex = 0.8
  )
  
  dev.off()
}

cat("\nSaved separate panel images:\n")
cat("panel_a_well_calibrated_matched_axes.png\n")
cat("panel_c_miscalibrated_free_axes.png\n")



#B: attenuated non-additive prediction score with matched axes

file_name <- "panel_b_nonadditive_attenuated_matched_axes.png"

png(
  filename = file_name,
  width = 6,
  height = 6,
  units = "in",
  res = 600
)

par(mar = c(4.5, 4.5, 1, 1))

x_full <- nonadd_prs_noisy
y_full <- Phen

x <- nonadd_prs_noisy[idx]
y <- Phen[idx]

fit_sample <- lm(y ~ x)

r2c <- r2_corr(y_full, x_full)
r2e <- r2_error(y_full, x_full)

axis_lim <- range(c(x, y), na.rm = TRUE)

plot(
  x, y,
  pch = 16,
  cex = 0.45,
  col = rgb(0, 0, 0, 0.2),
  xlab = "Prediction Score",
  ylab = "Observed Phenotype",
  main = "",
  xlim = axis_lim,
  ylim = axis_lim,
  asp = 1
)

#Perfect calibration line
abline(a = 0, b = 1, lty = 2, lwd = 2)

#Fitted auxiliary regression line
abline(fit_sample, col = "red", lwd = 2)

legend(
  "topleft",
  legend = c(
    bquote(R[corr]^2 == .(format(round(r2c, 2), nsmall = 2))),
    bquote(R[error]^2 == .(format(round(r2e, 2), nsmall = 2)))
  ),
  bty = "n",
  cex = 0.8,
  y.intersp = 1.6
)

legend(
  "bottomright",
  legend = c("Perfect calibration", "Fitted auxiliary regression"),
  lty = c(2, 1),
  col = c("black", "red"),
  lwd = 2,
  bty = "n",
  cex = 0.8
)

dev.off()

cat("\nSaved attenuated non-additive panel image:\n")
cat(file_name, "\n")
