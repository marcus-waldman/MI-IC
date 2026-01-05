library(MASS)
source("simulations/000-premise-very-basic-Y/mvn_analytics.R")

set.seed(123)
p <- 3
N <- 50
mu <- rep(0, p)
Sigma <- toeplitz(0.5^(0:(p-1)))

# First test on complete data (should match mvn_score)
Y_complete <- mvrnorm(N, mu, Sigma)

S_complete_analytic <- mvn_score(Y_complete, mu, Sigma)
S_Q_complete <- compute_Q_score(Y_complete, mu, Sigma)

cat("Complete data test:\n")
cat("mvn_score:     ", round(S_complete_analytic, 4), "\n")
cat("compute_Q_score:", round(S_Q_complete, 4), "\n")
cat("Difference:    ", round(S_complete_analytic - S_Q_complete, 6), "\n\n")

# Now test with missing data
Y_miss <- Y_complete
Y_miss[1:10, 3] <- NA

S_Q_analytic <- compute_Q_score(Y_miss, mu, Sigma)

# Numerical gradient
eps <- 1e-6
theta <- c(mu, vech(Sigma))
Q_len <- length(theta)
S_Q_numerical <- numeric(Q_len)

for (j in 1:Q_len) {
  theta_plus <- theta
  theta_minus <- theta

  theta_plus[j] <- theta_plus[j] + eps
  theta_minus[j] <- theta_minus[j] - eps

  if (j <= p) {
    mu_plus <- theta_plus[1:p]
    mu_minus <- theta_minus[1:p]
    Sigma_plus <- Sigma
    Sigma_minus <- Sigma
  } else {
    mu_plus <- mu
    mu_minus <- mu
    Sigma_plus <- vech_inv(theta_plus[(p+1):Q_len], p)
    Sigma_minus <- vech_inv(theta_minus[(p+1):Q_len], p)
  }

  Q_plus <- compute_Q_analytic(Y_miss, mu_plus, Sigma_plus)
  Q_minus <- compute_Q_analytic(Y_miss, mu_minus, Sigma_minus)

  S_Q_numerical[j] <- (Q_plus - Q_minus) / (2 * eps)
}

cat("Missing data test:\n")
cat("Analytic:  ", round(S_Q_analytic, 4), "\n")
cat("Numerical: ", round(S_Q_numerical, 4), "\n")
cat("Difference:", round(S_Q_analytic - S_Q_numerical, 4), "\n")
cat("Max diff:  ", max(abs(S_Q_analytic - S_Q_numerical)), "\n")
