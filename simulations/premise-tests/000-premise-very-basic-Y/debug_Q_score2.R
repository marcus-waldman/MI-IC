library(MASS)
source("simulations/000-premise-very-basic-Y/mvn_analytics.R")

set.seed(123)
p <- 3
N <- 50
mu <- rep(0, p)
Sigma <- toeplitz(0.5^(0:(p-1)))

Y_complete <- mvrnorm(N, mu, Sigma)
Y_miss <- Y_complete
Y_miss[1:10, 3] <- NA

# Compute E[S_mu] manually
sum_E_Y <- rep(0, p)
sum_E_YYt <- matrix(0, p, p)

for (i in 1:N) {
  obs_idx <- which(!is.na(Y_miss[i, ]))
  mis_idx <- which(is.na(Y_miss[i, ]))

  if (length(mis_idx) == 0) {
    Y_i <- Y_miss[i, ]
    sum_E_Y <- sum_E_Y + Y_i
    sum_E_YYt <- sum_E_YYt + outer(Y_i, Y_i)
  } else {
    Y_o <- Y_miss[i, obs_idx]

    # Partition
    mu_o <- mu[obs_idx]
    mu_m <- mu[mis_idx]
    Sigma_oo <- Sigma[obs_idx, obs_idx, drop = FALSE]
    Sigma_mm <- Sigma[mis_idx, mis_idx, drop = FALSE]
    Sigma_mo <- Sigma[mis_idx, obs_idx, drop = FALSE]
    Sigma_om <- Sigma[obs_idx, mis_idx, drop = FALSE]

    # Conditional moments
    Sigma_oo_inv <- solve(Sigma_oo)
    mu_m_given_o <- mu_m + Sigma_mo %*% Sigma_oo_inv %*% (Y_o - mu_o)
    Sigma_m_given_o <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om

    # E[Y_i]
    E_Y_i <- rep(0, p)
    E_Y_i[obs_idx] <- Y_o
    E_Y_i[mis_idx] <- as.vector(mu_m_given_o)
    sum_E_Y <- sum_E_Y + E_Y_i

    # E[Y_i Y_i']
    E_YYt_i <- matrix(0, p, p)
    E_YYt_i[obs_idx, obs_idx] <- outer(Y_o, Y_o)
    E_YYt_i[obs_idx, mis_idx] <- outer(Y_o, as.vector(mu_m_given_o))
    E_YYt_i[mis_idx, obs_idx] <- t(E_YYt_i[obs_idx, mis_idx])
    E_YYt_i[mis_idx, mis_idx] <- Sigma_m_given_o + outer(as.vector(mu_m_given_o), as.vector(mu_m_given_o))
    sum_E_YYt <- sum_E_YYt + E_YYt_i
  }
}

E_Y_bar <- sum_E_Y / N
E_S_mu <- sum_E_YYt / N - outer(mu, E_Y_bar) - outer(E_Y_bar, mu) + outer(mu, mu)

cat("E_S_mu:\n")
print(round(E_S_mu, 4))

# Score matrix before vech conversion
Sigma_inv <- solve(Sigma)
score_Sigma_mat <- (N / 2) * Sigma_inv %*% (E_S_mu - Sigma) %*% Sigma_inv

cat("\nScore matrix (before vech):\n")
print(round(score_Sigma_mat, 4))

# Convert using D_p'
D_p <- duplication_matrix(p)
score_vech_Sigma <- as.vector(t(D_p) %*% as.vector(score_Sigma_mat))

cat("\nScore vech(Sigma) (after D_p'):\n")
print(round(score_vech_Sigma, 4))

# Compare with numerical for Sigma components only
theta <- c(mu, vech(Sigma))
eps <- 1e-6
numerical_Sigma_scores <- numeric(6)

for (k in 1:6) {
  j <- k + 3  # Skip the first 3 (mu components)

  theta_plus <- theta
  theta_minus <- theta
  theta_plus[j] <- theta_plus[j] + eps
  theta_minus[j] <- theta_minus[j] - eps

  Sigma_plus <- vech_inv(theta_plus[4:9], p)
  Sigma_minus <- vech_inv(theta_minus[4:9], p)

  Q_plus <- compute_Q_analytic(Y_miss, mu, Sigma_plus)
  Q_minus <- compute_Q_analytic(Y_miss, mu, Sigma_minus)

  numerical_Sigma_scores[k] <- (Q_plus - Q_minus) / (2 * eps)
}

cat("\nNumerical Sigma scores:\n")
print(round(numerical_Sigma_scores, 4))

cat("\nDifference:\n")
print(round(score_vech_Sigma - numerical_Sigma_scores, 4))
