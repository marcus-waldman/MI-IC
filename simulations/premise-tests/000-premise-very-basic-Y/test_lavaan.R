library(lavaan)
library(MASS)

set.seed(12345)
N <- 100
p <- 3
mu <- rep(0, p)
Sigma <- toeplitz(0.5^(0:(p-1)))
Y <- mvrnorm(N, mu, Sigma)
colnames(Y) <- paste0("Y", 1:p)

# Generate model syntax
var_names <- paste0("Y", 1:p)
lines <- character()
for (i in 1:p) {
  covs <- paste(var_names[i:p], collapse = " + ")
  lines <- c(lines, paste0(var_names[i], " ~~ ", covs))
}
model_syntax <- paste(lines, collapse = "\n")
cat("Model syntax:\n", model_syntax, "\n\n")

# Fit
fit <- sem(model_syntax, data = as.data.frame(Y), meanstructure = TRUE)

# Check what lavInspect returns
cat("Coefficient names:\n")
print(names(coef(fit)))

cat("\nlavInspect est:\n")
est <- lavInspect(fit, "est")
print(est)

cat("\nType of est:\n")
print(class(est))

cat("\nNames in est:\n")
print(names(est))
