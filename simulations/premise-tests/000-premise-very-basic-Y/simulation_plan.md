# Simulation Plan: Term Decomposition Verification for Saturated MVN

## 1. Objective

Empirically verify the Taylor expansion of imputation bias to identify the source of the observed 2×tr(RIV) discrepancy. The theoretical derivation predicts imputation bias ≈ tr(RIV), but empirical results show bias ≈ 2×tr(RIV).

## 2. Data Generating Mechanism

Generate complete data from a multivariate normal distribution with saturated (unstructured) covariance:

$$
Y_i \stackrel{iid}{\sim} \mathcal{N}_p(\mu^*, \Sigma^*), \quad i = 1, \ldots, N \tag{1}
$$

where:
- $\mu^* \in \mathbb{R}^p$ is the true mean vector
- $\Sigma^* \in \mathbb{R}^{p \times p}$ is the true (positive definite) covariance matrix
- $N$ is the sample size
- $p$ is the dimension

## 3. Missingness Mechanism

Impose missingness on the complete data $Y$ to obtain observed data $Y_{obs}$ and missing data $Y_{mis}$:

$$
Y = (Y_{obs}, Y_{mis}) \tag{2}
$$

Use a monotone MCAR pattern with missing rate $\gamma \in (0,1)$.

## 4. Parameter Space

For the saturated MVN model, the parameter vector is:

$$
\theta = (\mu, \text{vech}(\Sigma)) \in \mathbb{R}^Q \tag{3}
$$

where $\text{vech}(\cdot)$ denotes the half-vectorization (lower triangle including diagonal), and the total number of parameters is:

$$
Q = p + \frac{p(p+1)}{2} \tag{4}
$$

## 5. Imputation Procedure

Use Amelia to generate $M$ imputations from:

$$
Y_{mis}^{(m)} \sim f(Y_{mis} \mid Y_{obs}, \hat{\theta}_{MLE}), \quad m = 1, \ldots, M \tag{5}
$$

where $\hat{\theta}_{MLE}$ is the observed-data MLE obtained via EM algorithm (Amelia with `boot.type = "none"`).

Each imputation yields a completed dataset:

$$
Y^{(m)} = (Y_{obs}, Y_{mis}^{(m)}) \tag{6}
$$

## 6. Complete-Data Log-Likelihood

For the saturated MVN model, the complete-data log-likelihood is:

$$
\ell_{com}(\mu, \Sigma \mid Y) = -\frac{Np}{2}\log(2\pi) - \frac{N}{2}\log|\Sigma| - \frac{N}{2}\text{tr}(\Sigma^{-1} S_\mu) \tag{7}
$$

where the sample covariance centered at $\mu$ is:

$$
S_\mu = \frac{1}{N} \sum_{i=1}^{N} (Y_i - \mu)(Y_i - \mu)^\top \tag{8}
$$

## 7. Q-Function (Averaged Over Imputations)

The MI analogue of the Q-function is:

$$
\bar{Q}_{MI}(\theta) = \frac{1}{M} \sum_{m=1}^{M} \ell_{com}(\theta \mid Y^{(m)}) \tag{9}
$$

## 8. Taylor Expansion of Imputation Bias

The total imputation bias at the MLE is:

$$
\text{Total Bias} = \bar{Q}_{MI}(\hat{\theta}) - \ell_{com}(\hat{\theta}) \tag{10}
$$

This can be decomposed via Taylor expansion around the true parameter $\theta^*$:

$$
\bar{Q}_{MI}(\hat{\theta}) - \ell_{com}(\hat{\theta}) = \text{Term}_1 + \text{Term}_2 + \text{Term}_3 + O(N^{-3/2}) \tag{11}
$$

### Term 1: Log-likelihood difference at true parameters

$$
\text{Term}_1 = \bar{Q}_{MI}(\theta^*) - \ell_{com}(\theta^*) \tag{12}
$$

**Theoretical prediction:** Term 1 ≈ 0 (to first order)

### Term 2: Score difference × parameter error

$$
\text{Term}_2 = (\hat{\theta} - \theta^*)^\top \left[ \bar{S}_{MI}(\theta^*) - S_{com}(\theta^*) \right] \tag{13}
$$

where:
- $S_{com}(\theta) = \nabla_\theta \ell_{com}(\theta)$ is the complete-data score
- $\bar{S}_{MI}(\theta) = \frac{1}{M} \sum_{m=1}^{M} S_{com}^{(m)}(\theta)$ is the averaged score over imputations

**Theoretical prediction:** Term 2 ≈ tr(RIV)

### Term 3: Hessian quadratic form

$$
\text{Term}_3 = -\frac{1}{2} (\hat{\theta} - \theta^*)^\top \left[ \bar{J}_{MI}(\theta^*) - J_{com}(\theta^*) \right] (\hat{\theta} - \theta^*) \tag{14}
$$

where:
- $J_{com}(\theta) = -\nabla^2_\theta \ell_{com}(\theta)$ is the observed information (negative Hessian)
- $\bar{J}_{MI}(\theta) = \frac{1}{M} \sum_{m=1}^{M} J_{com}^{(m)}(\theta)$ is the averaged information

**Theoretical prediction:** Term 3 ≈ 0 (claimed negligible, $O(N^{-3/2})$)

## 9. Analytic Formulas for MVN

### 9.1 Score with respect to μ

$$
\frac{\partial \ell}{\partial \mu} = N \Sigma^{-1} (\bar{Y} - \mu) \tag{15}
$$

where $\bar{Y} = \frac{1}{N} \sum_{i=1}^{N} Y_i$ is the sample mean.

### 9.2 Score with respect to Σ (matrix form)

$$
\frac{\partial \ell}{\partial \Sigma} = -\frac{N}{2} \Sigma^{-1} + \frac{N}{2} \Sigma^{-1} S_\mu \Sigma^{-1} = \frac{N}{2} \Sigma^{-1} (S_\mu - \Sigma) \Sigma^{-1} \tag{16}
$$

### 9.3 Score with respect to vech(Σ)

Using the duplication matrix $D_p$ that maps $\text{vech}(A)$ to $\text{vec}(A)$:

$$
\frac{\partial \ell}{\partial \text{vech}(\Sigma)} = D_p^\top \text{vec}\left( \frac{\partial \ell}{\partial \Sigma} \right) \tag{17}
$$

### 9.4 Fisher Information for μ

$$
\mathcal{I}_{\mu\mu} = N \Sigma^{-1} \tag{18}
$$

### 9.5 Fisher Information for vech(Σ)

$$
\mathcal{I}_{\text{vech}(\Sigma)} = \frac{N}{2} D_p^\top (\Sigma^{-1} \otimes \Sigma^{-1}) D_p \tag{19}
$$

where $\otimes$ denotes the Kronecker product.

### 9.6 Cross-Information (μ and Σ)

$$
\mathcal{I}_{\mu, \text{vech}(\Sigma)} = 0 \tag{20}
$$

The mean and covariance parameters are orthogonal in the MVN model.

### 9.7 Complete-Data Observed Information (Hessian)

For MVN with known structure, the observed and expected information coincide:

$$
J_{com}(\theta) = \mathcal{I}_{com}(\theta) = \begin{pmatrix} N\Sigma^{-1} & 0 \\ 0 & \frac{N}{2} D_p^\top (\Sigma^{-1} \otimes \Sigma^{-1}) D_p \end{pmatrix} \tag{21}
$$

## 10. Rubin's Variance Components

For each imputed dataset $m$, compute the MLE:

$$
\hat{\theta}^{(m)} = (\bar{Y}^{(m)}, \text{vech}(S^{(m)})) \tag{22}
$$

where $S^{(m)} = \frac{1}{N} \sum_{i=1}^{N} (Y_i^{(m)} - \bar{Y}^{(m)})(Y_i^{(m)} - \bar{Y}^{(m)})^\top$.

### 10.1 Within-Imputation Variance

$$
\bar{W} = \frac{1}{M} \sum_{m=1}^{M} \text{Var}(\hat{\theta}^{(m)}) = \frac{1}{M} \sum_{m=1}^{M} \mathcal{I}_{com}^{-1}(\hat{\theta}^{(m)}) \tag{23}
$$

### 10.2 Between-Imputation Variance

$$
B = \frac{1}{M-1} \sum_{m=1}^{M} (\hat{\theta}^{(m)} - \bar{\theta})(\hat{\theta}^{(m)} - \bar{\theta})^\top \tag{24}
$$

where $\bar{\theta} = \frac{1}{M} \sum_{m=1}^{M} \hat{\theta}^{(m)}$.

### 10.3 Total Variance and RIV

$$
T = \bar{W} + \left(1 + \frac{1}{M}\right) B \tag{25}
$$

$$
\text{RIV} = \left(1 + \frac{1}{M}\right) \bar{W}^{-1} B \tag{26}
$$

$$
\text{tr}(\text{RIV}) = \text{trace}\left[ \left(1 + \frac{1}{M}\right) \bar{W}^{-1} B \right] \tag{27}
$$

## 11. Experimental Design

### Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| $N$ | 2000 | Large sample to minimize Var($\hat{\theta}$) |
| $p$ | 5 | Moderate dimension, Q = 20 parameters |
| $\gamma$ | 0.6 | High missingness to maximize tr(RIV) |
| $M$ | 100 | Sufficient imputations for stable estimates |
| $n_{reps}$ | Start with 1 | Single replication may suffice with large N |

### True Parameters

- $\mu^* = \mathbf{0}_p$
- $\Sigma^* = $ Toeplitz or CS structure with $\rho = 0.5$

## 12. Quantities to Compute

For each replication:

1. **Total Bias** (Eq. 10): $\bar{Q}_{MI}(\hat{\theta}) - \ell_{com}(\hat{\theta})$
2. **Term 1** (Eq. 12): $\bar{Q}_{MI}(\theta^*) - \ell_{com}(\theta^*)$
3. **Term 2** (Eq. 13): $(\hat{\theta} - \theta^*)^\top [\bar{S}_{MI}(\theta^*) - S_{com}(\theta^*)]$
4. **Term 3** (Eq. 14): $-\frac{1}{2}(\hat{\theta} - \theta^*)^\top [\bar{J}_{MI}(\theta^*) - J_{com}(\theta^*)](\hat{\theta} - \theta^*)$
5. **Taylor Sum**: Term 1 + Term 2 + Term 3
6. **Residual**: Total Bias - Taylor Sum
7. **tr(RIV)** (Eq. 27)

## 13. Expected Outcomes

| Scenario | Term 1 | Term 2 | Term 3 | Implication |
|----------|--------|--------|--------|-------------|
| Theory correct | ≈ 0 | ≈ tr(RIV) | ≈ 0 | Bias = tr(RIV) |
| Term 3 not negligible | ≈ 0 | ≈ tr(RIV) | ≈ tr(RIV) | Bias = 2×tr(RIV) |
| Term 2 error | ≈ 0 | ≈ 2×tr(RIV) | ≈ 0 | Derivation error in Step 14-16 |
| Term 1 contributes | ≈ tr(RIV) | ≈ tr(RIV) | ≈ 0 | First-order term not zero |

## 14. Validation Checks

1. **Taylor expansion validity**: |Residual| / |Total Bias| < 5%
2. **Stable ratios**: Term_k / tr(RIV) should be consistent across replications
3. **Sign consistency**: Total Bias should be positive (overfitting direction)
