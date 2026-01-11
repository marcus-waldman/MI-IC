# Analytic Q-Function for MVN with Missing Data

## Setup

Complete data: $Y_i \sim \mathcal{N}_p(\mu, \Sigma)$, $i = 1, \ldots, N$

Partition each observation into observed and missing components:
$$
Y_i = \begin{pmatrix} Y_{obs,i} \\ Y_{mis,i} \end{pmatrix}, \quad
\mu = \begin{pmatrix} \mu_{obs} \\ \mu_{mis} \end{pmatrix}, \quad
\Sigma = \begin{pmatrix} \Sigma_{obs} & \Sigma_{obs,mis} \\ \Sigma_{mis,obs} & \Sigma_{mis} \end{pmatrix}
$$

## Complete-Data Log-Likelihood

$$
\ell_{com}(\theta | Y) = -\frac{Np}{2}\log(2\pi) - \frac{N}{2}\log|\Sigma| - \frac{1}{2}\sum_{i=1}^{N}(Y_i - \mu)^\top \Sigma^{-1} (Y_i - \mu) \tag{1}
$$

## Conditional Distribution of Missing Data

Given observed data $Y_{obs,i}$ and parameters $\theta' = (\mu', \Sigma')$:
$$
Y_{mis,i} | Y_{obs,i}, \theta' \sim \mathcal{N}(\tilde{\mu}_i, \Sigma_{mis|obs}) \tag{2}
$$

where:
$$
\tilde{\mu}_i = \mu'_{mis} + \Sigma'_{mis,obs} (\Sigma'_{obs})^{-1} (Y_{obs,i} - \mu'_{obs}) \tag{3}
$$

$$
\Sigma_{mis|obs} = \Sigma'_{mis} - \Sigma'_{mis,obs} (\Sigma'_{obs})^{-1} \Sigma'_{obs,mis} \tag{4}
$$

**Key insight**: $\Sigma_{mis|obs}$ does not depend on $i$ — it's the same for all observations with the same missing pattern.

## Q-Function Definition

$$
Q(\theta | \theta') = E\left[\ell_{com}(\theta | Y) \,\Big|\, Y_{obs}, \theta'\right] \tag{5}
$$

## Deriving the Q-Function

The expectation affects only terms involving $Y_{mis}$. Using the identity for quadratic forms of Gaussian vectors:
$$
E[X^\top A X] = \mu_X^\top A \mu_X + \text{tr}(A \cdot \text{Var}(X)) \tag{6}
$$

For each observation $i$:
$$
E\left[(Y_i - \mu)^\top \Sigma^{-1} (Y_i - \mu) \,\Big|\, Y_{obs,i}, \theta'\right] = (\hat{Y}_i - \mu)^\top \Sigma^{-1} (\hat{Y}_i - \mu) + \text{tr}(\Sigma^{-1} C_i) \tag{7}
$$

where:
- $\hat{Y}_i = E[Y_i | Y_{obs,i}, \theta'] = (Y_{obs,i}, \tilde{\mu}_i)^\top$ is the conditionally imputed observation
- $C_i = \text{Cov}(Y_i | Y_{obs,i}, \theta')$ is block matrix:
$$
C_i = \begin{pmatrix} 0 & 0 \\ 0 & \Sigma_{mis|obs} \end{pmatrix} \tag{8}
$$

## The Q-Function Explicitly

$$
Q(\theta | \theta') = -\frac{Np}{2}\log(2\pi) - \frac{N}{2}\log|\Sigma| - \frac{1}{2}\sum_{i=1}^{N}\left[(\hat{Y}_i - \mu)^\top \Sigma^{-1} (\hat{Y}_i - \mu) + \text{tr}(\Sigma^{-1} C_i)\right] \tag{9}
$$

## Difference: Q-Function minus Complete-Data Log-Likelihood

$$
Q(\theta | \theta') - \ell_{com}(\theta) = -\frac{1}{2}\sum_{i=1}^{N}\left[\text{tr}(\Sigma^{-1} C_i) + (\hat{Y}_i - \mu)^\top \Sigma^{-1} (\hat{Y}_i - \mu) - (Y_i - \mu)^\top \Sigma^{-1} (Y_i - \mu)\right] \tag{10}
$$

Let $e_i = \hat{Y}_i - Y_i$ be the imputation error. Then:
$$
(\hat{Y}_i - \mu)^\top \Sigma^{-1} (\hat{Y}_i - \mu) = (Y_i - \mu + e_i)^\top \Sigma^{-1} (Y_i - \mu + e_i)
$$
$$
= (Y_i - \mu)^\top \Sigma^{-1} (Y_i - \mu) + 2e_i^\top \Sigma^{-1} (Y_i - \mu) + e_i^\top \Sigma^{-1} e_i \tag{11}
$$

Substituting:
$$
Q(\theta | \theta') - \ell_{com}(\theta) = -\frac{1}{2}\sum_{i=1}^{N}\left[\text{tr}(\Sigma^{-1} C_i) + 2e_i^\top \Sigma^{-1} (Y_i - \mu) + e_i^\top \Sigma^{-1} e_i\right] \tag{12}
$$

## Special Case: Evaluating at True Parameters

When $\theta = \theta' = \theta^*$ (true parameters), and we're computing the **expected** difference over the sampling distribution:

The imputation $\hat{Y}_i$ is drawn from the correct conditional distribution, so:
- $E[e_i] = 0$ (unbiased imputation)
- $E[e_i e_i^\top] = C_i = \text{blockdiag}(0, \Sigma_{mis|obs})$

Taking expectations over the data generating process:

$$
E[Q(\theta^* | \theta^*) - \ell_{com}(\theta^*)] = -\frac{1}{2}\sum_{i=1}^{N}\left[\text{tr}(\Sigma^{*-1} C_i) + \text{tr}(\Sigma^{*-1} C_i)\right] \tag{13}
$$

$$
= -\sum_{i=1}^{N} \text{tr}(\Sigma^{*-1} C_i) \tag{14}
$$

## Simplification for Homogeneous Missing Pattern

If all $N_{mis}$ incomplete observations have the same missing pattern (missing the same $p_{mis}$ variables), and $N_{obs}$ observations are complete:

$$
E[Q - \ell_{com}] = -N_{mis} \cdot \text{tr}(\Sigma^{*-1} C) \tag{15}
$$

where $C$ has $\Sigma_{mis|obs}$ in the missing-missing block.

Using block matrix properties:
$$
\text{tr}(\Sigma^{-1} C) = \text{tr}\left(\Sigma^{-1}_{mis|obs,\text{block}} \cdot \Sigma_{mis|obs}\right) \tag{16}
$$

where $\Sigma^{-1}_{mis|obs,\text{block}}$ is the missing-missing block of $\Sigma^{-1}$.

By the Schur complement identity:
$$
(\Sigma^{-1})_{mis,mis} = (\Sigma_{mis|obs})^{-1} \tag{17}
$$

Therefore:
$$
\text{tr}(\Sigma^{-1} C) = \text{tr}((\Sigma_{mis|obs})^{-1} \Sigma_{mis|obs}) = \text{tr}(I_{p_{mis}}) = p_{mis} \tag{18}
$$

## Key Result

$$
\boxed{E[Q(\theta^* | \theta^*) - \ell_{com}(\theta^*)] = -N_{mis} \cdot p_{mis}} \tag{19}
$$

For our simulation with missing rate $\gamma$ and monotone pattern, this becomes:
$$
E[\text{Term 1}] \approx -\gamma \cdot N \cdot \bar{p}_{mis} \tag{20}
$$

where $\bar{p}_{mis}$ is the average number of missing variables per incomplete observation.

## Connection to Log-Determinant (User's Conjecture)

The user's formula $-p - \frac{1}{2}\log|\Sigma_{mis|obs}|$ might relate to the **entropy** contribution.

The differential entropy of $Y_{mis} | Y_{obs}$ is:
$$
H(Y_{mis} | Y_{obs}) = \frac{p_{mis}}{2}\log(2\pi e) + \frac{1}{2}\log|\Sigma_{mis|obs}| \tag{21}
$$

However, the Q-function difference (Eq. 19) is determined by the **dimension** $p_{mis}$, not directly by $\log|\Sigma_{mis|obs}|$.

The $\log|\Sigma_{mis|obs}|$ term would appear if we were computing:
$$
\log f(Y_{mis} | Y_{obs}, \theta) = -\frac{p_{mis}}{2}\log(2\pi) - \frac{1}{2}\log|\Sigma_{mis|obs}| - \frac{1}{2}(Y_{mis} - \tilde{\mu})^\top \Sigma_{mis|obs}^{-1} (Y_{mis} - \tilde{\mu}) \tag{22}
$$

**Question for verification**: Is the user's formula perhaps for a different quantity, such as the observed-data log-likelihood contribution from missing data, rather than the Q-function?

## Implications for Term 1

If Eq. (19) is correct, then:
- Term 1 should scale with **number of missing values**, not with tr(RIV)
- This explains why Term 1 ≈ 0 when averaged across replications (after centering)
- The scaling is $O(N)$, not $O(1)$ like tr(RIV)

This suggests Term 1 may have a deterministic component that cancels out, leaving only sampling variability around zero mean.
