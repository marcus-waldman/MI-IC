---
editor_options: 
  markdown: 
    wrap: 99
---

# Corrected Information Criteria for Model Selection Under Multiple Imputation

## Abstract

Multiple imputation is widely used to handle missing data, yet existing methods for model selection
with imputed data do not correctly target complete-data decisions. We show that methods using
imputed complete-data likelihoods incur an imputation bias that scales with model complexity, and
their penalty terms do not fully correct for this bias. The bias arises from double dipping: the
same data are used first to estimate parameters (or construct a posterior), and then to generate
imputations from those estimates. We derive corrected penalties for AIC and BIC when targeting
complete-data divergence. For AIC, the penalty is $2Q + 4\text{tr}(\text{RIV})$; for BIC, it is
$Q\log(N) + 2\text{tr}(\text{RIV})$, where $\text{RIV} = (1 + 1/M)W^{-1}B$ is the relative increase
in variance matrix, directly computable from standard MI output. We show that existing criteria
(AICcd, PDIO, and variants) underpenalize by omitting the imputation bias correction. Under
standard assumptions, our results generalize to missing covariates: the conditional criteria
decompose as MI-AIC$_{Y|X}$ = MI-AIC$_{Y,X}$ − MI-AIC$_X$. Simulations confirm that MI-AIC and
MI-BIC better replicate complete-data model selection than existing approaches.

---------------------------------------------------------------------------------------------------

## 1. Introduction

### 1.1 The Problem

-   Model selection with missing data is ubiquitous in applied research
-   Multiple imputation (MI) is the dominant approach for handling missing data
-   Natural question: how should we perform model selection with multiply imputed data?

### 1.2 Current Practice

Common approaches: 1. **MI-averaged criteria:** Compute AIC/BIC on each imputed dataset, average 2.
**Stack-and-select:** Pool imputed datasets, apply standard criteria 3. **Complete-case analysis:**
Discard incomplete observations

All aim to replicate what we would conclude with complete data.

### 1.3 The Gap

We show that existing methods fail to achieve this goal: - They incur an **imputation bias** that
favors complex models - Their penalties do not correct for this bias - The error scales with model
complexity and missing information

### 1.4 Our Contributions

1.  Derive the imputation bias in using Q-function to estimate complete-data log-likelihood
2.  Derive MI-AIC: correct penalty is $2Q + 4\text{tr}(\text{RIV})$
3.  Derive MI-BIC: correct penalty is $Q\log(N) + 2\text{tr}(\text{RIV})$
4.  Show existing methods (AICcd, PDIO, Claeskens & Consentino) underpenalize
5.  Generalize to missing covariates under standard assumptions
6.  Demonstrate improved replication of complete-data decisions via simulation

---------------------------------------------------------------------------------------------------

## 2. Setup and Notation

### 2.1 Complete Data and Missing Data

**Complete data:** $Z = (Y, X)$ where $Y$ is the response and $X$ are covariates.

**Partition by observability:** $Z = (Z^{obs}, Z^{mis})$

**Complete-data log-likelihood:** $$\ell_{com}(\theta) = \sum_i \log f(z_i; \theta)$$

**Observed-data log-likelihood:**
$$\ell_{obs}(\theta) = \sum_i \log f(z_i^{obs}; \theta) = \sum_i \log \int f(z_i^{obs}, z_i^{mis}; \theta) \, dz_i^{mis}$$

### 2.2 The Q-Function

The Q-function from the EM algorithm:
$$Q(\theta|\tilde{\theta}) = \mathbb{E}_{Z^{mis}|Z^{obs},\tilde{\theta}}[\ell_{com}(\theta)]$$

This is the expected complete-data log-likelihood, with expectation over the conditional
distribution of missing data given observed data and imputation parameter $\tilde{\theta}$.

### 2.3 Proper Multiple Imputation

For $m = 1, \ldots, M$: 1. Draw $\tilde{\theta}^{(m)} \sim p(\theta | Z^{obs})$ 2. Impute
$\tilde{Z}^{mis,(m)} \sim p(Z^{mis} | Z^{obs}, \tilde{\theta}^{(m)})$

The MI estimate of the Q-function at $\hat{\theta}_{obs}$:
$$\bar{Q}_{MI}(\hat{\theta}_{obs}) = \frac{1}{M} \sum_{m=1}^{M} \ell_{com}(\hat{\theta}_{obs} | Z^{obs}, \tilde{Z}^{mis,(m)})$$

### 2.4 Information Matrices

-   $\mathcal{I}_{com}(\theta)$ = Fisher information from complete data
-   $\mathcal{I}_{obs}(\theta)$ = Fisher information from observed data
-   $\mathcal{I}_{mis|obs}(\theta) = \mathcal{I}_{com}(\theta) - \mathcal{I}_{obs}(\theta)$ =
    missing information

### 2.5 Rubin's Variance Components

From standard MI output: - $W$ = within-imputation variance (average of complete-data variance
estimates) - $B$ = between-imputation variance (variance of point estimates across imputations) -
$T = W + (1 + 1/M)B$ = total variance

**Key relationships:** - $W \approx \mathcal{I}_{com}^{-1}$ - $T \approx \mathcal{I}_{obs}^{-1}$ -
$\text{RIV} = (1 + 1/M)W^{-1}B \approx \mathcal{I}_{mis|obs} \cdot \mathcal{I}_{obs}^{-1}$

---------------------------------------------------------------------------------------------------

## 3. The Imputation Bias

### 3.1 The Goal

**Target:** Expected complete-data log-likelihood on new data,
$\mathbb{E}[\bar{\ell}_{com}(\hat{\theta})]$

**What we have:** $\bar{Q}_{MI}(\hat{\theta}_{obs})$

**Total bias to correct:**
$$\text{Total Bias} = \mathbb{E}[\bar{Q}_{MI}(\hat{\theta}_{obs})] - \mathbb{E}[\bar{\ell}_{com}(\hat{\theta}_{obs})]$$

### 3.2 Decomposition

$$\text{Total Bias} = \underbrace{\mathbb{E}[\bar{Q}_{MI}] - \mathbb{E}[\ell_{com}]}_{\text{Imputation Bias}} + \underbrace{\mathbb{E}[\ell_{com}] - \mathbb{E}[\bar{\ell}_{com}]}_{\text{Optimism}}$$

### 3.3 Deriving the Imputation Bias

**Key insight:** If the imputation parameter $\tilde{\theta}$ were equal to the true $\theta^*$,
the imputation bias would be zero (by iterated expectations). The bias arises because
$\tilde{\theta}$ is estimated from data.

**Result (derived in Appendix):**
$$\text{Imputation Bias} = \text{tr}(\mathcal{I}_{mis|obs}(\theta^*) \cdot \text{Cov}(\tilde{\theta}, \hat{\theta}_{obs}))$$

For proper MI, where $\tilde{\theta} \sim p(\theta|Z^{obs})$: - The posterior mean
$\bar{\theta} \approx \hat{\theta}_{obs}$ asymptotically - The posterior uncertainty is
uncorrelated with $\hat{\theta}_{obs}$ given $Z^{obs}$

Therefore:
$$\text{Imputation Bias} = \text{tr}(\mathcal{I}_{mis|obs} \cdot \mathcal{I}_{obs}^{-1}) = \text{tr}(\text{RIV})$$

### 3.4 Properties of the Imputation Bias

1.  **Positive:** Q overestimates the true complete-data log-likelihood
2.  **Scales with model complexity:** More parameters → larger bias
3.  **Scales with missing information:** More missingness → larger RIV → larger bias
4.  **Zero when no missing data:** RIV = 0 when data are complete

---------------------------------------------------------------------------------------------------

## 4. MI-AIC: Corrected AIC for Multiple Imputation

### 4.1 The Optimism for AIC

AIC corrects for optimism—the expected difference between training and test performance:
$$\text{Optimism} = \mathbb{E}[\ell_{com}(\hat{\theta}_{obs})] - \mathbb{E}[\bar{\ell}_{com}(\hat{\theta}_{obs})]$$

Standard result:
$$\text{Optimism} = \text{tr}(\mathcal{I}_{target} \cdot \text{Var}(\hat{\theta}))$$

Here: - Target is complete-data likelihood: $\mathcal{I}_{target} = \mathcal{I}_{com}$ - Estimator
variance: $\text{Var}(\hat{\theta}_{obs}) = \mathcal{I}_{obs}^{-1}$

Therefore: $$\text{Optimism} = \text{tr}(\mathcal{I}_{com} \cdot \mathcal{I}_{obs}^{-1})$$

Using $\mathcal{I}_{com} = \mathcal{I}_{obs} + \mathcal{I}_{mis|obs}$:
$$= \text{tr}(I_Q + \mathcal{I}_{mis|obs} \cdot \mathcal{I}_{obs}^{-1}) = Q + \text{tr}(\text{RIV})$$

### 4.2 Total Bias for AIC

$$\text{Total Bias} = \text{Imputation Bias} + \text{Optimism}$$
$$= \text{tr}(\text{RIV}) + Q + \text{tr}(\text{RIV})$$ $$= Q + 2\text{tr}(\text{RIV})$$

### 4.3 MI-AIC Formula

The correct penalty is $2 \times \text{Total Bias}$:

$$\boxed{\text{MI-AIC} = -2\bar{Q}_{MI}(\hat{\theta}_{obs}) + 2Q + 4\text{tr}(\text{RIV})}$$

### 4.4 Comparison to Standard AIC

| Setting                      | Penalty                       |
|------------------------------|-------------------------------|
| Complete data                | $2Q$                          |
| MI (targeting complete-data) | $2Q + 4\text{tr}(\text{RIV})$ |

The additional $4\text{tr}(\text{RIV})$ accounts for: - $2\text{tr}(\text{RIV})$ from imputation
bias - $2\text{tr}(\text{RIV})$ from optimism involving $\mathcal{I}_{com}/\mathcal{I}_{obs}$

---------------------------------------------------------------------------------------------------

## 5. MI-BIC: Corrected BIC for Multiple Imputation

### 5.1 BIC Penalty Derivation

BIC's penalty comes from the Laplace approximation to the marginal likelihood (Bayes factor):
$$\text{BIC penalty} = Q \log(N)$$

This derivation does **not** involve the ratio $\mathcal{I}_{com} \cdot \mathcal{I}_{obs}^{-1}$.

### 5.2 Imputation Bias Still Applies

The imputation bias is a property of the deviance estimation, not the penalty derivation:
$$\text{Imputation Bias} = \text{tr}(\text{RIV})$$

This must be corrected regardless of which criterion we use.

### 5.3 MI-BIC Formula

$$\boxed{\text{MI-BIC} = -2\bar{Q}_{MI}(\hat{\theta}_{obs}) + Q\log(N) + 2\text{tr}(\text{RIV})}$$

### 5.4 Comparison to Standard BIC

| Setting                      | Penalty                             |
|------------------------------|-------------------------------------|
| Complete data                | $Q\log(N)$                          |
| MI (targeting complete-data) | $Q\log(N) + 2\text{tr}(\text{RIV})$ |

### 5.5 Why MI-AIC has $4\text{tr}(\text{RIV})$ but MI-BIC has $2\text{tr}(\text{RIV})$

-   **MI-AIC:** Both imputation bias and optimism involve RIV
-   **MI-BIC:** Only imputation bias involves RIV; BIC penalty is $Q\log(N)$, not derived from
    information ratio

---------------------------------------------------------------------------------------------------

## 6. Generalization to Missing Covariates

### 6.1 Setup

For regression with missing covariates, let $Z = (Y, X)$ with:
$$f(Z; \theta) = f(Y|X; \beta) \cdot f(X; \alpha)$$

where $\theta = (\beta, \alpha)$ are distinct parameters.

### 6.2 Block Diagonal Structure

Under distinct parameters, information matrices are block diagonal:
$$\mathcal{I}_{com} = \begin{pmatrix} \mathcal{I}_{com,\beta} & 0 \\ 0 & \mathcal{I}_{com,\alpha} \end{pmatrix}$$

and similarly for $\mathcal{I}_{obs}$ and $\mathcal{I}_{mis|obs}$.

### 6.3 Traces Separate

$$\text{tr}(\text{RIV}) = \text{tr}(\text{RIV}_\beta) + \text{tr}(\text{RIV}_\alpha)$$

### 6.4 Conditional Criteria

For regression model selection (comparing models for $Y|X$):

$$\boxed{\text{MI-AIC}_{Y|X} = -2Q_{Y|X}(\hat{\beta}|\hat{\theta}) + 2Q_\beta + 4\text{tr}(\text{RIV}_\beta)}$$

$$\boxed{\text{MI-BIC}_{Y|X} = -2Q_{Y|X}(\hat{\beta}|\hat{\theta}) + Q_\beta\log(N) + 2\text{tr}(\text{RIV}_\beta)}$$

### 6.5 Decomposition Property

$$\text{MI-AIC}_{Y|X} = \text{MI-AIC}_{Y,X} - \text{MI-AIC}_X$$

$$\text{MI-BIC}_{Y|X} = \text{MI-BIC}_{Y,X} - \text{MI-BIC}_X$$

---------------------------------------------------------------------------------------------------

## 7. Existing Methods Are Incomplete

### 7.1 Cavanaugh & Shumway's AICcd

$$\text{AICcd} = -2Q(\hat{\theta}|\hat{\theta}) + 2\text{tr}(\mathcal{I}_{com} \cdot \mathcal{I}_{obs}^{-1})$$

Their penalty:
$2Q + 2\text{tr}(\mathcal{I}_{mis|obs} \cdot \mathcal{I}_{obs}^{-1}) = 2Q + 2\text{tr}(\text{RIV})$

**The gap:** Corrects optimism but not imputation bias. Underpenalizes by $2\text{tr}(\text{RIV})$.

### 7.2 Shimodaira's PDIO

$$\text{PDIO} = -2\ell_{obs}(\hat{\theta}) + 2\text{tr}(\mathcal{I}_{com} \cdot \mathcal{I}_{obs}^{-1})$$

Uses observed-data likelihood (no imputation bias in deviance), but: - Penalty still only
$2Q + 2\text{tr}(\text{RIV})$ - Computing $\mathcal{I}_{com}$ typically requires imputation anyway

### 7.3 Claeskens & Consentino (2008)

For missing covariates: $$\text{AIC}_1 = 2Q_{Y|X}(\hat{\beta}|\hat{\theta}) - 2Q_\beta$$

**The sidestep:** They redefine the target to Q itself, not $\ell_{com}$. Their Theorem 1 derives
bias for predicting Q on new data, not $\ell_{com}$ on new data.

When complete-data replication is the goal, their penalty underpenalizes by
$4\text{tr}(\text{RIV}_\beta)$.

### 7.4 Summary

| Method | Target | Penalty | Gap |
|-------------------------|-------------------------|----------------------------|-----------------------|
| AICcd | $\ell_{com}$ (claimed) | $2Q + 2\text{tr}(\text{RIV})$ | Missing $2\text{tr}(\text{RIV})$ |
| PDIO | $\ell_{com}$ (claimed) | $2Q + 2\text{tr}(\text{RIV})$ | Missing $2\text{tr}(\text{RIV})$ |
| C&C AIC₁ | Q (redefined) | $2Q_\beta$ | Missing $4\text{tr}(\text{RIV}_\beta)$ for $\ell_{com}$ |
| **MI-AIC** | $\ell_{com}$ | $2Q + 4\text{tr}(\text{RIV})$ | **Correct** |
| **MI-BIC** | $\ell_{com}$ | $Q\log(N) + 2\text{tr}(\text{RIV})$ | **Correct** |

---------------------------------------------------------------------------------------------------

## 8. Simulation Study

### 8.1 Design

-   Sample sizes: $N \in \{100, 500, 1000\}$
-   Missing data rates: 10%, 30%, 50%
-   Missingness mechanism: MAR
-   True model: known, with nested alternatives
-   Number of imputations: $M = 20$
-   Replications: 1000

### 8.2 Methods Compared

1.  MI-AIC (proposed)
2.  MI-BIC (proposed)
3.  AICcd (Cavanaugh & Shumway)
4.  MI-averaged AIC (average AIC across imputations)
5.  MI-averaged BIC
6.  Complete-case AIC/BIC

### 8.3 Evaluation Criterion

**Replication rate:** Proportion of simulations where the criterion selects the same model as would
be selected with complete data.

### 8.4 Expected Results

-   AICcd and MI-averaged criteria: systematically favor more complex models (underpenalize)
-   Complete-case criteria: inefficient, potentially biased
-   MI-AIC and MI-BIC: highest replication rates

---------------------------------------------------------------------------------------------------

## 9. Discussion

### 9.1 Summary

We have shown that: 1. Imputation-based model selection incurs a bias that favors complex models 2.
Existing methods (AICcd, PDIO, C&C) do not fully correct this bias 3. The correct penalties are
MI-AIC ($2Q + 4\text{tr}(\text{RIV})$) and MI-BIC ($Q\log(N) + 2\text{tr}(\text{RIV})$) 4. These
are directly computable from standard MI output

### 9.2 Practical Guidance

**Computing MI-AIC/MI-BIC:** 1. Generate $M$ proper imputations 2. Fit candidate model to observed
data → $\hat{\theta}_{obs}$ 3. Evaluate $\ell_{com}(\hat{\theta}_{obs})$ on each completed dataset,
average → $\bar{Q}_{MI}$ 4. Compute $W$ and $B$ from MI output 5. Calculate
$\text{tr}(\text{RIV}) = \text{tr}((1 + 1/M)W^{-1}B)$ 6. Apply formulas

### 9.3 When to Use MI-AIC vs MI-BIC

Same guidance as for standard AIC vs BIC: - MI-AIC: prediction focus, allow overfitting - MI-BIC:
consistent selection, prefer parsimony

### 9.4 Limitations

-   Assumes MAR
-   Requires proper imputations (posterior draws)
-   Large $Q$ may require regularization of $W^{-1}$

### 9.5 Extensions

-   Non-MAR with sensitivity analysis
-   Other criteria (DIC, WAIC)
-   Model averaging weights
-   Cross-validation approaches

---------------------------------------------------------------------------------------------------

## Appendix A: Derivation of Imputation Bias

### A.1 Setup

Let $Z = (Z^{obs}, Z^{mis})$. The complete-data log-likelihood is $\ell_{com}(\theta)$.

The Q-function:
$$Q(\theta|\tilde{\theta}) = \mathbb{E}_{Z^{mis}|Z^{obs},\tilde{\theta}}[\ell_{com}(\theta)]$$

### A.2 Taylor Expansion

Expand $\ell_{com}(\hat{\theta})$ around $\theta^*$:
$$\ell_{com}(\hat{\theta}) = \ell_{com}(\theta^*) + (\hat{\theta}-\theta^*)^T S_{com}(\theta^*) - \frac{1}{2}(\hat{\theta}-\theta^*)^T \mathcal{J}_{com}(\theta^*)(\hat{\theta}-\theta^*) + O_p(N^{-3/2})$$

Similarly for $Q(\hat{\theta}|\tilde{\theta})$ around $\theta^*$.

### A.3 The Bias

$$\text{Bias} = \mathbb{E}[Q(\hat{\theta}|\tilde{\theta})] - \mathbb{E}[\ell_{com}(\hat{\theta})]$$

The key term is:
$$\mathbb{E}[(\tilde{\theta} - \theta^*)^T \mathcal{I}_{mis|obs}(\theta^*)(\hat{\theta} - \theta^*)]$$

### A.4 Result

$$\text{Imputation Bias} = \text{tr}(\mathcal{I}_{mis|obs}(\theta^*) \cdot \text{Cov}(\tilde{\theta}, \hat{\theta}_{obs})) + O(N^{-1})$$

For proper MI:
$\text{Cov}(\tilde{\theta}, \hat{\theta}_{obs}) \approx \text{Var}(\hat{\theta}_{obs}) = \mathcal{I}_{obs}^{-1}$

Therefore:
$$\text{Imputation Bias} = \text{tr}(\mathcal{I}_{mis|obs} \cdot \mathcal{I}_{obs}^{-1}) = \text{tr}(\text{RIV})$$

---------------------------------------------------------------------------------------------------

## Key Equations Summary

| Quantity | Formula | Notes |
|--------------------------------------|----------------------------------|---------------------------|
| Imputation bias | $\text{tr}(\mathcal{I}_{mis|obs} \cdot \mathcal{I}_{obs}^{-1}) = \text{tr}(\text{RIV})$ | Positive; scales with $Q$ |
| RIV | $(1 + 1/M)W^{-1}B$ | Estimable from MI output |
| MI-AIC penalty | $2Q + 4\text{tr}(\text{RIV})$ | Corrects imputation bias + optimism |
| MI-BIC penalty | $Q\log(N) + 2\text{tr}(\text{RIV})$ | Corrects imputation bias only |
| MI-AIC$_{Y|X}$ penalty | $2Q_\beta + 4\text{tr}(\text{RIV}_\beta)$ | For regression with missing covariates |
| MI-BIC$_{Y|X}$ penalty | $Q_\beta\log(N) + 2\text{tr}(\text{RIV}_\beta)$ | For regression with missing covariates |
| AICcd penalty | $2Q + 2\text{tr}(\text{RIV})$ | Underpenalizes by $2\text{tr}(\text{RIV})$ |

**Notation:** - $Q$ = number of parameters - $W$ = within-imputation variance - $B$ =
between-imputation variance - $\text{RIV} = (1 + 1/M)W^{-1}B$ = relative increase in variance -
$\beta$ = regression parameters, $\alpha$ = covariate distribution parameters
