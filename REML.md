# PAGssGBLUP: Restricted maximum likelihood and AI-REML

Part of the methodology of [PAGssGBLUP](../README.md). Symbols are defined in [Notation](../README.md#notation), and references are listed in [References](../README.md#references).

← [Sections 8–9](genomic-information-and-mme.md) · [References](../README.md#references) →

---

## 10. Restricted maximum likelihood

The MME need $`\alpha = \sigma^2_e/\sigma^2_a`$, so the variance components $`\boldsymbol{\theta} = (\sigma^2_a, \sigma^2_e)^\top`$ must be estimated first. This is done by restricted (or residual) maximum likelihood, REML (Patterson and Thompson, 1971).

### 10.1 Why REML and not ML

Under the model, $`\mathbf{y} \sim N(\mathbf{Xb}, \mathbf{V})`$, with

```math
\mathbf{V} = \text{Var}[\mathbf{y}] = \text{Var}[\mathbf{Wa} + \mathbf{e}] = \mathbf{W}\mathbf{H}\mathbf{W}^\top\sigma^2_a + \mathbf{I}\sigma^2_e
```

because $`\mathbf{a}`$ and $`\mathbf{e}`$ are independent. $`\mathbf{W}\mathbf{H}\mathbf{W}^\top`$ is the block of $`\mathbf{H}`$ for the phenotyped animals, in the order of $`\mathbf{y}`$.

Maximum likelihood (ML) maximizes the likelihood of $`\mathbf{y}`$ as if $`\mathbf{b}`$ were known. It ignores that $`\mathbf{b}`$ is estimated from the same data, so the variances are biased downwards. For example, in $`\mathbf{y} = \mathbf{Xb} + \mathbf{e}`$, ML gives $`\hat{\sigma}^2_e = \text{SSE}/N`$ instead of $`\text{SSE}/(N-p)`$.

REML uses only the part of the data that is free of $`\mathbf{b}`$. Let $`\mathbf{S}`$ be any $`N \times (N-p)`$ matrix of full column rank with $`\mathbf{S}^\top\mathbf{X} = \mathbf{0}`$ (for example, $`N - p`$ linearly independent columns of $`\mathbf{I} - \mathbf{X}(\mathbf{X}^\top\mathbf{X})^{-1}\mathbf{X}^\top`$). Each element of $`\mathbf{S}^\top\mathbf{y}`$ is a linear combination of the records in which the fixed effects cancel out (an *error contrast*):

```math
\mathbf{S}^\top\mathbf{y} = \mathbf{S}^\top\mathbf{Xb} + \mathbf{S}^\top(\mathbf{Wa} + \mathbf{e}) = \mathbf{S}^\top(\mathbf{Wa} + \mathbf{e}) \sim N\left(\mathbf{0}, \mathbf{S}^\top\mathbf{V}\mathbf{S}\right)
```

REML is ML applied to $`\mathbf{S}^\top\mathbf{y}`$ instead of $`\mathbf{y}`$. Its distribution does not depend on $`\mathbf{b}`$, and the $`p`$ degrees of freedom used by the fixed effects are accounted for ($`N - p`$ contrasts instead of $`N`$ records). The result does not depend on the choice of $`\mathbf{S}`$, and $`\mathbf{S}`$ is never computed: the REML log-likelihood is written with $`\mathbf{V}`$ and $`\mathbf{X}`$ only.

### 10.2 REML log-likelihood

```math
L(\boldsymbol{\theta}) = -\frac{1}{2}\left[\log\left|\mathbf{V}\right| + \log\left|\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\right| + \mathbf{y}^\top\mathbf{P}\mathbf{y}\right] + \text{constant}
```

```math
\mathbf{P} = \mathbf{V}^{-1} - \mathbf{V}^{-1}\mathbf{X}\left(\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\right)^{-1}\mathbf{X}^\top\mathbf{V}^{-1}
```

$`\mathbf{P}\mathbf{X} = \mathbf{0}`$, so $`\mathbf{P}`$ removes the fixed effects. Also, $`\mathbf{P}\mathbf{y} = \mathbf{V}^{-1}(\mathbf{y} - \mathbf{X}\hat{\mathbf{b}})`$, where $`\hat{\mathbf{b}}`$ is the generalized least squares solution, $`\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\hat{\mathbf{b}} = \mathbf{X}^\top\mathbf{V}^{-1}\mathbf{y}`$. The three terms of $`L`$ are:

1. $`\mathbf{y}^\top\mathbf{P}\mathbf{y} = (\mathbf{y} - \mathbf{X}\hat{\mathbf{b}})^\top\mathbf{V}^{-1}(\mathbf{y} - \mathbf{X}\hat{\mathbf{b}})`$: a weighted sum of squares of residuals.
2. $`\log\lvert\mathbf{V}\rvert`$: depends only on the variance structure.
3. $`\log\lvert\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\rvert`$: a penalty for estimating $`\mathbf{b}`$. ML does not have this term.

### 10.3 Maximizing L

$`L`$ has no closed-form maximum, so $`\boldsymbol{\theta}`$ is found iteratively from the derivatives of $`L`$. $`\mathbf{V}`$ is linear in $`\boldsymbol{\theta}`$:

```math
\mathbf{V}_a = \frac{\partial\mathbf{V}}{\partial\sigma^2_a} = \mathbf{W}\mathbf{H}\mathbf{W}^\top,
\qquad
\mathbf{V}_e = \frac{\partial\mathbf{V}}{\partial\sigma^2_e} = \mathbf{I}
```

**Score** (first derivatives), for $`i \in \lbrace a, e \rbrace`$:

```math
\frac{\partial L}{\partial\theta_i} = -\frac{1}{2}\left[\text{tr}\left(\mathbf{P}\mathbf{V}_i\right) - \mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{y}\right]
```

At the maximum the score is zero: the quadratic form of the data, $`\mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{y}`$, equals its expected value, $`\text{tr}(\mathbf{P}\mathbf{V}_i)`$.

**Information matrices** ($`2 \times 2`$), element $`(i, j)`$:

| Information | Element $`(i, j)`$ | Method |
|:---|:---|:---|
| Observed, $`-\partial^2 L/\partial\theta_i\partial\theta_j`$ | $`\mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j\mathbf{P}\mathbf{y} - \frac{1}{2}\text{tr}(\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j)`$ | Newton-Raphson |
| Expected (Fisher) | $`\frac{1}{2}\text{tr}(\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j)`$ | Fisher scoring |
| Average (AI) | $`\frac{1}{2}\mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j\mathbf{P}\mathbf{y}`$ | AI-REML (Gilmour et al., 1995) |

The AI matrix is the average of the observed and the expected information. It has no trace terms, so it is cheap to compute. Each iteration updates

```math
\boldsymbol{\theta}^{(t+1)} = \boldsymbol{\theta}^{(t)} + \left(\mathbf{AI}^{(t)}\right)^{-1}\frac{\partial L}{\partial\boldsymbol{\theta}}\bigg|_{\boldsymbol{\theta}^{(t)}}
```

At convergence, $`\text{Var}[\hat{\boldsymbol{\theta}}] \approx \mathbf{AI}^{-1}`$, which gives the standard errors.

### 10.4 Computing L and its derivatives from the MME

When $`N`$ is large, $`\mathbf{V}`$ ($`N \times N`$) is too large to invert, and everything is obtained from the MME. Let $`\mathbf{T} = \begin{pmatrix}\mathbf{X} & \mathbf{W}\end{pmatrix}`$ and

```math
\mathbf{C} =
\begin{pmatrix}
\mathbf{X}^\top\mathbf{X} & \mathbf{X}^\top\mathbf{W} \\
\mathbf{W}^\top\mathbf{X} & \mathbf{W}^\top\mathbf{W} + \mathbf{H}^{-1}\alpha
\end{pmatrix}
```

The key identity is

```math
\mathbf{P} = \frac{1}{\sigma^2_e}\left(\mathbf{I} - \mathbf{T}\mathbf{C}^{-1}\mathbf{T}^\top\right)
```

With the MME solutions $`\hat{\mathbf{b}}`$, $`\hat{\mathbf{a}}`$ and the residuals $`\hat{\mathbf{e}} = \mathbf{y} - \mathbf{X}\hat{\mathbf{b}} - \mathbf{W}\hat{\mathbf{a}}`$:

```math
\mathbf{P}\mathbf{y} = \frac{\hat{\mathbf{e}}}{\sigma^2_e},
\qquad
\mathbf{y}^\top\mathbf{P}\mathbf{y} = \frac{\mathbf{y}^\top\hat{\mathbf{e}}}{\sigma^2_e}
```

**Log-likelihood** (e.g., Meyer, 1989), with $`p`$ the number of columns of $`\mathbf{X}`$ (full rank). $`\log\lvert\mathbf{H}\rvert`$ does not depend on $`\boldsymbol{\theta}`$; it is written explicitly so that this $`L`$ equals the one in 10.2, not only up to a constant:

```math
-2L = (N - p - K)\log\sigma^2_e + K\log\sigma^2_a + \log\left|\mathbf{H}\right| + \log\left|\mathbf{C}\right| + \frac{\mathbf{y}^\top\hat{\mathbf{e}}}{\sigma^2_e} + \text{constant}
```

**Score.** $`\mathbf{C}^{22}`$ is the block of $`\mathbf{C}^{-1}`$ for $`\mathbf{a}`$, with $`\text{Var}[\mathbf{a} - \hat{\mathbf{a}}] = \mathbf{C}^{22}\sigma^2_e`$ (prediction error variance, PEV):

```math
\frac{\partial L}{\partial\sigma^2_a} =
\frac{1}{2}\left[\frac{\hat{\mathbf{a}}^\top\mathbf{H}^{-1}\hat{\mathbf{a}}}{\sigma^4_a} - \frac{K}{\sigma^2_a} + \frac{\text{tr}\left(\mathbf{H}^{-1}\mathbf{C}^{22}\right)\sigma^2_e}{\sigma^4_a}\right]
```

```math
\frac{\partial L}{\partial\sigma^2_e} =
\frac{1}{2}\left[\frac{\hat{\mathbf{e}}^\top\hat{\mathbf{e}}}{\sigma^4_e} - \frac{N - p - K}{\sigma^2_e} - \frac{\text{tr}\left(\mathbf{H}^{-1}\mathbf{C}^{22}\right)}{\sigma^2_a}\right]
```

**AI matrix.** The vectors $`\mathbf{V}_i\mathbf{P}\mathbf{y}`$ are

```math
\mathbf{V}_a\mathbf{P}\mathbf{y} = \frac{\mathbf{W}\hat{\mathbf{a}}}{\sigma^2_a},
\qquad
\mathbf{V}_e\mathbf{P}\mathbf{y} = \frac{\hat{\mathbf{e}}}{\sigma^2_e}
```

With $`\mathbf{F} = \begin{pmatrix}\mathbf{W}\hat{\mathbf{a}}/\sigma^2_a & \hat{\mathbf{e}}/\sigma^2_e\end{pmatrix} \in \mathbb{R}^{N \times 2}`$, $`\mathbf{AI} = \frac{1}{2}\mathbf{F}^\top\mathbf{P}\mathbf{F}`$, and by the key identity

```math
\mathbf{AI} = \frac{1}{2\sigma^2_e}\left(\mathbf{F}^\top\mathbf{F} - \mathbf{F}^\top\mathbf{T}\mathbf{C}^{-1}\mathbf{T}^\top\mathbf{F}\right)
```

That is, the MME are solved with the two columns of $`\mathbf{F}`$ in place of $`\mathbf{y}`$ (right-hand side $`\mathbf{T}^\top\mathbf{F}`$).

### 10.5 EM-REML

Setting both scores to zero and rearranging gives (Mrode, 2014, Eqns 17.5 and 17.6)

```math
\sigma^2_e = \frac{\mathbf{y}^\top\hat{\mathbf{e}}}{N - p},
\qquad
\sigma^2_a = \frac{\hat{\mathbf{a}}^\top\mathbf{H}^{-1}\hat{\mathbf{a}} + \text{tr}\left(\mathbf{H}^{-1}\mathbf{C}^{22}\right)\sigma^2_e}{K}
```

Used as an update (right-hand sides computed with the current $`\boldsymbol{\theta}`$), this is the EM algorithm (Dempster et al., 1977). $`\sigma^2_e`$ comes from the residuals, and $`\sigma^2_a`$ from the predicted breeding values plus their PEV. EM estimates always stay positive, but convergence is slow. AI-REML converges in a few iterations, but an update can give a negative variance. A common strategy is to start with a few EM iterations, or to use an EM step when the AI step leaves the parameter space.

### 10.6 AI-REML algorithm

1. **Start values**, e.g., from a guess of heritability $`h^2_0`$ and the phenotypic variance $`\sigma^2_y`$: $`\sigma^2_a = h^2_0\sigma^2_y`$ and $`\sigma^2_e = (1 - h^2_0)\sigma^2_y`$.
2. **Compute what $`L`$, the score and AI need.** Either $`\mathbf{V}`$ and $`\mathbf{P}`$ (small $`N`$), or the MME solutions $`\hat{\mathbf{b}}`$, $`\hat{\mathbf{a}}`$, $`\hat{\mathbf{e}}`$ and $`\mathbf{C}^{22}`$ with $`\alpha = \sigma^2_e/\sigma^2_a`$ (large $`N`$).
3. **Compute** $`L`$, the score and the AI matrix.
4. **Update** $`\boldsymbol{\theta}`$. A variance cannot be negative. If the update makes it negative, it is set to a small positive value (the boundary), and it stays there while its score is negative. An EM step is an alternative.
5. **Repeat** steps 2–4 until convergence (changes in $`L`$ and $`\boldsymbol{\theta}`$ below a small tolerance).
6. **Report** $`\hat{\sigma}^2_a`$, $`\hat{\sigma}^2_e`$, their standard errors from $`\mathbf{AI}^{-1}`$, and $`\hat{h}^2 = \hat{\sigma}^2_a/(\hat{\sigma}^2_a + \hat{\sigma}^2_e)`$ with its standard error by the delta method.

**Note on Mrode (2014, Chapter 17).** Two formulas there seem to have misprints. $`\mathbf{P}`$ is printed without the first $`\mathbf{V}^{-1}`$ in the second term, and in Eqns 17.2 and 17.3 some signs and powers of the variances differ from the ones above. The expressions here are consistent with Eqns 17.5 and 17.6 and were checked numerically (Section 11.4).

---

## 11. AI-REML implementation

### 11.1 Which route

Here there are far fewer records $`N`$ than MME equations $`K + p`$ (the pedigree has many more animals than phenotypes). Factorizing $`\mathbf{V}`$ ($`N \times N`$) is then much cheaper than factorizing $`\mathbf{C}`$ ($`(K+p) \times (K+p)`$). $`L`$, the score and the AI matrix are therefore computed directly with the formulas of [10.3](#103-maximizing-l), one line of code per formula. The MME route of [10.4](#104-computing-l-and-its-derivatives-from-the-mme) gives the same numbers and is used as a check. With many records, the MME route with sparse matrices is the one to use.

REML needs only $`\mathbf{y}`$, $`\mathbf{X}_{\text{ord}}`$ and $`\mathbf{W}\mathbf{H}\mathbf{W}^\top`$: $`\mathbf{V}`$ depends on $`\mathbf{H}`$ only through the phenotyped animals. In the code, `Hyy = H[y_ids, y_ids]` is $`\mathbf{W}\mathbf{H}\mathbf{W}^\top`$, in the order of $`\mathbf{y}`$.

### 11.2 Steps 2 and 3: `reml_terms()`

For a given $`\boldsymbol{\theta}`$:

| Quantity | Formula | Computation |
|:---|:---|:---|
| $`\mathbf{V}`$ | $`\mathbf{W}\mathbf{H}\mathbf{W}^\top\sigma^2_a + \mathbf{I}\sigma^2_e`$ | `sa * Hyy + se * diag(N)` |
| $`\mathbf{V}^{-1}`$, $`\log\lvert\mathbf{V}\rvert`$ | Cholesky $`\mathbf{V} = \mathbf{U}^\top\mathbf{U}`$; $`\log\lvert\mathbf{V}\rvert = 2\sum_i\log u_{ii}`$ | `chol()`, `chol2inv()` |
| $`\mathbf{P}`$ | $`\mathbf{V}^{-1} - \mathbf{V}^{-1}\mathbf{X}(\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}^\top\mathbf{V}^{-1}`$ | `solve(XViX, t(ViX))` |
| $`L`$ | $`-\frac{1}{2}\left[\log\lvert\mathbf{V}\rvert + \log\lvert\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\rvert + \mathbf{y}^\top\mathbf{P}\mathbf{y}\right]`$ | `determinant()` for the second term |
| $`\mathbf{V}_a\mathbf{P}\mathbf{y}`$, $`\mathbf{V}_e\mathbf{P}\mathbf{y}`$ | $`\mathbf{W}\mathbf{H}\mathbf{W}^\top\mathbf{P}\mathbf{y}`$, $`\mathbf{P}\mathbf{y}`$ | |
| Score | $`-\frac{1}{2}\left[\text{tr}(\mathbf{P}\mathbf{V}_i) - \mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{y}\right]`$ | $`\text{tr}(\mathbf{P}\mathbf{V}_a) = \sum_{ij}p_{ij}(\mathbf{W}\mathbf{H}\mathbf{W}^\top)_{ij}`$ (both symmetric), avoiding the matrix product |
| AI | $`\frac{1}{2}\mathbf{F}^\top\mathbf{P}\mathbf{F}`$, $`\mathbf{F} = (\mathbf{V}_a\mathbf{P}\mathbf{y}, \mathbf{V}_e\mathbf{P}\mathbf{y})`$ | `crossprod(Fm, P %*% Fm)` |

Here $`p_{ij}`$ is element $`(i, j)`$ of $`\mathbf{P}`$ (not an allele frequency).

### 11.3 Steps 1, 4, 5 and 6: `ai_reml()`

1. **Start values.** $`\sigma^2_y`$ is the residual variance of $`\mathbf{y}`$ after fitting only the fixed effects by ordinary least squares, $`\sigma^2_y = \text{SSE}/(N - p)`$. The start heritability is $`h^2_0 = 0.3`$.
2. **Update.** $`\boldsymbol{\theta} \leftarrow \boldsymbol{\theta} + \mathbf{AI}^{-1}\,\partial L/\partial\boldsymbol{\theta}`$. If the update makes a variance $`\le 0`$, it is set to $`\epsilon = 10^{-8}\sigma^2_y`$ (the boundary). While a variance is at the boundary and its score is negative ($`L`$ would only increase if the variance could be negative), it stays fixed. Only the other variance is updated, with its own element of $`\mathbf{AI}`$ (an active-set update).
3. **Convergence.** Stop when the largest change in $`\boldsymbol{\theta}`$, relative to $`\sigma^2_y`$, and the change in $`L`$ are both below $`10^{-8}`$.
4. **Estimates and standard errors.** $`\text{Var}[\hat{\boldsymbol{\theta}}] \approx \mathbf{AI}^{-1}`$. For $`\hat{h}^2`$, the delta method gives

```math
\text{Var}\left[\hat{h}^2\right] \approx \left(\nabla h^2\right)^\top\mathbf{AI}^{-1}\left(\nabla h^2\right),
\qquad
\nabla h^2 = \frac{\partial h^2}{\partial\boldsymbol{\theta}} = \frac{1}{\left(\sigma^2_a + \sigma^2_e\right)^2}\begin{pmatrix}\sigma^2_e \\ -\sigma^2_a\end{pmatrix}
```

If a variance is at the boundary, the large-sample theory does not apply. Its standard error, and that of $`\hat{h}^2`$, are not reported (`NA`), and the standard error of the other variance uses only its own element of $`\mathbf{AI}`$.

**Profile log-likelihood.** For each $`h^2`$ on a grid, $`L`$ is maximized over the total variance $`\sigma^2_a + \sigma^2_e`$:

```math
L_{\text{prof}}(h^2) = \max_{\sigma^2_a + \sigma^2_e}\; L\left(h^2\left(\sigma^2_a + \sigma^2_e\right),\; \left(1 - h^2\right)\left(\sigma^2_a + \sigma^2_e\right)\right)
```

Plotting $`L_{\text{prof}}`$ against $`h^2`$ shows the shape of the likelihood around $`\hat{h}^2`$, including when the estimate is at the boundary.

### 11.4 Checks

1. **The MME route gives the same $`L`$ and score.** At the start values, $`L`$ and the score are recomputed with the formulas of [10.4](#104-computing-l-and-its-derivatives-from-the-mme). This includes $`\log\lvert\mathbf{H}\rvert`$, so both values of $`L`$ are equal, not only up to a constant. The code builds $`\mathbf{W}`$, $`\mathbf{T}`$, $`\mathbf{C}`$, $`\mathbf{C}^{-1}`$, $`\hat{\mathbf{a}}`$, $`\hat{\mathbf{e}}`$, $`\text{tr}(\mathbf{H}^{-1}\mathbf{C}^{22})`$ and $`\hat{\mathbf{a}}^\top\mathbf{H}^{-1}\hat{\mathbf{a}}`$, and checks agreement within $`10^{-6}`$.
2. **A general-purpose optimizer finds the same maximum.** `optim()` (L-BFGS-B, without derivatives, with $`\sigma^2_a, \sigma^2_e \ge 0`$) maximizes $`L`$ directly. AI-REML must reach at least the same $`L`$.
3. **The code recovers known variances.** 100 data sets are simulated with the same $`\mathbf{X}`$ and $`\mathbf{H}`$, with $`\sigma^2_a = 6`$ and $`\sigma^2_e = 9`$ ($`h^2 = 0.4`$). If $`\mathbf{z} \sim N(\mathbf{0}, \mathbf{I})`$ and $`\mathbf{U}_H`$ is the Cholesky factor of $`\mathbf{W}\mathbf{H}\mathbf{W}^\top`$ ($`\mathbf{W}\mathbf{H}\mathbf{W}^\top = \mathbf{U}_H^\top\mathbf{U}_H`$), then $`\mathbf{U}_H^\top\mathbf{z} \sim N(\mathbf{0}, \mathbf{W}\mathbf{H}\mathbf{W}^\top)`$. So

```math
\mathbf{y}_{\text{sim}} = \mathbf{X}_{\text{ord}}\mathbf{1}_p + \sigma_a\,\mathbf{U}_H^\top\mathbf{z}_1 + \sigma_e\,\mathbf{z}_2,
\qquad \mathbf{z}_1, \mathbf{z}_2 \sim N(\mathbf{0}, \mathbf{I}_N)
```

The fixed effects do not change the REML estimates, so all are set to 1. The mean of the 100 estimates should be close to the true values (REML is not exactly unbiased in small samples).

---

## 12. Conclusion

There are few phenotypes. With the current data, the REML estimate of $`\sigma^2_a`$ for `delta` goes to the boundary, so $`\hat{h}^2 \simeq 0`$. At that point $`\mathbf{V} = \mathbf{I}\sigma^2_e`$, and $`\hat{\sigma}^2_e`$ is the residual variance of the model with only the fixed effects. The profile log-likelihood decreases as $`h^2`$ increases.

With this amount of data, $`\sigma^2_a`$ (and therefore $`h^2`$) cannot be estimated reliably. Subsequent predictions of additive genetic values, and their propagation to SNP effects, would be unreliable and misleading.

At the same time, the checks in [11.4](#114-checks) show that the pipeline itself is correct. The two routes agree, a general-purpose optimizer reaches the same maximum, and the method recovers known variance components in simulated data with the same structure.

---

← [Sections 8–9](genomic-information-and-mme.md) · [References](../README.md#references) →
