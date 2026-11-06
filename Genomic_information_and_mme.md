# PAGssGBLUP: Genomic information and mixed model equations

Part of the methodology of [PAGssGBLUP](../README.md). Symbols are defined in [Notation](../README.md#notation), and references are listed in [References](../README.md#references).

← [Sections 1–7](../README.md) · [Sections 10–12](reml.md) →

---

## 8. Genomic information

### 8.1 The pedigree propagates genomic information

The pedigree relationships are the bridge from genotyped to non-genotyped animals. With $`\mathbb{E}[\mathbf{a}_1] = \mathbf{0}`$ and $`\mathbb{E}[\mathbf{a}_2] = \mathbf{0}`$:

```math
\mathbb{E}[\mathbf{a}_1 \mid \mathbf{a}_2] =
\mathbb{E}[\mathbf{a}_1] +
\Sigma_{\mathbf{a}_1\mathbf{a}_2}\Sigma_{\mathbf{a}_2\mathbf{a}_2}^{-1}\left(\mathbf{a}_2 - \mathbb{E}[\mathbf{a}_2]\right)
=
\Sigma_{\mathbf{a}_1\mathbf{a}_2}\Sigma_{\mathbf{a}_2\mathbf{a}_2}^{-1}\mathbf{a}_2
```

```math
\Sigma_{\mathbf{a}_1\mathbf{a}_2} = \text{Cov}\left(\mathbf{a}_1, \mathbf{a}_2^\top\right) = \mathbf{A}_{12}\sigma^2_a,
\qquad
\Sigma_{\mathbf{a}_2\mathbf{a}_2} = \text{Cov}\left(\mathbf{a}_2, \mathbf{a}_2^\top\right) = \mathbf{A}_{22}\sigma^2_a
```

These are the pedigree (co)variances: the link between $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$ comes from the pedigree. The genomic information enters through $`\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a`$ ([8.5](#85-variance-of-a)).

The additive effects of the genotyped animals are expressed through the markers:

```math
\mathbf{a}_2 = \mathbf{Zg}, \qquad \mathbf{Z} \in \mathbb{R}^{K_g \times J}, \qquad \mathbf{g} \in \mathbb{R}^{J \times 1}
```

The rows of $`\mathbf{Z}`$ are the animals of $`\mathbf{a}_2`$ ($`K_g = n_{a_2}`$). Genotyped animals without phenotypes (e.g., genotyped sires) could also be part of $`\mathbf{a}_2`$, although there are none in this data set.

### 8.2 Marker effects

The marker effects are independent with a common variance $`\sigma^2_g`$:

```math
g_j \sim N(0, \sigma^2_g) \;\; \forall j \in \lbrace 1, \dots, J \rbrace,
\qquad
\text{Cov}\left[g_j, g_{j'}\right] = 0 \;\; \forall j \neq j'
\qquad \Longleftrightarrow \qquad
\text{Var}[\mathbf{g}] = \mathbf{I}\sigma^2_g
```

Then

```math
\begin{aligned}
\text{Var}[\mathbf{Zg}]
&= \text{Cov}\left[\mathbf{Zg}, \mathbf{g}^\top\mathbf{Z}^\top\right]
= \mathbf{Z}\,\text{Cov}\left[\mathbf{g}, \mathbf{g}^\top\right]\mathbf{Z}^\top
= \mathbf{Z}\,\text{Var}[\mathbf{g}]\,\mathbf{Z}^\top \\
&= \mathbf{Z}\mathbf{I}\sigma^2_g\mathbf{Z}^\top
= \mathbf{Z}\mathbf{Z}^\top\sigma^2_g
\end{aligned}
```

### 8.3 Genomic relationship matrix G

To link $`\sigma^2_g`$ with $`\sigma^2_a`$, the dosage $`m_{ij}`$ is treated as a random variable: the genotype of a random animal from a population in HWE, with allele frequency $`p_j`$ and $`q_j = 1 - p_j`$.

```math
\begin{aligned}
\mathbb{E}\left[m_{ij}\right] &= 2p_j
\quad \Longrightarrow \quad \left(\mathbb{E}\left[m_{ij}\right]\right)^2 = 4p_j^2 \\[4pt]
\mathbb{E}\left[m_{ij}^2\right]
&= 0^2 \times (1-p_j)^2 + 1^2 \times 2(1-p_j)p_j + 2^2 \times p_j^2
= 2p_j - 2p_j^2 + 4p_j^2
= 2p_j + 2p_j^2 \\[4pt]
\text{Var}\left[m_{ij}\right]
&= \mathbb{E}\left[m_{ij}^2\right] - \left(\mathbb{E}\left[m_{ij}\right]\right)^2
= 2p_j + 2p_j^2 - 4p_j^2
= 2p_j(1 - p_j)
= 2p_jq_j
\end{aligned}
```

The centered dosage $`z_{ij} = m_{ij} - 2p_j`$ only shifts $`m_{ij}`$ by a constant, so $`\mathbb{E}[z_{ij}] = 0`$ and $`\text{Var}[z_{ij}] = \mathbb{E}[z_{ij}^2] = 2p_jq_j`$.

The additive genetic value of animal $`i`$ (element $`i`$ of $`\mathbf{Zg}`$) is

```math
a_i = \sum_{j=1}^{J} z_{ij}g_j
```

The marker effects are assumed independent of the genotypes. Because the $`g_j`$ have mean zero and are uncorrelated, terms of different markers are uncorrelated, even if their dosages are correlated (linkage disequilibrium):

```math
\text{Cov}\left[z_{ij}g_j, z_{ij'}g_{j'}\right] = \mathbb{E}\left[z_{ij}z_{ij'}\right]\mathbb{E}\left[g_jg_{j'}\right] = 0 \;\; \forall j \neq j'
```

```math
\text{Var}[a_i]
= \sum_{j=1}^{J}\text{Var}\left[z_{ij}g_j\right]
= \sum_{j=1}^{J}\mathbb{E}\left[z_{ij}^2\right]\mathbb{E}\left[g_j^2\right]
= \sum_{j=1}^{J} 2p_jq_j\sigma^2_g
= 2\sum_{j=1}^{J}p_jq_j\;\sigma^2_g
```

With $`\sigma^2_a = \text{Var}[a_i]`$:

```math
\sigma^2_a = 2\sum_{j=1}^{J}p_jq_j\;\sigma^2_g
\quad \Longrightarrow \quad
\sigma^2_g = \frac{\sigma^2_a}{2\sum_{j=1}^{J}p_jq_j}
```

With $`\mathbf{Z}`$ fixed (the observed genotypes):

```math
\text{Var}[\mathbf{a}_2] = \text{Var}[\mathbf{Zg}]
= \mathbf{Z}\mathbf{Z}^\top\sigma^2_g
= \frac{\mathbf{Z}\mathbf{Z}^\top}{2\sum_{j=1}^{J}p_jq_j}\,\sigma^2_a
= \mathbf{G}\sigma^2_a,
\qquad
\mathbf{G} = \frac{\mathbf{Z}\mathbf{Z}^\top}{2\sum_{j=1}^{J}p_jq_j} \in \mathbb{R}^{K_g \times K_g}
```

This is the genomic relationship matrix of VanRaden (2008). It connects the marker model ($`\mathbf{g}`$, $`\sigma^2_g`$) with the animal model ($`\mathbf{a}`$, $`\sigma^2_a`$). The $`p_j`$ should be base-population frequencies. Because those are unknown, the frequencies of the genotyped animals are used (the same $`p_j`$ used to center $`\mathbf{Z}`$). For genotyped animals, $`\mathbf{G}`$ replaces $`\mathbf{A}_{22}`$.

**Computation.** $`\mathbf{Z}\mathbf{Z}^\top`$ is computed with `tcrossprod(Z)` (one BLAS call, without forming $`\mathbf{Z}^\top`$), and divided by $`2\sum_j p_jq_j`$.

**Checks.**

1. $`\mathbf{G}`$ is symmetric.
2. Its rows and columns follow $`\mathbf{a}_2`$, the same order as $`\mathbf{A}_{22}`$.
3. Every row adds up to zero: each column of $`\mathbf{Z}`$ has mean zero, so $`\mathbf{Z}^\top\mathbf{1} = \mathbf{0}`$ and $`\mathbf{G}\mathbf{1} = \mathbf{Z}(\mathbf{Z}^\top\mathbf{1})/(2\sum_j p_jq_j) = \mathbf{0}`$.

**$`\mathbf{G}`$ is singular.** Since $`\mathbf{G}\mathbf{1} = \mathbf{0}`$, $`\text{rank}(\mathbf{G}) \le K_g - 1`$ and its smallest eigenvalue is zero (up to rounding). $`\mathbf{G}^{-1}`$ does not exist, but $`\mathbf{H}^{-1}`$ needs it.

### 8.4 Blending G

```math
\mathbf{G}_w = (1 - w)\,\mathbf{G} + w\,\mathbf{A}_{22}, \qquad w = 0.05
```

($`w = 0.05`$ is a common default, e.g., in BLUPF90.)

**Meaning of $`w`$.** The markers may not capture all the additive genetic variance. The additive effect of a genotyped animal is split into a part explained by the markers and a residual polygenic part $`\mathbf{u}`$, which follows the pedigree and is independent of $`\mathbf{g}`$:

```math
\mathbf{a}_2 = \mathbf{Zg} + \mathbf{u},
\qquad
\text{Var}[\mathbf{Zg}] = (1-w)\,\mathbf{G}\sigma^2_a,
\qquad
\text{Var}[\mathbf{u}] = w\,\mathbf{A}_{22}\sigma^2_a
```

```math
\Longrightarrow \quad
\text{Var}[\mathbf{a}_2] = (1-w)\,\mathbf{G}\sigma^2_a + w\,\mathbf{A}_{22}\sigma^2_a = \mathbf{G}_w\sigma^2_a
```

The markers explain a proportion $`1 - w`$ of $`\sigma^2_a`$ and the pedigree the remaining $`w`$. With $`w = 0`$, $`\mathbf{a}_2 = \mathbf{Zg}`$ as in [8.3](#83-genomic-relationship-matrix-g).

$`\mathbf{G}_w`$ is positive definite. Since $`\mathbf{G}`$ is positive semi-definite, the smallest eigenvalue of $`\mathbf{G}_w`$ is at least $`w`$ times the smallest eigenvalue of $`\mathbf{A}_{22}`$ (checked in the code). Its inverse is computed from the Cholesky factor (`chol2inv(chol(Gw))`). This is faster and more stable than `solve()` for a symmetric positive definite matrix, and `chol()` stops if $`\mathbf{G}_w`$ is not positive definite. The code checks that $`\mathbf{G}_w\mathbf{G}_w^{-1} = \mathbf{I}`$.

**Effect of blending on what follows.** The derivation of $`\text{Var}[\mathbf{a}]`$ below uses only that $`\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a`$ with $`\mathbf{G}`$ symmetric, not how $`\mathbf{G}`$ is built. Every step therefore holds with $`\mathbf{G}_w`$ in place of $`\mathbf{G}`$. **From here on, $`\mathbf{G}`$ in the formulas stands for the matrix actually used, $`\mathbf{G}_w`$.** Since

```math
\mathbf{G}_w - \mathbf{A}_{22} = (1-w)\,\mathbf{G} + w\,\mathbf{A}_{22} - \mathbf{A}_{22} = (1-w)\left(\mathbf{G} - \mathbf{A}_{22}\right)
```

blending shrinks the genomic correction of the non-genotyped animals by $`(1-w)`$:

```math
\text{Var}[\mathbf{a}_1] =
\left[\mathbf{A}_{11} + (1-w)\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\left(\mathbf{G} - \mathbf{A}_{22}\right)\mathbf{A}_{22}^{-1}\mathbf{A}_{21}\right]\sigma^2_a
```

In the same way, the whole matrix is $`\mathbf{H}_w = (1-w)\,\mathbf{H} + w\,\mathbf{A}`$: a weighted average of $`\mathbf{H}`$ built with the unblended $`\mathbf{G}`$ and the pedigree matrix $`\mathbf{A}`$. In $`\mathbf{H}^{-1}`$, $`\mathbf{G}^{-1}`$ is $`\mathbf{G}_w^{-1}`$.

Blending fixes the singularity of $`\mathbf{G}`$, but not the difference in scale between $`\mathbf{G}`$ and $`\mathbf{A}_{22}`$. Tuning $`\mathbf{G}`$ to $`\mathbf{A}_{22}`$ is a separate step and is not done here.

### 8.5 Variance of a

**Naive approach (incorrect).** Replace $`\mathbf{A}_{22}`$ by $`\mathbf{G}`$ and keep the other blocks of $`\mathbf{A}`$:

```math
\text{Var}\left[\begin{pmatrix}\mathbf{a}_1 \\ \mathbf{a}_2\end{pmatrix}\right]
\overset{?}{=}
\begin{pmatrix}
\mathbf{A}_{11} & \mathbf{A}_{12} \\
\mathbf{A}_{21} & \mathbf{G}
\end{pmatrix}\sigma^2_a
=
\left[\mathbf{A} +
\begin{pmatrix}
\mathbf{0}_{11} & \mathbf{0}_{12} \\
\mathbf{0}_{21} & \mathbf{G} - \mathbf{A}_{22}
\end{pmatrix}\right]\sigma^2_a
```

This is wrong. The non-genotyped animals are related to the genotyped ones through the pedigree, so when $`\text{Var}[\mathbf{a}_2]`$ changes from $`\mathbf{A}_{22}\sigma^2_a`$ to $`\mathbf{G}\sigma^2_a`$, then $`\text{Var}[\mathbf{a}_1]`$ and $`\text{Cov}[\mathbf{a}_1, \mathbf{a}_2^\top]`$ must change too. The correct blocks are derived below.

**Decomposition of $`\mathbf{a}_1`$.** From [8.1](#81-the-pedigree-propagates-genomic-information):

```math
\mathbb{E}[\mathbf{a}_1 \mid \mathbf{a}_2]
= \mathbf{A}_{12}\sigma^2_a\left(\mathbf{A}_{22}\sigma^2_a\right)^{-1}\mathbf{a}_2
= \mathbf{A}_{12}\sigma^2_a\mathbf{A}_{22}^{-1}\frac{1}{\sigma^2_a}\mathbf{a}_2
= \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{a}_2
= \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{Zg}
```

Then $`\mathbf{a}_1`$ is its conditional expectation plus a deviation $`\boldsymbol{\varepsilon}`$:

```math
\mathbf{a}_1 = \mathbb{E}[\mathbf{a}_1 \mid \mathbf{a}_2] + \boldsymbol{\varepsilon}
= \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{a}_2 + \boldsymbol{\varepsilon}
= \mathbf{D}\mathbf{a}_2 + \boldsymbol{\varepsilon}
\quad \Longrightarrow \quad
\boldsymbol{\varepsilon} = \mathbf{a}_1 - \mathbf{D}\mathbf{a}_2,
\qquad
\mathbf{D} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1} \in \mathbb{R}^{n_{a_1} \times n_{a_2}}
```

$`\mathbf{D}\mathbf{a}_2`$ is the part of $`\mathbf{a}_1`$ predicted from the genotyped animals through the pedigree. $`\boldsymbol{\varepsilon}`$ is the part that $`\mathbf{a}_2`$ cannot explain (e.g., the Mendelian sampling of the non-genotyped animals). Note that $`\mathbf{D}^\top = (\mathbf{A}_{22}^{-1})^\top\mathbf{A}_{12}^\top = \mathbf{A}_{22}^{-1}\mathbf{A}_{21}`$, because $`\mathbf{A}_{22}^{-1}`$ is symmetric and $`\mathbf{A}_{12}^\top = \mathbf{A}_{21}`$.

**Key assumption (Legarra et al., 2009).** The genomic information changes only the distribution of $`\mathbf{a}_2`$: $`\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a`$ instead of $`\mathbf{A}_{22}\sigma^2_a`$. The link between $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$ ($`\mathbf{D}`$ and $`\boldsymbol{\varepsilon}`$) is still the one given by the pedigree. $`\boldsymbol{\varepsilon}`$ and its variance are therefore computed with the pedigree (co)variances, written $`\text{Var}_{\mathbf{A}}`$ and $`\text{Cov}_{\mathbf{A}}`$:

```math
\text{Var}_{\mathbf{A}}[\mathbf{a}_1] = \mathbf{A}_{11}\sigma^2_a, \quad
\text{Var}_{\mathbf{A}}[\mathbf{a}_2] = \mathbf{A}_{22}\sigma^2_a, \quad
\text{Cov}_{\mathbf{A}}\left[\mathbf{a}_1, \mathbf{a}_2^\top\right] = \mathbf{A}_{12}\sigma^2_a, \quad
\text{Cov}_{\mathbf{A}}\left[\mathbf{a}_2, \mathbf{a}_1^\top\right] = \mathbf{A}_{21}\sigma^2_a
```

#### Covariance between $`\mathbf{D}\mathbf{a}_2`$ and $`\boldsymbol{\varepsilon}`$

```math
\begin{aligned}
\text{Cov}_{\mathbf{A}}\left[\mathbf{a}_2, \boldsymbol{\varepsilon}^\top\right]
&= \text{Cov}_{\mathbf{A}}\left[\mathbf{a}_2, \left(\mathbf{a}_1 - \mathbf{D}\mathbf{a}_2\right)^\top\right]
= \text{Cov}_{\mathbf{A}}\left[\mathbf{a}_2, \mathbf{a}_1^\top\right] - \text{Cov}_{\mathbf{A}}\left[\mathbf{a}_2, \mathbf{a}_2^\top\right]\mathbf{D}^\top \\
&= \mathbf{A}_{21}\sigma^2_a - \mathbf{A}_{22}\sigma^2_a\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
= \mathbf{A}_{21}\sigma^2_a - \mathbf{A}_{21}\sigma^2_a
= \mathbf{0} \in \mathbb{R}^{n_{a_2} \times n_{a_1}}
\end{aligned}
```

By the key assumption, $`\boldsymbol{\varepsilon}`$ stays independent of $`\mathbf{a}_2`$ when $`\text{Var}[\mathbf{a}_2]`$ changes to $`\mathbf{G}\sigma^2_a`$. Then

```math
\text{Cov}\left[\mathbf{D}\mathbf{a}_2, \boldsymbol{\varepsilon}^\top\right]
= \text{Cov}\left[\boldsymbol{\varepsilon}, \left(\mathbf{D}\mathbf{a}_2\right)^\top\right]^\top
= \mathbf{0} \in \mathbb{R}^{n_{a_1} \times n_{a_1}}
\qquad \Longrightarrow \qquad
\text{Var}[\mathbf{a}_1] = \text{Var}[\mathbf{D}\mathbf{a}_2 + \boldsymbol{\varepsilon}] = \text{Var}[\mathbf{D}\mathbf{a}_2] + \text{Var}[\boldsymbol{\varepsilon}]
```

#### Variance of $`\mathbf{D}\mathbf{a}_2`$

Here $`\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a`$ (genomic information):

```math
\begin{aligned}
\text{Var}[\mathbf{D}\mathbf{a}_2]
&= \text{Cov}\left[\mathbf{D}\mathbf{a}_2, \left(\mathbf{D}\mathbf{a}_2\right)^\top\right]
= \mathbf{D}\,\text{Cov}\left[\mathbf{a}_2, \mathbf{a}_2^\top\right]\mathbf{D}^\top
= \mathbf{D}\,\text{Var}[\mathbf{a}_2]\,\mathbf{D}^\top
= \mathbf{D}\mathbf{G}\sigma^2_a\mathbf{D}^\top \\
&= \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\left(\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\right)^\top
= \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
\end{aligned}
```

#### Variance of $`\boldsymbol{\varepsilon}`$

Here all (co)variances are the pedigree ones (key assumption):

```math
\begin{aligned}
\text{Var}[\boldsymbol{\varepsilon}]
&= \text{Var}_{\mathbf{A}}\left[\mathbf{a}_1 - \mathbf{D}\mathbf{a}_2\right] \\
&= \text{Var}_{\mathbf{A}}[\mathbf{a}_1]
+ \text{Var}_{\mathbf{A}}[-\mathbf{D}\mathbf{a}_2]
+ \text{Cov}_{\mathbf{A}}\left[\mathbf{a}_1, \left(-\mathbf{D}\mathbf{a}_2\right)^\top\right]
+ \text{Cov}_{\mathbf{A}}\left[-\mathbf{D}\mathbf{a}_2, \mathbf{a}_1^\top\right]
\end{aligned}
```

The four terms:

```math
\begin{aligned}
\text{Var}_{\mathbf{A}}[\mathbf{a}_1] &= \mathbf{A}_{11}\sigma^2_a \\[4pt]
\text{Var}_{\mathbf{A}}[-\mathbf{D}\mathbf{a}_2]
&= \text{Cov}_{\mathbf{A}}\left[-\mathbf{D}\mathbf{a}_2, \mathbf{a}_2^\top\mathbf{D}^\top(-1)\right]
= (-1)\,\mathbf{D}\,\text{Var}_{\mathbf{A}}[\mathbf{a}_2]\,\mathbf{D}^\top(-1)
= \mathbf{D}\mathbf{A}_{22}\sigma^2_a\mathbf{D}^\top \\
&= \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{22}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
= \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} \\[4pt]
\text{Cov}_{\mathbf{A}}\left[\mathbf{a}_1, \mathbf{a}_2^\top\mathbf{D}^\top(-1)\right]
&= \text{Cov}_{\mathbf{A}}\left[\mathbf{a}_1, \mathbf{a}_2^\top\right]\mathbf{D}^\top(-1)
= -\sigma^2_a\,\mathbf{A}_{12}\mathbf{D}^\top
= -\sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} \\[4pt]
\text{Cov}_{\mathbf{A}}\left[-\mathbf{D}\mathbf{a}_2, \mathbf{a}_1^\top\right]
&= (-1)\,\mathbf{D}\,\text{Cov}_{\mathbf{A}}\left[\mathbf{a}_2, \mathbf{a}_1^\top\right]
= -\sigma^2_a\,\mathbf{D}\mathbf{A}_{21}
= -\sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
\end{aligned}
```

Adding them:

```math
\text{Var}[\boldsymbol{\varepsilon}]
= \mathbf{A}_{11}\sigma^2_a
+ \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
- \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
- \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
= \left(\mathbf{A}_{11} - \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}\right)\sigma^2_a
```

#### Variance of $`\mathbf{a}_1`$

```math
\begin{aligned}
\text{Var}[\mathbf{a}_1]
&= \text{Var}[\mathbf{D}\mathbf{a}_2] + \text{Var}[\boldsymbol{\varepsilon}]
= \sigma^2_a\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
+ \left(\mathbf{A}_{11} - \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}\right)\sigma^2_a \\
&= \left(\mathbf{A}_{11}
+ \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}
- \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{22}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}\right)\sigma^2_a \\
&= \left[\mathbf{A}_{11} + \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\left(\mathbf{G} - \mathbf{A}_{22}\right)\mathbf{A}_{22}^{-1}\mathbf{A}_{21}\right]\sigma^2_a
= \mathbf{H}_{11}\sigma^2_a
\end{aligned}
```

The second line uses $`\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{22}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}`$ to factor the last two terms. Interpretation:

- $`\mathbf{A}_{11}`$: pedigree relationships among non-genotyped animals.
- $`\mathbf{G} - \mathbf{A}_{22}`$: what the markers add to (or remove from) the pedigree relationships among genotyped animals.
- $`\mathbf{A}_{12}\mathbf{A}_{22}^{-1}`$: propagates this difference to the non-genotyped animals through their pedigree relationships with the genotyped animals.
- If $`\mathbf{G} = \mathbf{A}_{22}`$, then $`\text{Var}[\mathbf{a}_1] = \mathbf{A}_{11}\sigma^2_a`$ (pedigree only).

#### Variance of $`\mathbf{a}_2`$

```math
\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a
```

#### Covariance between $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$

```math
\text{Cov}\left[\mathbf{a}_1, \mathbf{a}_2^\top\right]
= \text{Cov}\left[\mathbf{D}\mathbf{a}_2 + \boldsymbol{\varepsilon}, \mathbf{a}_2^\top\right]
= \text{Cov}\left[\mathbf{D}\mathbf{a}_2, \mathbf{a}_2^\top\right] + \text{Cov}\left[\boldsymbol{\varepsilon}, \mathbf{a}_2^\top\right]
```

```math
\text{Cov}\left[\mathbf{D}\mathbf{a}_2, \mathbf{a}_2^\top\right]
= \mathbf{D}\,\text{Cov}\left[\mathbf{a}_2, \mathbf{a}_2^\top\right]
= \mathbf{D}\,\text{Var}[\mathbf{a}_2]
= \mathbf{D}\mathbf{G}\sigma^2_a
= \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\sigma^2_a
```

```math
\begin{aligned}
\text{Cov}_{\mathbf{A}}\left[\boldsymbol{\varepsilon}, \mathbf{a}_2^\top\right]
&= \text{Cov}_{\mathbf{A}}\left[\mathbf{a}_1 - \mathbf{D}\mathbf{a}_2, \mathbf{a}_2^\top\right]
= \text{Cov}_{\mathbf{A}}\left[\mathbf{a}_1, \mathbf{a}_2^\top\right] + \text{Cov}_{\mathbf{A}}\left[-\mathbf{D}\mathbf{a}_2, \mathbf{a}_2^\top\right] \\
&= \sigma^2_a\mathbf{A}_{12} - \sigma^2_a\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{22}
= \sigma^2_a\mathbf{A}_{12} - \sigma^2_a\mathbf{A}_{12}
= \mathbf{0} \in \mathbb{R}^{n_{a_1} \times n_{a_2}}
\end{aligned}
```

This is the transpose of $`\text{Cov}_{\mathbf{A}}[\mathbf{a}_2, \boldsymbol{\varepsilon}^\top] = \mathbf{0}`$, and by the key assumption it also holds with genomic information. Then

```math
\text{Cov}\left[\mathbf{a}_1, \mathbf{a}_2^\top\right] = \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\sigma^2_a + \mathbf{0} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\sigma^2_a
```

```math
\text{Cov}\left[\mathbf{a}_2, \mathbf{a}_1^\top\right]
= \left(\text{Cov}\left[\mathbf{a}_1, \mathbf{a}_2^\top\right]\right)^\top
= \left(\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}\sigma^2_a\right)^\top
= \mathbf{G}^\top\left(\mathbf{A}_{22}^{-1}\right)^\top\mathbf{A}_{12}^\top\sigma^2_a
= \mathbf{G}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}\sigma^2_a
```

because $`\mathbf{G}`$ and $`\mathbf{A}_{22}^{-1}`$ are symmetric and $`\mathbf{A}_{12}^\top = \mathbf{A}_{21}`$.

#### Full expression: the H matrix

```math
\text{Var}[\mathbf{a}] =
\text{Var}\left[\begin{pmatrix}\mathbf{a}_1 \\ \mathbf{a}_2\end{pmatrix}\right]
=
\begin{pmatrix}
\mathbf{A}_{11} + \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\left(\mathbf{G} - \mathbf{A}_{22}\right)\mathbf{A}_{22}^{-1}\mathbf{A}_{21} &
\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G} \\
\mathbf{G}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} &
\mathbf{G}
\end{pmatrix}\sigma^2_a
=
\begin{pmatrix}
\mathbf{H}_{11} & \mathbf{H}_{12} \\
\mathbf{H}_{21} & \mathbf{H}_{22}
\end{pmatrix}\sigma^2_a
= \mathbf{H}\sigma^2_a
```

```math
\mathbf{H}_{11} = \mathbf{A}_{11} + \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\left(\mathbf{G} - \mathbf{A}_{22}\right)\mathbf{A}_{22}^{-1}\mathbf{A}_{21},
\qquad
\mathbf{H}_{12} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G},
\qquad
\mathbf{H}_{21} = \mathbf{G}\mathbf{A}_{22}^{-1}\mathbf{A}_{21},
\qquad
\mathbf{H}_{22} = \mathbf{G}
```

$`\mathbf{H}_{21} = \mathbf{H}_{12}^\top`$, so $`\mathbf{H}`$ is symmetric. $`\mathbf{H}`$ is the relationship matrix that combines pedigree and genomic relationships (Legarra et al., 2009), with the animals in the order of $`\mathbf{A}_{\text{ord}}`$. In compact form, with $`\mathbf{D} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1}`$:

```math
\mathbf{H} = \mathbf{A} + \begin{pmatrix}\mathbf{D} \\ \mathbf{I}\end{pmatrix}\left(\mathbf{G} - \mathbf{A}_{22}\right)\begin{pmatrix}\mathbf{D} \\ \mathbf{I}\end{pmatrix}^\top
```

Multiplying out gives the four blocks above.

### 8.6 Inverse of H

The MME need $`\mathbf{H}^{-1}`$, not $`\mathbf{H}`$. Aguilar et al. (2010) showed that

```math
\mathbf{H}^{-1} = \mathbf{A}^{-1} +
\begin{pmatrix}
\mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{G}^{-1} - \mathbf{A}_{22}^{-1}
\end{pmatrix}
```

Only the block of the genotyped animals changes. Three inverses are needed:

1. $`\mathbf{A}^{-1}`$: all animals, in the order of $`\mathbf{A}_{\text{ord}}`$.
2. $`\mathbf{A}_{22}^{-1}`$: the inverse of the block $`\mathbf{A}_{22}`$. This is **not** the 22 block of $`\mathbf{A}^{-1}`$. It is the plain pedigree block and is **not** blended.
3. $`\mathbf{G}^{-1} = \mathbf{G}_w^{-1}`$: the inverse of the blended matrix ([8.4](#84-blending-g)).

$`\mathbf{A}`$ and $`\mathbf{A}_{22}`$ are symmetric positive definite and are inverted from their Cholesky factors. Here $`\mathbf{A}^{-1}`$ is obtained by inverting the dense $`\mathbf{A}`$, which is fine for this pedigree. For large pedigrees, it is built directly from the pedigree (Henderson's rules), because it is sparse.

**Checks.**

1. $`\mathbf{H}^{-1}`$ is symmetric, with the animals in the order of $`\mathbf{A}_{\text{ord}}`$.
2. The rows of the non-genotyped animals are identical to those of $`\mathbf{A}^{-1}`$.
3. $`\mathbf{H}`$ is built from the compact form in [8.5](#85-variance-of-a) with $`\mathbf{G} = \mathbf{G}_w`$. The code checks that $`\mathbf{H}_{22} = \mathbf{G}_w`$ and that $`\max\lvert\mathbf{H}^{-1}\mathbf{H} - \mathbf{I}\rvert < 10^{-8}`$. $`\mathbf{H}`$ is kept, because REML uses it.

---

## 9. Mixed model equations

The model is $`\mathbf{y} = \mathbf{Xb} + \mathbf{Wa} + \mathbf{e}`$ with

```math
\mathbf{a} \sim N\left(\mathbf{0}, \mathbf{H}\sigma^2_a\right), \quad \mathbf{H}\sigma^2_a = \mathbf{Q} \in \mathbb{R}^{K \times K},
\qquad
\mathbf{e} \sim N\left(\mathbf{0}, \mathbf{I}\sigma^2_e\right), \quad \mathbf{I}\sigma^2_e = \mathbf{R} \in \mathbb{R}^{N \times N},
\qquad
\text{Cov}\left[\mathbf{a}, \mathbf{e}^\top\right] = \mathbf{0}
```

where $`K`$ is the number of animals in $`\mathbf{a}`$ (all animals in the pedigree) and $`N`$ the number of phenotypes. The MME (Henderson, 1984) are

```math
\begin{pmatrix}
\mathbf{X}^\top\mathbf{R}^{-1}\mathbf{X} & \mathbf{X}^\top\mathbf{R}^{-1}\mathbf{W} \\
\mathbf{W}^\top\mathbf{R}^{-1}\mathbf{X} & \mathbf{W}^\top\mathbf{R}^{-1}\mathbf{W} + \mathbf{H}^{-1}\frac{1}{\sigma^2_a}
\end{pmatrix}
\begin{pmatrix} \widehat{\mathbf{b}} \\ \widehat{\mathbf{a}} \end{pmatrix}
=
\begin{pmatrix} \mathbf{X}^\top\mathbf{R}^{-1}\mathbf{y} \\ \mathbf{W}^\top\mathbf{R}^{-1}\mathbf{y} \end{pmatrix}
```

with $`\mathbf{H}^{-1}\frac{1}{\sigma^2_a} = \left(\mathbf{H}\sigma^2_a\right)^{-1} = \mathbf{Q}^{-1}`$ and $`\mathbf{R}^{-1} = \left(\mathbf{I}\sigma^2_e\right)^{-1} = \mathbf{I}\frac{1}{\sigma^2_e}`$. Then

```math
\begin{pmatrix}
\mathbf{X}^\top\frac{1}{\sigma^2_e}\mathbf{X} & \mathbf{X}^\top\frac{1}{\sigma^2_e}\mathbf{W} \\
\mathbf{W}^\top\frac{1}{\sigma^2_e}\mathbf{X} & \mathbf{W}^\top\frac{1}{\sigma^2_e}\mathbf{W} + \mathbf{H}^{-1}\frac{1}{\sigma^2_a}
\end{pmatrix}
\begin{pmatrix} \widehat{\mathbf{b}} \\ \widehat{\mathbf{a}} \end{pmatrix}
=
\begin{pmatrix} \mathbf{X}^\top\frac{1}{\sigma^2_e}\mathbf{y} \\ \mathbf{W}^\top\frac{1}{\sigma^2_e}\mathbf{y} \end{pmatrix}
```

Multiplying both sides by $`\sigma^2_e`$:

```math
\begin{pmatrix}
\mathbf{X}^\top\mathbf{X} & \mathbf{X}^\top\mathbf{W} \\
\mathbf{W}^\top\mathbf{X} & \mathbf{W}^\top\mathbf{W} + \mathbf{H}^{-1}\alpha
\end{pmatrix}
\begin{pmatrix} \widehat{\mathbf{b}} \\ \widehat{\mathbf{a}} \end{pmatrix}
=
\begin{pmatrix} \mathbf{X}^\top\mathbf{y} \\ \mathbf{W}^\top\mathbf{y} \end{pmatrix},
\qquad
\alpha = \frac{\sigma^2_e}{\sigma^2_a} = \frac{1 - h^2}{h^2},
\qquad
h^2 = \frac{\sigma^2_a}{\sigma^2_a + \sigma^2_e}
```

$`\mathbf{W} \in \mathbb{R}^{N \times K}`$ has one 1 per row, in the column of the animal with the record. The columns of animals without a record ($`A_{100}`$, $`A_{101}`$) are zero.

---

← [Sections 1–7](../README.md) · [Sections 10–12](reml.md) →
