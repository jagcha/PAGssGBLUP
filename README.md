# PAGssGBLUP

Single-step genomic BLUP (ssGBLUP) analysis of the change in pregnancy-associated glycoproteins (PAGs) between day 16 and day 26 in dairy cows, with variance components estimated by average-information REML (AI-REML).

The full analysis is in [`code/PAGssGBLUP.Rmd`](code/PAGssGBLUP.Rmd). This README and the two pages in [`docs/`](docs/) describe the theory and the methodology applied there, step by step and in the same notation. The text is split across three pages because GitHub renders only a limited number of math expressions per page. Section names in *italics* refer to sections of the Rmd.

> **Data are private and are not included in this repository.** The Rmd can be read as a worked example. To run it, the data files described in [How to run](#how-to-run) are needed.

---

## Contents

1. [Repository structure](#repository-structure)
2. [How to run](#how-to-run)
3. [Notation](#notation)
4. [Phenotype data](#1-phenotype-data)
5. [Working model](#2-working-model)
6. [Genotype data](#3-genotype-data)
7. [Pedigree](#4-pedigree)
8. [Numerator relationship matrix A](#5-numerator-relationship-matrix-a)
9. [Reorder and partition A by sets of animals](#6-reorder-and-partition-a-by-sets-of-animals)
10. [Align data with the ordered A](#7-align-data-with-the-ordered-a)
11. [Genomic information](docs/genomic-information-and-mme.md#8-genomic-information)
12. [Mixed model equations](docs/genomic-information-and-mme.md#9-mixed-model-equations)
13. [Restricted maximum likelihood](docs/reml.md#10-restricted-maximum-likelihood)
14. [AI-REML implementation](docs/reml.md#11-ai-reml-implementation)
15. [Conclusion](docs/reml.md#12-conclusion)
16. [References](#references)

---

## Repository structure

```
PAGssGBLUP/
├── README.md
├── .gitignore          # excludes data/, rendered HTML and R session files
├── code/
│   └── PAGssGBLUP.Rmd  # full analysis
└── docs/
    ├── genomic-information-and-mme.md  # sections 8–9
    └── reml.md                         # sections 10–12
```

## How to run

The Rmd reads its inputs with paths relative to `code/`, so the data must be in a `data/` folder next to `code/` (not included here):

```
PAGssGBLUP/
├── code/PAGssGBLUP.Rmd
└── data/
    ├── Main_Data.xlsx          # phenotypes: one row per cow
    ├── imputed_genotypes.txt   # genotypes: ID, two numeric fields, one string of digits (one per SNP)
    ├── Keys_to_ID.txt          # map from genotype ID (Key) to pedigree ID (ID); first line blank
    └── 5gen_ped.txt            # pedigree, CDCB Format 1 (fixed width, 99 characters per line)
```

Required R packages: `openxlsx`, `readr`, `data.table`, `HardyWeinberg`, `nadiv` (used only as an independent check of $`\mathbf{A}`$), and `rmarkdown`/`knitr` to knit. Knit the Rmd from RStudio, or run its chunks with `code/` as the working directory.

---

## Notation

Matrices are bold upper case, vectors bold lower case, scalars in italics. $`^\top`$ is the transpose.

| Symbol | Meaning |
|:---|:---|
| $`N = n_{y_1} + n_{y_2}`$ | Number of phenotypes (one per cow), split in non-genotyped ($`n_{y_1}`$) and genotyped ($`n_{y_2}`$) cows |
| $`K = n_{a_1} + n_{a_2}`$ | Number of animals in the pedigree (levels of $`\mathbf{a}`$), split in non-genotyped ($`n_{a_1}`$) and genotyped ($`n_{a_2}`$) |
| $`K_g`$ | Number of genotyped animals after quality control; $`K_g = n_{a_2}`$ |
| $`J`$ | Number of markers (SNPs) after quality control |
| $`p`$ | Number of fixed-effect parameters (columns of $`\mathbf{X}`$); `n_b` in the code |
| $`p_j`$, $`q_j = 1 - p_j`$ | Frequency of the alternative allele of marker $`j`$, and of the other allele; `p` in the code |
| $`\mathbf{y}`$, $`\mathbf{b}`$, $`\mathbf{a}`$, $`\mathbf{e}`$ | Phenotypes, fixed effects, additive genetic effects, residuals |
| $`\mathbf{X}`$, $`\mathbf{W}`$ | Design matrices of $`\mathbf{b}`$ and $`\mathbf{a}`$ |
| $`\mathbf{M}`$, $`\mathbf{Z}`$ | Raw SNP dosages, and dosages centered by $`2p_j`$ |
| $`\mathbf{g}`$ | Additive effects of the markers |
| $`\mathbf{A}`$, $`\mathbf{G}`$, $`\mathbf{H}`$ | Pedigree, genomic, and combined (pedigree + genomic) relationship matrices |
| $`\mathbf{G}_w`$ | Blended genomic matrix, $`(1-w)\mathbf{G} + w\mathbf{A}_{22}`$ |
| $`\sigma^2_a`$, $`\sigma^2_e`$, $`\sigma^2_g`$ | Additive genetic, residual, and per-marker variance |
| $`\alpha = \sigma^2_e/\sigma^2_a`$ | Variance ratio in the MME |
| $`h^2 = \sigma^2_a/(\sigma^2_a + \sigma^2_e)`$ | Heritability |
| $`\boldsymbol{\theta} = (\sigma^2_a, \sigma^2_e)^\top`$ | Variance components estimated by REML |
| $`\mathbf{V}`$, $`\mathbf{P}`$ | Variance of $`\mathbf{y}`$, and REML projection matrix |
| $`L`$ | REML log-likelihood |
| $`A_{100}, A_{110}, A_{101}, A_{111}`$ | Sets of animals (pedigree, phenotype, genotype; 1 = yes, 0 = no) |

Subscript 1 always refers to animals **without** genotype, and subscript 2 to animals **with** genotype.

Three symbols are used in two places, and the context separates them. $`p`$ (number of columns of $`\mathbf{X}`$) is not $`p_j`$ (allele frequency). $`F_i`$ (inbreeding coefficient) is not $`\mathbf{F}`$ (the $`N \times 2`$ matrix in the AI matrix). $`\mathbf{P}`$ (projection matrix) is not $`\mathbb{P}`$ (probability). The Rmd also uses $`\alpha`$ for the Hardy-Weinberg threshold and $`\mathbf{g}`$ for the gradient of $`h^2`$; here these are written $`\alpha_{\text{HWE}}`$ and $`\nabla h^2`$ to avoid confusion with the variance ratio and the marker effects.

---

## 1. Phenotype data

Each row of the phenotype file is one cow. The cow ID used everywhere is `'HO'` + `EID`, which is the ID of the cow in the pedigree; it links phenotypes, pedigree and genotypes. There are no repeated records (one record per cow; checked in the code).

**Trait.** The change in PAGs between day 16 and day 26:

```math
\text{delta} = \text{PAGs}_{d26} - \text{PAGs}_{d16}
```

**Pregnancy status.** From the pregnancy diagnoses at day 32 (`PIA32`) and day 60 (`PIA60`):

| Status | Day 32 | Day 60 | Meaning |
|:---:|:---:|:---:|:---|
| 0 | not pregnant | not pregnant | not pregnant at either diagnosis |
| 1 | not pregnant | pregnant | pregnancy confirmed only at day 60 |
| 2 | pregnant | not pregnant | observed pregnancy loss between day 32 and day 60 |
| 3 | pregnant | pregnant | pregnant at both diagnoses |

Status 1 has no observations in this data set, so `status` has three observed levels: 0, 2 and 3. `PARITY` (lactation number) has three levels: 1, 2 and 3. Both are treated as fixed classification effects.

---

## 2. Working model

```math
\begin{pmatrix} \mathbf{y}_1 \\ \mathbf{y}_2 \end{pmatrix} =
\mathbf{X}\mathbf{b} +
\begin{pmatrix} \mathbf{W}_1 & \mathbf{0} \\ \mathbf{0} & \mathbf{W}_2 \end{pmatrix}
\begin{pmatrix} \mathbf{a}_1 \\ \mathbf{a}_2 \end{pmatrix} +
\mathbf{e}
\qquad \Longleftrightarrow \qquad
\mathbf{y} = \mathbf{Xb} + \mathbf{Wa} + \mathbf{e}
```

where $`\mathbf{y}_1 \in \mathbb{R}^{n_{y_1} \times 1}`$ holds the phenotypes of cows with pedigree but without genotype, and $`\mathbf{y}_2 \in \mathbb{R}^{n_{y_2} \times 1}`$ the phenotypes of cows with pedigree and genotype. The random effects are independent, with

```math
\mathbf{a} \sim N\left(\mathbf{0}, \mathbf{H}\sigma^2_a\right), \qquad
\mathbf{e} \sim N\left(\mathbf{0}, \mathbf{I}\sigma^2_e\right), \qquad
\text{Cov}\left[\mathbf{a}, \mathbf{e}^\top\right] = \mathbf{0}
```

where $`\mathbf{H}`$ combines pedigree and genomic relationships ([Section 8](docs/genomic-information-and-mme.md#8-genomic-information)).

### 2.1 Fixed effects

$`\mathbf{X} \in \mathbb{R}^{N \times p}`$ and $`\mathbf{b} \in \mathbb{R}^{p \times 1}`$. The fixed effects are an intercept $`\mu`$ and the factors `status` and `PARITY`:

```math
\mathbf{Xb} =
\begin{pmatrix} \mathbf{1}_N & \mathbf{X}_1 & \mathbf{X}_2 \end{pmatrix}
\begin{pmatrix} \mu \\ \mathbf{b}_1 \\ \mathbf{b}_2 \end{pmatrix} =
\mathbf{1}_N \mu + \mathbf{X}_1 \mathbf{b}_1 + \mathbf{X}_2 \mathbf{b}_2,
\qquad p = 1 + n_{b_1} + n_{b_2}
```

where $`\mathbf{b}_1`$ holds the effects of the status levels and $`\mathbf{b}_2`$ the effects of the parity levels.

**Why a constraint is needed.** With an intercept plus one column per level of each factor, the level columns of each factor add up to $`\mathbf{1}_N`$. Then $`\mathbf{X}`$ would not have full column rank, and $`\mathbf{X}^\top\mathbf{X}`$, the fixed-effects block of the mixed model equations (MME), would be singular. ($`\mathbf{X}`$ is rectangular and is never inverted. What must be non-singular is $`\mathbf{X}^\top\mathbf{X}`$ and, with it, the coefficient matrix of the MME.) One constraint per factor is therefore imposed.

**Sum-to-zero constraint.** The effects of the levels of each factor add up to zero. The intercept $`\mu`$ is then the unweighted mean of the level means, and each level effect is a deviation from $`\mu`$. ($`\mu`$ equals the raw mean of the data only with balanced data, which is not the case here.) A factor with $`c`$ levels needs $`c - 1`$ columns in $`\mathbf{X}`$.

**Status** (levels 0, 2, 3; $`n_{b_1} = 2`$). Level 2 (observed loss) is the absorbed level, chosen because it has few observations:

```math
\mathbf{b}_1 = \begin{pmatrix} b_{1_0} \\ b_{1_3} \end{pmatrix},
\qquad
x_{i,j_1=1} =
\begin{cases}
1 & \text{status } 0 \\
-1 & \text{status } 2 \\
0 & \text{status } 3
\end{cases}
\qquad
x_{i,j_1=2} =
\begin{cases}
0 & \text{status } 0 \\
-1 & \text{status } 2 \\
1 & \text{status } 3
\end{cases}
```

**Parity** (levels 1, 2, 3; $`n_{b_2} = 2`$). Level 3 is the absorbed level:

```math
\mathbf{b}_2 = \begin{pmatrix} b_{2_1} \\ b_{2_2} \end{pmatrix},
\qquad
x_{i,j_2=1} =
\begin{cases}
1 & \text{parity } 1 \\
0 & \text{parity } 2 \\
-1 & \text{parity } 3
\end{cases}
\qquad
x_{i,j_2=2} =
\begin{cases}
0 & \text{parity } 1 \\
1 & \text{parity } 2 \\
-1 & \text{parity } 3
\end{cases}
```

Here $`x_{i,j_1}`$ is the element in row $`i \in \lbrace 1, \dots, N \rbrace`$ and column $`j_1`$ of $`\mathbf{X}_1`$ (and likewise for $`\mathbf{X}_2`$). Then $`p = 1 + 2 + 2 = 5`$ and

```math
\mathbf{X} = \begin{pmatrix} \mathbf{1}_N & \mathbf{X}_1 & \mathbf{X}_2 \end{pmatrix} \in \mathbb{R}^{N \times 5},
\qquad
\mathbf{b} = \begin{pmatrix} \mu \\ b_{1_0} \\ b_{1_3} \\ b_{2_1} \\ b_{2_2} \end{pmatrix}
```

| Column of $`\mathbf{X}`$ | Name in the code | Effect |
|:---:|:---|:---|
| 1 | `mu` | $`\mu`$ |
| 2 | `status0` | $`b_{1_0}`$ |
| 3 | `status3` | $`b_{1_3}`$ |
| 4 | `parity1` | $`b_{2_1}`$ |
| 5 | `parity2` | $`b_{2_2}`$ |

The code checks that $`\text{rank}(\mathbf{X}) = p`$, so $`\mathbf{X}^\top\mathbf{X}`$ is invertible. The rows of $`\mathbf{X}`$ are named with the animal ID (one record per animal).

**Recovering the absorbed levels.** The constraints give

```math
\hat{b}_{1_2} = -\hat{b}_{1_0} - \hat{b}_{1_3},
\qquad
\hat{b}_{2_3} = -\hat{b}_{2_1} - \hat{b}_{2_2}
```

and, as linear combinations of estimated parameters,

```math
\text{Var}(\hat{b}_{1_2}) = \text{Var}(\hat{b}_{1_0}) + \text{Var}(\hat{b}_{1_3}) + 2\,\text{Cov}(\hat{b}_{1_0}, \hat{b}_{1_3}),
\qquad
\text{Var}(\hat{b}_{2_3}) = \text{Var}(\hat{b}_{2_1}) + \text{Var}(\hat{b}_{2_2}) + 2\,\text{Cov}(\hat{b}_{2_1}, \hat{b}_{2_2})
```

with the variances and covariances taken from the inverse of the MME coefficient matrix. The choice of the absorbed level does not change $`\hat{\mu}`$ or the level effects.

### 2.2 Random effects

```math
\mathbf{a} = \begin{pmatrix} \mathbf{a}_1 \\ \mathbf{a}_2 \end{pmatrix} \in \mathbb{R}^{K \times 1},
\qquad K = n_{a_1} + n_{a_2}
```

Each level of $`\mathbf{a}`$ is the additive genetic effect of one animal. Every animal belongs to one of four sets. The three digits of the set name indicate, in this order, pedigree, phenotype and genotype (1 = yes, 0 = no):

| Set | Pedigree | Phenotype | Genotype | Belongs to | Record in $`\mathbf{y}`$ |
|:---:|:---:|:---:|:---:|:---:|:---:|
| $`A_{100}`$ | Yes | No  | No  | $`\mathbf{a}_1`$ | No |
| $`A_{110}`$ | Yes | Yes | No  | $`\mathbf{a}_1`$ | Yes ($`\mathbf{y}_1`$) |
| $`A_{101}`$ | Yes | No  | Yes | $`\mathbf{a}_2`$ | No |
| $`A_{111}`$ | Yes | Yes | Yes | $`\mathbf{a}_2`$ | Yes ($`\mathbf{y}_2`$) |

- "Genotyped" means genotyped **and kept after quality control**. An animal removed by quality control is treated as non-genotyped.
- Animals in $`A_{100}`$ and $`A_{101}`$ have all-zero columns in $`\mathbf{W}`$. They stay in $`\mathbf{a}`$ because they connect the phenotyped animals through the relationship matrix.
- All animals are coded as having pedigree. An animal without pedigree is assumed unrelated to the others and is added to the pedigree as a founder (unknown sire and dam).

The pedigree relationship matrix is therefore the only matrix that connects all animals.

---

## 3. Genotype data

The genotypes are processed **before** the pedigree. Quality control can remove animals (low call rate), and a removed animal is treated as non-genotyped. The genotyped animals used to extend the pedigree and to define the sets must therefore be the ones that survive quality control.

### 3.1 From M to Z

$`\mathbf{M} \in \mathbb{R}^{K_g \times J}`$ holds the raw dosages $`m_{ij} \in \lbrace 0, 1, 2 \rbrace`$: the number of copies of the alternative allele of marker $`j \in \lbrace 1, \dots, J \rbrace`$ in animal $`i \in \lbrace 1, \dots, K_g \rbrace`$. Column $`j`$ is $`\mathbf{m}_j = (m_{1j}, \dots, m_{K_g j})^\top`$.

Under Hardy-Weinberg equilibrium (HWE), with allele frequency $`p_j`$, $`\mathbb{P}(m_{ij}=0) = (1-p_j)^2`$, $`\mathbb{P}(m_{ij}=1) = 2(1-p_j)p_j`$ and $`\mathbb{P}(m_{ij}=2) = p_j^2`$. Then

```math
\begin{aligned}
\mathbb{E}\left[m_{ij}\right]
&= 0 \times \mathbb{P}(m_{ij}=0) + 1 \times \mathbb{P}(m_{ij}=1) + 2 \times \mathbb{P}(m_{ij}=2) \\
&= 0\,(1-p_j)^2 + 1\left[2(1-p_j)p_j\right] + 2p_j^2 \\
&= 2p_j - 2p_j^2 + 2p_j^2 \\
&= 2p_j
\end{aligned}
```

Each column is centered by its expected dosage, so every column of $`\mathbf{Z}`$ has mean zero:

```math
\mathbf{z}_j = \mathbf{m}_j - 2p_j,
\qquad
\mathbf{Z} = \begin{pmatrix} \mathbf{z}_1 & \mathbf{z}_2 & \cdots & \mathbf{z}_J \end{pmatrix} \in \mathbb{R}^{K_g \times J}
```

### 3.2 Reading M

- The genotype file has, per line, an animal ID, two numeric fields, and one string of digits (one per SNP). Only the ID and the string are read (`data.table::fread`, as text).
- Genotype IDs are recoded to pedigree IDs with `Keys_to_ID.txt`. IDs without a key keep their original ID.
- Each string is converted to integers with `utf8ToInt(x) - 48` (48 is the character code of `'0'`). Rows of $`\mathbf{M}`$ are animals (pedigree IDs) and columns are SNPs in file order.

### 3.3 Quality control

Besides 0, 1 and 2, the file contains codes 3, 4 and 5 (one or both alleles unknown). These are collapsed into code 6 (unknown). The steps are applied in this order:

| Step | Rule |
|:---|:---|
| Minor allele frequency | $`p_j`$ from known genotypes only, $`p_j = (n_{1j} + 2n_{2j}) / \left[2(n_{0j} + n_{1j} + n_{2j})\right]`$, where $`n_{kj}`$ is the number of animals with dosage $`k`$ at marker $`j`$. $`\text{MAF}_j = \min(p_j, 1-p_j)`$. Keep markers with $`\text{MAF}_j > 0.05`$; also remove markers with no known genotype |
| Call rate per animal | Proportion of the animal's markers with a known genotype. Keep animals with call rate $`\ge 0.90`$. Removed animals are non-genotyped in the rest of the analysis |
| Hardy-Weinberg equilibrium | Exact test (Wigginton et al., 2005; `HardyWeinberg::HWExactStats`) on known genotypes. Bonferroni threshold $`\alpha_{\text{HWE}} = 0.05/J`$, with $`J`$ the number of markers tested. Remove markers with $`p`$-value $`\le \alpha_{\text{HWE}}`$ |
| Imputation | Recompute $`p_j`$ from known genotypes of the surviving animals and markers. Replace each missing dosage of marker $`j`$ by its expected value $`2p_j`$ |

After imputation, every column mean of $`\mathbf{M}`$ equals $`2p_j`$ (checked). Then $`\mathbf{Z}`$ is built by subtracting $`2p_j`$ from column $`j`$ (checked: all column means of $`\mathbf{Z}`$ are 0). The same $`p_j`$ are used later in $`\mathbf{G}`$.

---

## 4. Pedigree

The pedigree file follows CDCB Format 1: 99 characters per line, no header, no delimiter. Only these fields are read:

| Characters | Field | Description |
|:---:|:---|:---|
| 3–19 | `animal` | Animal ID (17 characters) |
| 20–36 | `sire` | Sire ID |
| 37–53 | `dam` | Dam ID |
| 71–78 | `birth` | Birth date `YYYYMMDD`, kept as text (some dates end with `E`) |

Processing steps:

1. **Missing parents.** Blank parents, or parent IDs that do not have 17 characters (e.g., dams coded only as `HO`), are set to `'0'` (unknown).
2. **Parents not listed as animals.** Known sires and dams that do not appear in the `animal` column are added as founders (both parents `'0'`) at the top of the pedigree. The code first checks that no ID is both a sire and a dam.
3. **Genotyped animals not in the pedigree** (after quality control, with pedigree IDs) are added as founders.
4. **Phenotyped cows not in the pedigree** are added as founders, following the assumption in [2.2](#22-random-effects).
5. **Final checks.** Every genotyped animal, phenotyped cow and known parent is an animal of the pedigree, and no animal is repeated.

---

## 5. Numerator relationship matrix A

$`\mathbf{A} \in \mathbb{R}^{K \times K}`$ contains the additive genetic relationships between all animals. It is built with the tabular method (Henderson, 1976; Mrode, 2014, Chapter 2).

**Sorting.** The tabular method needs every parent before its offspring. Each animal gets a generation number: founders are generation 0, and any other animal is 1 + the largest generation of its known parents. Sorting by generation puts parents before offspring. If, in one pass, no new animal can be assigned a generation, the pedigree has a loop (an animal is its own ancestor), and the code stops. Animals are then coded as integers by their row in the sorted pedigree (sire and dam codes are 0 if unknown), and the code checks that every parent has a smaller code than its offspring.

**Tabular rules.** Animal $`i`$ with sire $`s`$ and dam $`d`$ is processed row by row, for all $`j < i`$:

| Parents of $`i`$ | Off-diagonal ($`j < i`$) | Diagonal |
|:---|:---|:---|
| Both known | $`a_{ij} = a_{ji} = 0.5\,(a_{js} + a_{jd})`$ | $`a_{ii} = 1 + 0.5\,a_{sd}`$ |
| One known, parent $`p`$ | $`a_{ij} = a_{ji} = 0.5\,a_{jp}`$ | $`a_{ii} = 1`$ |
| Both unknown | $`a_{ij} = a_{ji} = 0`$ | $`a_{ii} = 1`$ |

The inbreeding coefficient of animal $`i`$ is $`F_i = 0.5\,a_{sd}`$ (0 if a parent is unknown), so $`a_{ii} = 1 + F_i`$. In `make_A()`, row $`i`$ is computed at once for all $`j < i`$.

**Checks.**

1. `make_A()` reproduces the exact $`\mathbf{A}`$ of the 6-animal pedigree of Mrode (2014, Chapter 2).
2. $`\mathbf{A}`$ is symmetric.
3. Founders have $`a_{ii} = 1`$ and are unrelated to each other, so their block of $`\mathbf{A}`$ is an identity matrix.
4. $`a_{ii} = 1 + F_i`$ for every animal.
5. $`\mathbf{A}`$ is identical to the one from an independent implementation, `nadiv::makeA()`.

$`\mathbf{A}`$ is stored dense, so memory grows with $`K^2`$. That is fine for this pedigree. For much larger pedigrees, $`\mathbf{A}^{-1}`$ is built directly, because it is sparse.

---

## 6. Reorder and partition A by sets of animals

To keep $`\mathbf{A}`$, $`\mathbf{G}`$ and $`\mathbf{y}`$ aligned, the rows and columns of $`\mathbf{A}`$ are reordered by set, in the order $`A_{100}, A_{110}, A_{101}, A_{111}`$. Non-genotyped animals ($`\mathbf{a}_1`$) come first and genotyped animals ($`\mathbf{a}_2`$) last. Within each group, phenotyped animals come last, so $`\mathbf{y}_1`$ holds the phenotypes of $`A_{110}`$ and $`\mathbf{y}_2`$ those of $`A_{111}`$. Non-genotyped animals keep the order of the sorted pedigree. Genotyped animals keep the order of the genotype file, which is also the row order of $`\mathbf{M}`$ and $`\mathbf{Z}`$. The reordered matrix is $`\mathbf{A}_{\text{ord}}`$ (from here on, $`\mathbf{A}`$ always refers to this order):

```math
\mathbf{A} =
\begin{pmatrix}
\mathbf{A}_{11} & \mathbf{A}_{12} \\
\mathbf{A}_{21} & \mathbf{A}_{22}
\end{pmatrix}
=
\left(
\begin{array}{cc|cc}
\mathbf{A}_{11_{00}} & \mathbf{A}_{11_{01}} & \mathbf{A}_{12_{00}} & \mathbf{A}_{12_{01}} \\
\mathbf{A}_{11_{10}} & \mathbf{A}_{11_{11}} & \mathbf{A}_{12_{10}} & \mathbf{A}_{12_{11}} \\ \hline
\mathbf{A}_{21_{00}} & \mathbf{A}_{21_{01}} & \mathbf{A}_{22_{00}} & \mathbf{A}_{22_{01}} \\
\mathbf{A}_{21_{10}} & \mathbf{A}_{21_{11}} & \mathbf{A}_{22_{10}} & \mathbf{A}_{22_{11}}
\end{array}
\right)
```

Here $`\mathbf{A}_{11} \in \mathbb{R}^{n_{a_1} \times n_{a_1}}`$ holds the relationships among non-genotyped animals, $`\mathbf{A}_{22} \in \mathbb{R}^{n_{a_2} \times n_{a_2}}`$ those among genotyped animals, and $`\mathbf{A}_{12} = \mathbf{A}_{21}^\top`$ those between the two groups. The four rows (and columns) of sub-blocks correspond to $`A_{100}, A_{110}, A_{101}, A_{111}`$.

**Notation.** In $`\mathbf{A}_{kl_{rc}}`$:

- $`k`$ (rows) and $`l`$ (columns) give the block: 1 = without genotype, 2 = with genotype.
- $`r`$ (rows) and $`c`$ (columns) give the sub-block within that block: 0 = without phenotype, 1 = with phenotype.

For example, $`\mathbf{A}_{12_{01}}`$ has rows from $`A_{100}`$ and columns from $`A_{111}`$. "Between set $`X`$ and set $`Y`$" always means rows from $`X`$ and columns from $`Y`$.

| Sub-block | Rows | Columns | Sub-block | Rows | Columns |
|:---:|:---:|:---:|:---:|:---:|:---:|
| $`\mathbf{A}_{11_{00}}`$ | $`A_{100}`$ | $`A_{100}`$ | $`\mathbf{A}_{12_{00}}`$ | $`A_{100}`$ | $`A_{101}`$ |
| $`\mathbf{A}_{11_{01}}`$ | $`A_{100}`$ | $`A_{110}`$ | $`\mathbf{A}_{12_{01}}`$ | $`A_{100}`$ | $`A_{111}`$ |
| $`\mathbf{A}_{11_{10}}`$ | $`A_{110}`$ | $`A_{100}`$ | $`\mathbf{A}_{12_{10}}`$ | $`A_{110}`$ | $`A_{101}`$ |
| $`\mathbf{A}_{11_{11}}`$ | $`A_{110}`$ | $`A_{110}`$ | $`\mathbf{A}_{12_{11}}`$ | $`A_{110}`$ | $`A_{111}`$ |
| $`\mathbf{A}_{22_{00}}`$ | $`A_{101}`$ | $`A_{101}`$ | $`\mathbf{A}_{21_{00}}`$ | $`A_{101}`$ | $`A_{100}`$ |
| $`\mathbf{A}_{22_{01}}`$ | $`A_{101}`$ | $`A_{111}`$ | $`\mathbf{A}_{21_{01}}`$ | $`A_{101}`$ | $`A_{110}`$ |
| $`\mathbf{A}_{22_{10}}`$ | $`A_{111}`$ | $`A_{101}`$ | $`\mathbf{A}_{21_{10}}`$ | $`A_{111}`$ | $`A_{100}`$ |
| $`\mathbf{A}_{22_{11}}`$ | $`A_{111}`$ | $`A_{111}`$ | $`\mathbf{A}_{21_{11}}`$ | $`A_{111}`$ | $`A_{110}`$ |

**Symmetry.** Only ten of the sixteen sub-blocks are distinct. The other six are transposes:

```math
\mathbf{A}_{11_{10}} = \mathbf{A}_{11_{01}}^\top, \qquad
\mathbf{A}_{22_{10}} = \mathbf{A}_{22_{01}}^\top
```

```math
\mathbf{A}_{21} = \mathbf{A}_{12}^\top =
\begin{pmatrix}
\mathbf{A}_{12_{00}}^\top & \mathbf{A}_{12_{10}}^\top \\
\mathbf{A}_{12_{01}}^\top & \mathbf{A}_{12_{11}}^\top
\end{pmatrix}
\quad \Rightarrow \quad
\mathbf{A}_{21_{00}} = \mathbf{A}_{12_{00}}^\top, \quad
\mathbf{A}_{21_{01}} = \mathbf{A}_{12_{10}}^\top, \quad
\mathbf{A}_{21_{10}} = \mathbf{A}_{12_{01}}^\top, \quad
\mathbf{A}_{21_{11}} = \mathbf{A}_{12_{11}}^\top
```

Transposing $`\mathbf{A}_{12}`$ swaps its two off-diagonal sub-blocks. The code checks all six relations.

**Empty sets.** An empty set gives sub-blocks with zero rows or columns, which drop out (`drop = FALSE` keeps them as matrices). In this data set $`A_{101}`$ is empty (all genotyped animals have a phenotype), so

```math
\mathbf{A}_{22} = \mathbf{A}_{22_{11}}, \qquad
\mathbf{A}_{12} = \begin{pmatrix} \mathbf{A}_{12_{01}} \\ \mathbf{A}_{12_{11}} \end{pmatrix}, \qquad
\mathbf{A}_{21} = \begin{pmatrix} \mathbf{A}_{21_{10}} & \mathbf{A}_{21_{11}} \end{pmatrix}
```

The 16 sub-blocks are bound into the 4 blocks, and the 4 blocks into $`\mathbf{A}_{\text{ord}}`$. The code checks that each block equals $`\mathbf{A}`$ indexed by the IDs of $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$, and that $`\mathbf{A}_{\text{ord}}`$ has the same animals and relationships as $`\mathbf{A}`$.

---

## 7. Align data with the ordered A

- $`\mathbf{y} = (\mathbf{y}_1^\top, \mathbf{y}_2^\top)^\top`$ follows the order of the phenotyped animals in $`\mathbf{A}_{\text{ord}}`$: first $`A_{110}`$, then $`A_{111}`$. Each value is taken by animal ID.
- The rows of $`\mathbf{X}`$ are reordered ($`\mathbf{X}_{\text{ord}}`$) so that row $`i`$ of $`\mathbf{X}_{\text{ord}}`$ and element $`i`$ of $`\mathbf{y}`$ belong to the same animal. The code checks that each row keeps its animal's status and parity coding.
- The rows of $`\mathbf{Z}`$ are reordered as $`\mathbf{a}_2`$ (first $`A_{101}`$, then $`A_{111}`$), so $`\mathbf{Z}`$, and later $`\mathbf{G}`$, follow the genotyped animals in $`\mathbf{A}_{\text{ord}}`$. With $`A_{101}`$ empty, the order does not change, but the step keeps the code correct in general.

---

Sections 8–9 continue in [Genomic information and mixed model equations](docs/genomic-information-and-mme.md), and sections 10–12 in [Restricted maximum likelihood and AI-REML](docs/reml.md).

---

## References

- Aguilar, I., Misztal, I., Johnson, D. L., Legarra, A., Tsuruta, S., and Lawlor, T. J. (2010). Hot topic: A unified approach to utilize phenotypic, full pedigree, and genomic information for genetic evaluation of Holstein final score. *Journal of Dairy Science*, 93, 743–752.
- Dempster, A. P., Laird, N. M., and Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society, Series B*, 39, 1–38.
- Gilmour, A. R., Thompson, R., and Cullis, B. R. (1995). Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models. *Biometrics*, 51, 1440–1450.
- Henderson, C. R. (1976). A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values. *Biometrics*, 32, 69–83.
- Henderson, C. R. (1984). *Applications of Linear Models in Animal Breeding*. University of Guelph, Guelph, ON.
- Legarra, A., Aguilar, I., and Misztal, I. (2009). A relationship matrix including full pedigree and genomic information. *Journal of Dairy Science*, 92, 4656–4663.
- Meyer, K. (1989). Restricted maximum likelihood to estimate variance components for animal models with several random effects using a derivative-free algorithm. *Genetics Selection Evolution*, 21, 317–340.
- Mrode, R. A. (2014). *Linear Models for the Prediction of Animal Breeding Values* (3rd ed.). CABI, Wallingford, UK.
- Patterson, H. D., and Thompson, R. (1971). Recovery of inter-block information when block sizes are unequal. *Biometrika*, 58, 545–554.
- VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. *Journal of Dairy Science*, 91, 4414–4423.
- Wigginton, J. E., Cutler, D. J., and Abecasis, G. R. (2005). A note on exact tests of Hardy-Weinberg equilibrium. *American Journal of Human Genetics*, 76, 887–893.
