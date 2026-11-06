# PAGssGBLUP

Single-step genomic BLUP (ssGBLUP) analysis of the change in pregnancy-associated glycoproteins (PAGs) between day 16 and day 26 in dairy cows, with variance components estimated by average-information REML (AI-REML).

The full analysis is in [`code/PAGssGBLUP.Rmd`](code/PAGssGBLUP.Rmd). This README describes the theory and the methodology applied there, step by step and in the same notation. Section names in *italics* refer to sections of the Rmd.

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
11. [Genomic information](#8-genomic-information)
12. [Mixed model equations](#9-mixed-model-equations)
13. [Restricted maximum likelihood](#10-restricted-maximum-likelihood)
14. [AI-REML implementation](#11-ai-reml-implementation)
15. [Conclusion](#12-conclusion)
16. [References](#references)

---

## Repository structure

```
PAGssGBLUP/
├── README.md
├── .gitignore          # excludes data/, rendered HTML and R session files
└── code/
    └── PAGssGBLUP.Rmd  # full analysis
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

Required R packages: `openxlsx`, `readr`, `data.table`, `HardyWeinberg`, `nadiv` (used only as an independent check of $\mathbf{A}$), and `rmarkdown`/`knitr` to knit. Knit the Rmd from RStudio, or run its chunks with `code/` as the working directory.

---

## Notation

Matrices are bold upper case, vectors bold lower case, scalars in italics. $^\top$ is the transpose.

| Symbol | Meaning |
|:---|:---|
| $N = n_{y_1} + n_{y_2}$ | Number of phenotypes (one per cow), split in non-genotyped ($n_{y_1}$) and genotyped ($n_{y_2}$) cows |
| $K = n_{a_1} + n_{a_2}$ | Number of animals in the pedigree (levels of $\mathbf{a}$), split in non-genotyped ($n_{a_1}$) and genotyped ($n_{a_2}$) |
| $K_g$ | Number of genotyped animals after quality control; $K_g = n_{a_2}$ |
| $J$ | Number of markers (SNPs) after quality control |
| $p$ | Number of fixed-effect parameters (columns of $\mathbf{X}$); `n_b` in the code |
| $p_j$, $q_j = 1 - p_j$ | Frequency of the alternative allele of marker $j$, and of the other allele; `p` in the code |
| $\mathbf{y}$, $\mathbf{b}$, $\mathbf{a}$, $\mathbf{e}$ | Phenotypes, fixed effects, additive genetic effects, residuals |
| $\mathbf{X}$, $\mathbf{W}$ | Design matrices of $\mathbf{b}$ and $\mathbf{a}$ |
| $\mathbf{M}$, $\mathbf{Z}$ | Raw SNP dosages, and dosages centered by $2p_j$ |
| $\mathbf{g}$ | Additive effects of the markers |
| $\mathbf{A}$, $\mathbf{G}$, $\mathbf{H}$ | Pedigree, genomic, and combined (pedigree + genomic) relationship matrices |
| $\mathbf{G}_w$ | Blended genomic matrix, $(1-w)\mathbf{G} + w\mathbf{A}_{22}$ |
| $\sigma^2_a$, $\sigma^2_e$, $\sigma^2_g$ | Additive genetic, residual, and per-marker variance |
| $\alpha = \sigma^2_e/\sigma^2_a$ | Variance ratio in the MME |
| $h^2 = \sigma^2_a/(\sigma^2_a + \sigma^2_e)$ | Heritability |
| $\boldsymbol{\theta} = (\sigma^2_a, \sigma^2_e)^\top$ | Variance components estimated by REML |
| $\mathbf{V}$, $\mathbf{P}$ | Variance of $\mathbf{y}$, and REML projection matrix |
| $L$ | REML log-likelihood |
| $A_{100}, A_{110}, A_{101}, A_{111}$ | Sets of animals (pedigree, phenotype, genotype; 1 = yes, 0 = no) |

Subscript 1 always refers to animals **without** genotype, and subscript 2 to animals **with** genotype.

Three symbols are used in two places, and the context separates them. $p$ (number of columns of $\mathbf{X}$) is not $p_j$ (allele frequency). $F_i$ (inbreeding coefficient) is not $\mathbf{F}$ (the $N \times 2$ matrix in the AI matrix). $\mathbf{P}$ (projection matrix) is not $\mathbb{P}$ (probability). The Rmd also uses $\alpha$ for the Hardy-Weinberg threshold and $\mathbf{g}$ for the gradient of $h^2$; here these are written $\alpha_{\text{HWE}}$ and $\nabla h^2$ to avoid confusion with the variance ratio and the marker effects.

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

where $\mathbf{y}_1 \in \mathbb{R}^{n_{y_1} \times 1}$ holds the phenotypes of cows with pedigree but without genotype, and $\mathbf{y}_2 \in \mathbb{R}^{n_{y_2} \times 1}$ the phenotypes of cows with pedigree and genotype. The random effects are independent, with

```math
\mathbf{a} \sim N\left(\mathbf{0}, \mathbf{H}\sigma^2_a\right), \qquad
\mathbf{e} \sim N\left(\mathbf{0}, \mathbf{I}\sigma^2_e\right), \qquad
\text{Cov}\left[\mathbf{a}, \mathbf{e}^\top\right] = \mathbf{0}
```

where $\mathbf{H}$ combines pedigree and genomic relationships ([Section 8](#8-genomic-information)).

### 2.1 Fixed effects

$\mathbf{X} \in \mathbb{R}^{N \times p}$ and $\mathbf{b} \in \mathbb{R}^{p \times 1}$. The fixed effects are an intercept $\mu$ and the factors `status` and `PARITY`:

```math
\mathbf{Xb} =
\begin{pmatrix} \mathbf{1}_N & \mathbf{X}_1 & \mathbf{X}_2 \end{pmatrix}
\begin{pmatrix} \mu \\ \mathbf{b}_1 \\ \mathbf{b}_2 \end{pmatrix} =
\mathbf{1}_N \mu + \mathbf{X}_1 \mathbf{b}_1 + \mathbf{X}_2 \mathbf{b}_2,
\qquad p = 1 + n_{b_1} + n_{b_2}
```

where $\mathbf{b}_1$ holds the effects of the status levels and $\mathbf{b}_2$ the effects of the parity levels.

**Why a constraint is needed.** With an intercept plus one column per level of each factor, the level columns of each factor add up to $\mathbf{1}_N$. Then $\mathbf{X}$ would not have full column rank, and $\mathbf{X}^\top\mathbf{X}$, the fixed-effects block of the mixed model equations (MME), would be singular. ($\mathbf{X}$ is rectangular and is never inverted. What must be non-singular is $\mathbf{X}^\top\mathbf{X}$ and, with it, the coefficient matrix of the MME.) One constraint per factor is therefore imposed.

**Sum-to-zero constraint.** The effects of the levels of each factor add up to zero. The intercept $\mu$ is then the unweighted mean of the level means, and each level effect is a deviation from $\mu$. ($\mu$ equals the raw mean of the data only with balanced data, which is not the case here.) A factor with $c$ levels needs $c - 1$ columns in $\mathbf{X}$.

**Status** (levels 0, 2, 3; $n_{b_1} = 2$). Level 2 (observed loss) is the absorbed level, chosen because it has few observations:

```math
\mathbf{b}_1 = \begin{pmatrix} b_{1_0} \\ b_{1_3} \end{pmatrix},
\qquad
x_{i,j_1=1} =
\begin{cases}
\phantom{-}1 & \text{status } 0 \\
-1 & \text{status } 2 \\
\phantom{-}0 & \text{status } 3
\end{cases}
\qquad
x_{i,j_1=2} =
\begin{cases}
\phantom{-}0 & \text{status } 0 \\
-1 & \text{status } 2 \\
\phantom{-}1 & \text{status } 3
\end{cases}
```

**Parity** (levels 1, 2, 3; $n_{b_2} = 2$). Level 3 is the absorbed level:

```math
\mathbf{b}_2 = \begin{pmatrix} b_{2_1} \\ b_{2_2} \end{pmatrix},
\qquad
x_{i,j_2=1} =
\begin{cases}
\phantom{-}1 & \text{parity } 1 \\
\phantom{-}0 & \text{parity } 2 \\
-1 & \text{parity } 3
\end{cases}
\qquad
x_{i,j_2=2} =
\begin{cases}
\phantom{-}0 & \text{parity } 1 \\
\phantom{-}1 & \text{parity } 2 \\
-1 & \text{parity } 3
\end{cases}
```

Here $x_{i,j_1}$ is the element in row $i \in \lbrace 1, \dots, N \rbrace$ and column $j_1$ of $\mathbf{X}_1$ (and likewise for $\mathbf{X}_2$). Then $p = 1 + 2 + 2 = 5$ and

```math
\mathbf{X} = \begin{pmatrix} \mathbf{1}_N & \mathbf{X}_1 & \mathbf{X}_2 \end{pmatrix} \in \mathbb{R}^{N \times 5},
\qquad
\mathbf{b} = \begin{pmatrix} \mu \\ b_{1_0} \\ b_{1_3} \\ b_{2_1} \\ b_{2_2} \end{pmatrix}
```

| Column of $\mathbf{X}$ | Name in the code | Effect |
|:---:|:---|:---|
| 1 | `mu` | $\mu$ |
| 2 | `status0` | $b_{1_0}$ |
| 3 | `status3` | $b_{1_3}$ |
| 4 | `parity1` | $b_{2_1}$ |
| 5 | `parity2` | $b_{2_2}$ |

The code checks that $\text{rank}(\mathbf{X}) = p$, so $\mathbf{X}^\top\mathbf{X}$ is invertible. The rows of $\mathbf{X}$ are named with the animal ID (one record per animal).

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

with the variances and covariances taken from the inverse of the MME coefficient matrix. The choice of the absorbed level does not change $\hat{\mu}$ or the level effects.

### 2.2 Random effects

```math
\mathbf{a} = \begin{pmatrix} \mathbf{a}_1 \\ \mathbf{a}_2 \end{pmatrix} \in \mathbb{R}^{K \times 1},
\qquad K = n_{a_1} + n_{a_2}
```

Each level of $\mathbf{a}$ is the additive genetic effect of one animal. Every animal belongs to one of four sets. The three digits of the set name indicate, in this order, pedigree, phenotype and genotype (1 = yes, 0 = no):

| Set | Pedigree | Phenotype | Genotype | Belongs to | Record in $\mathbf{y}$ |
|:---:|:---:|:---:|:---:|:---:|:---:|
| $A_{100}$ | Yes | No  | No  | $\mathbf{a}_1$ | No |
| $A_{110}$ | Yes | Yes | No  | $\mathbf{a}_1$ | Yes ($\mathbf{y}_1$) |
| $A_{101}$ | Yes | No  | Yes | $\mathbf{a}_2$ | No |
| $A_{111}$ | Yes | Yes | Yes | $\mathbf{a}_2$ | Yes ($\mathbf{y}_2$) |

- "Genotyped" means genotyped **and kept after quality control**. An animal removed by quality control is treated as non-genotyped.
- Animals in $A_{100}$ and $A_{101}$ have all-zero columns in $\mathbf{W}$. They stay in $\mathbf{a}$ because they connect the phenotyped animals through the relationship matrix.
- All animals are coded as having pedigree. An animal without pedigree is assumed unrelated to the others and is added to the pedigree as a founder (unknown sire and dam).

The pedigree relationship matrix is therefore the only matrix that connects all animals.

---

## 3. Genotype data

The genotypes are processed **before** the pedigree. Quality control can remove animals (low call rate), and a removed animal is treated as non-genotyped. The genotyped animals used to extend the pedigree and to define the sets must therefore be the ones that survive quality control.

### 3.1 From M to Z

$\mathbf{M} \in \mathbb{R}^{K_g \times J}$ holds the raw dosages $m_{ij} \in \lbrace 0, 1, 2 \rbrace$: the number of copies of the alternative allele of marker $j \in \lbrace 1, \dots, J \rbrace$ in animal $i \in \lbrace 1, \dots, K_g \rbrace$. Column $j$ is $\mathbf{m}_j = (m_{1j}, \dots, m_{K_g j})^\top$.

Under Hardy-Weinberg equilibrium (HWE), with allele frequency $p_j$, $\mathbb{P}(m_{ij}=0) = (1-p_j)^2$, $\mathbb{P}(m_{ij}=1) = 2(1-p_j)p_j$ and $\mathbb{P}(m_{ij}=2) = p_j^2$. Then

```math
\begin{aligned}
\mathbb{E}\left[m_{ij}\right]
&= 0 \times \mathbb{P}(m_{ij}=0) + 1 \times \mathbb{P}(m_{ij}=1) + 2 \times \mathbb{P}(m_{ij}=2) \\
&= 0\,(1-p_j)^2 + 1\left[2(1-p_j)p_j\right] + 2p_j^2 \\
&= 2p_j - 2p_j^2 + 2p_j^2 \\
&= 2p_j
\end{aligned}
```

Each column is centered by its expected dosage, so every column of $\mathbf{Z}$ has mean zero:

```math
\mathbf{z}_j = \mathbf{m}_j - 2p_j,
\qquad
\mathbf{Z} = \begin{pmatrix} \mathbf{z}_1 & \mathbf{z}_2 & \cdots & \mathbf{z}_J \end{pmatrix} \in \mathbb{R}^{K_g \times J}
```

### 3.2 Reading M

- The genotype file has, per line, an animal ID, two numeric fields, and one string of digits (one per SNP). Only the ID and the string are read (`data.table::fread`, as text).
- Genotype IDs are recoded to pedigree IDs with `Keys_to_ID.txt`. IDs without a key keep their original ID.
- Each string is converted to integers with `utf8ToInt(x) - 48` (48 is the character code of `'0'`). Rows of $\mathbf{M}$ are animals (pedigree IDs) and columns are SNPs in file order.

### 3.3 Quality control

Besides 0, 1 and 2, the file contains codes 3, 4 and 5 (one or both alleles unknown). These are collapsed into code 6 (unknown). The steps are applied in this order:

| Step | Rule |
|:---|:---|
| Minor allele frequency | $p_j$ from known genotypes only, $p_j = (n_{1j} + 2n_{2j}) / \left[2(n_{0j} + n_{1j} + n_{2j})\right]$, where $n_{kj}$ is the number of animals with dosage $k$ at marker $j$. $\text{MAF}_j = \min(p_j, 1-p_j)$. Keep markers with $\text{MAF}_j > 0.05$; also remove markers with no known genotype |
| Call rate per animal | Proportion of the animal's markers with a known genotype. Keep animals with call rate $\ge 0.90$. Removed animals are non-genotyped in the rest of the analysis |
| Hardy-Weinberg equilibrium | Exact test (Wigginton et al., 2005; `HardyWeinberg::HWExactStats`) on known genotypes. Bonferroni threshold $\alpha_{\text{HWE}} = 0.05/J$, with $J$ the number of markers tested. Remove markers with $p$-value $\le \alpha_{\text{HWE}}$ |
| Imputation | Recompute $p_j$ from known genotypes of the surviving animals and markers. Replace each missing dosage of marker $j$ by its expected value $2p_j$ |

After imputation, every column mean of $\mathbf{M}$ equals $2p_j$ (checked). Then $\mathbf{Z}$ is built by subtracting $2p_j$ from column $j$ (checked: all column means of $\mathbf{Z}$ are 0). The same $p_j$ are used later in $\mathbf{G}$.

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

$\mathbf{A} \in \mathbb{R}^{K \times K}$ contains the additive genetic relationships between all animals. It is built with the tabular method (Henderson, 1976; Mrode, 2014, Chapter 2).

**Sorting.** The tabular method needs every parent before its offspring. Each animal gets a generation number: founders are generation 0, and any other animal is 1 + the largest generation of its known parents. Sorting by generation puts parents before offspring. If, in one pass, no new animal can be assigned a generation, the pedigree has a loop (an animal is its own ancestor), and the code stops. Animals are then coded as integers by their row in the sorted pedigree (sire and dam codes are 0 if unknown), and the code checks that every parent has a smaller code than its offspring.

**Tabular rules.** Animal $i$ with sire $s$ and dam $d$ is processed row by row, for all $j < i$:

| Parents of $i$ | Off-diagonal ($j < i$) | Diagonal |
|:---|:---|:---|
| Both known | $a_{ij} = a_{ji} = 0.5\,(a_{js} + a_{jd})$ | $a_{ii} = 1 + 0.5\,a_{sd}$ |
| One known, parent $p$ | $a_{ij} = a_{ji} = 0.5\,a_{jp}$ | $a_{ii} = 1$ |
| Both unknown | $a_{ij} = a_{ji} = 0$ | $a_{ii} = 1$ |

The inbreeding coefficient of animal $i$ is $F_i = 0.5\,a_{sd}$ (0 if a parent is unknown), so $a_{ii} = 1 + F_i$. In `make_A()`, row $i$ is computed at once for all $j < i$.

**Checks.**

1. `make_A()` reproduces the exact $\mathbf{A}$ of the 6-animal pedigree of Mrode (2014, Chapter 2).
2. $\mathbf{A}$ is symmetric.
3. Founders have $a_{ii} = 1$ and are unrelated to each other, so their block of $\mathbf{A}$ is an identity matrix.
4. $a_{ii} = 1 + F_i$ for every animal.
5. $\mathbf{A}$ is identical to the one from an independent implementation, `nadiv::makeA()`.

$\mathbf{A}$ is stored dense, so memory grows with $K^2$. That is fine for this pedigree. For much larger pedigrees, $\mathbf{A}^{-1}$ is built directly, because it is sparse.

---

## 6. Reorder and partition A by sets of animals

To keep $\mathbf{A}$, $\mathbf{G}$ and $\mathbf{y}$ aligned, the rows and columns of $\mathbf{A}$ are reordered by set, in the order $A_{100}, A_{110}, A_{101}, A_{111}$. Non-genotyped animals ($\mathbf{a}_1$) come first and genotyped animals ($\mathbf{a}_2$) last. Within each group, phenotyped animals come last, so $\mathbf{y}_1$ holds the phenotypes of $A_{110}$ and $\mathbf{y}_2$ those of $A_{111}$. Non-genotyped animals keep the order of the sorted pedigree. Genotyped animals keep the order of the genotype file, which is also the row order of $\mathbf{M}$ and $\mathbf{Z}$. The reordered matrix is $\mathbf{A}_{\text{ord}}$ (from here on, $\mathbf{A}$ always refers to this order):

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
\mathbf{A}_{11_{10}} & \mathbf{A}_{11_{11}} & \mathbf{A}_{12_{10}} & \mathbf{A}_{12_{11}} \\
\hline
\mathbf{A}_{21_{00}} & \mathbf{A}_{21_{01}} & \mathbf{A}_{22_{00}} & \mathbf{A}_{22_{01}} \\
\mathbf{A}_{21_{10}} & \mathbf{A}_{21_{11}} & \mathbf{A}_{22_{10}} & \mathbf{A}_{22_{11}}
\end{array}
\right)
```

Here $\mathbf{A}_{11} \in \mathbb{R}^{n_{a_1} \times n_{a_1}}$ holds the relationships among non-genotyped animals, $\mathbf{A}_{22} \in \mathbb{R}^{n_{a_2} \times n_{a_2}}$ those among genotyped animals, and $\mathbf{A}_{12} = \mathbf{A}_{21}^\top$ those between the two groups. The four rows (and columns) of sub-blocks correspond to $A_{100}, A_{110}, A_{101}, A_{111}$.

**Notation.** In $\mathbf{A}_{kl_{rc}}$:

- $k$ (rows) and $l$ (columns) give the block: 1 = without genotype, 2 = with genotype.
- $r$ (rows) and $c$ (columns) give the sub-block within that block: 0 = without phenotype, 1 = with phenotype.

For example, $\mathbf{A}_{12_{01}}$ has rows from $A_{100}$ and columns from $A_{111}$. "Between set $X$ and set $Y$" always means rows from $X$ and columns from $Y$.

| Sub-block | Rows | Columns | Sub-block | Rows | Columns |
|:---:|:---:|:---:|:---:|:---:|:---:|
| $\mathbf{A}_{11_{00}}$ | $A_{100}$ | $A_{100}$ | $\mathbf{A}_{12_{00}}$ | $A_{100}$ | $A_{101}$ |
| $\mathbf{A}_{11_{01}}$ | $A_{100}$ | $A_{110}$ | $\mathbf{A}_{12_{01}}$ | $A_{100}$ | $A_{111}$ |
| $\mathbf{A}_{11_{10}}$ | $A_{110}$ | $A_{100}$ | $\mathbf{A}_{12_{10}}$ | $A_{110}$ | $A_{101}$ |
| $\mathbf{A}_{11_{11}}$ | $A_{110}$ | $A_{110}$ | $\mathbf{A}_{12_{11}}$ | $A_{110}$ | $A_{111}$ |
| $\mathbf{A}_{22_{00}}$ | $A_{101}$ | $A_{101}$ | $\mathbf{A}_{21_{00}}$ | $A_{101}$ | $A_{100}$ |
| $\mathbf{A}_{22_{01}}$ | $A_{101}$ | $A_{111}$ | $\mathbf{A}_{21_{01}}$ | $A_{101}$ | $A_{110}$ |
| $\mathbf{A}_{22_{10}}$ | $A_{111}$ | $A_{101}$ | $\mathbf{A}_{21_{10}}$ | $A_{111}$ | $A_{100}$ |
| $\mathbf{A}_{22_{11}}$ | $A_{111}$ | $A_{111}$ | $\mathbf{A}_{21_{11}}$ | $A_{111}$ | $A_{110}$ |

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

Transposing $\mathbf{A}_{12}$ swaps its two off-diagonal sub-blocks. The code checks all six relations.

**Empty sets.** An empty set gives sub-blocks with zero rows or columns, which drop out (`drop = FALSE` keeps them as matrices). In this data set $A_{101}$ is empty (all genotyped animals have a phenotype), so

```math
\mathbf{A}_{22} = \mathbf{A}_{22_{11}}, \qquad
\mathbf{A}_{12} = \begin{pmatrix} \mathbf{A}_{12_{01}} \\ \mathbf{A}_{12_{11}} \end{pmatrix}, \qquad
\mathbf{A}_{21} = \begin{pmatrix} \mathbf{A}_{21_{10}} & \mathbf{A}_{21_{11}} \end{pmatrix}
```

The 16 sub-blocks are bound into the 4 blocks, and the 4 blocks into $\mathbf{A}_{\text{ord}}$. The code checks that each block equals $\mathbf{A}$ indexed by the IDs of $\mathbf{a}_1$ and $\mathbf{a}_2$, and that $\mathbf{A}_{\text{ord}}$ has the same animals and relationships as $\mathbf{A}$.

---

## 7. Align data with the ordered A

- $\mathbf{y} = (\mathbf{y}_1^\top, \mathbf{y}_2^\top)^\top$ follows the order of the phenotyped animals in $\mathbf{A}_{\text{ord}}$: first $A_{110}$, then $A_{111}$. Each value is taken by animal ID.
- The rows of $\mathbf{X}$ are reordered ($\mathbf{X}_{\text{ord}}$) so that row $i$ of $\mathbf{X}_{\text{ord}}$ and element $i$ of $\mathbf{y}$ belong to the same animal. The code checks that each row keeps its animal's status and parity coding.
- The rows of $\mathbf{Z}$ are reordered as $\mathbf{a}_2$ (first $A_{101}$, then $A_{111}$), so $\mathbf{Z}$, and later $\mathbf{G}$, follow the genotyped animals in $\mathbf{A}_{\text{ord}}$. With $A_{101}$ empty, the order does not change, but the step keeps the code correct in general.

---

## 8. Genomic information

### 8.1 The pedigree propagates genomic information

The pedigree relationships are the bridge from genotyped to non-genotyped animals. With $\mathbb{E}[\mathbf{a}_1] = \mathbf{0}$ and $\mathbb{E}[\mathbf{a}_2] = \mathbf{0}$:

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

These are the pedigree (co)variances: the link between $\mathbf{a}_1$ and $\mathbf{a}_2$ comes from the pedigree. The genomic information enters through $\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a$ ([8.5](#85-variance-of-a)).

The additive effects of the genotyped animals are expressed through the markers:

```math
\mathbf{a}_2 = \mathbf{Zg}, \qquad \mathbf{Z} \in \mathbb{R}^{K_g \times J}, \qquad \mathbf{g} \in \mathbb{R}^{J \times 1}
```

The rows of $\mathbf{Z}$ are the animals of $\mathbf{a}_2$ ($K_g = n_{a_2}$). Genotyped animals without phenotypes (e.g., genotyped sires) could also be part of $\mathbf{a}_2$, although there are none in this data set.

### 8.2 Marker effects

The marker effects are independent with a common variance $\sigma^2_g$:

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

To link $\sigma^2_g$ with $\sigma^2_a$, the dosage $m_{ij}$ is treated as a random variable: the genotype of a random animal from a population in HWE, with allele frequency $p_j$ and $q_j = 1 - p_j$.

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

The centered dosage $z_{ij} = m_{ij} - 2p_j$ only shifts $m_{ij}$ by a constant, so $\mathbb{E}[z_{ij}] = 0$ and $\text{Var}[z_{ij}] = \mathbb{E}[z_{ij}^2] = 2p_jq_j$.

The additive genetic value of animal $i$ (element $i$ of $\mathbf{Zg}$) is

```math
a_i = \sum_{j=1}^{J} z_{ij}g_j
```

The marker effects are assumed independent of the genotypes. Because the $g_j$ have mean zero and are uncorrelated, terms of different markers are uncorrelated, even if their dosages are correlated (linkage disequilibrium):

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

With $\sigma^2_a = \text{Var}[a_i]$:

```math
\sigma^2_a = 2\sum_{j=1}^{J}p_jq_j\;\sigma^2_g
\quad \Longrightarrow \quad
\sigma^2_g = \frac{\sigma^2_a}{2\sum_{j=1}^{J}p_jq_j}
```

With $\mathbf{Z}$ fixed (the observed genotypes):

```math
\text{Var}[\mathbf{a}_2] = \text{Var}[\mathbf{Zg}]
= \mathbf{Z}\mathbf{Z}^\top\sigma^2_g
= \frac{\mathbf{Z}\mathbf{Z}^\top}{2\sum_{j=1}^{J}p_jq_j}\,\sigma^2_a
= \mathbf{G}\sigma^2_a,
\qquad
\mathbf{G} = \frac{\mathbf{Z}\mathbf{Z}^\top}{2\sum_{j=1}^{J}p_jq_j} \in \mathbb{R}^{K_g \times K_g}
```

This is the genomic relationship matrix of VanRaden (2008). It connects the marker model ($\mathbf{g}$, $\sigma^2_g$) with the animal model ($\mathbf{a}$, $\sigma^2_a$). The $p_j$ should be base-population frequencies. Because those are unknown, the frequencies of the genotyped animals are used (the same $p_j$ used to center $\mathbf{Z}$). For genotyped animals, $\mathbf{G}$ replaces $\mathbf{A}_{22}$.

**Computation.** $\mathbf{Z}\mathbf{Z}^\top$ is computed with `tcrossprod(Z)` (one BLAS call, without forming $\mathbf{Z}^\top$), and divided by $2\sum_j p_jq_j$.

**Checks.**

1. $\mathbf{G}$ is symmetric.
2. Its rows and columns follow $\mathbf{a}_2$, the same order as $\mathbf{A}_{22}$.
3. Every row adds up to zero: each column of $\mathbf{Z}$ has mean zero, so $\mathbf{Z}^\top\mathbf{1} = \mathbf{0}$ and $\mathbf{G}\mathbf{1} = \mathbf{Z}(\mathbf{Z}^\top\mathbf{1})/(2\sum_j p_jq_j) = \mathbf{0}$.

**$\mathbf{G}$ is singular.** Since $\mathbf{G}\mathbf{1} = \mathbf{0}$, $\text{rank}(\mathbf{G}) \le K_g - 1$ and its smallest eigenvalue is zero (up to rounding). $\mathbf{G}^{-1}$ does not exist, but $\mathbf{H}^{-1}$ needs it.

### 8.4 Blending G

```math
\mathbf{G}_w = (1 - w)\,\mathbf{G} + w\,\mathbf{A}_{22}, \qquad w = 0.05
```

($w = 0.05$ is a common default, e.g., in BLUPF90.)

**Meaning of $w$.** The markers may not capture all the additive genetic variance. The additive effect of a genotyped animal is split into a part explained by the markers and a residual polygenic part $\mathbf{u}$, which follows the pedigree and is independent of $\mathbf{g}$:

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

The markers explain a proportion $1 - w$ of $\sigma^2_a$ and the pedigree the remaining $w$. With $w = 0$, $\mathbf{a}_2 = \mathbf{Zg}$ as in [8.3](#83-genomic-relationship-matrix-g).

$\mathbf{G}_w$ is positive definite. Since $\mathbf{G}$ is positive semi-definite, the smallest eigenvalue of $\mathbf{G}_w$ is at least $w$ times the smallest eigenvalue of $\mathbf{A}_{22}$ (checked in the code). Its inverse is computed from the Cholesky factor (`chol2inv(chol(Gw))`). This is faster and more stable than `solve()` for a symmetric positive definite matrix, and `chol()` stops if $\mathbf{G}_w$ is not positive definite. The code checks that $\mathbf{G}_w\mathbf{G}_w^{-1} = \mathbf{I}$.

**Effect of blending on what follows.** The derivation of $\text{Var}[\mathbf{a}]$ below uses only that $\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a$ with $\mathbf{G}$ symmetric, not how $\mathbf{G}$ is built. Every step therefore holds with $\mathbf{G}_w$ in place of $\mathbf{G}$. **From here on, $\mathbf{G}$ in the formulas stands for the matrix actually used, $\mathbf{G}_w$.** Since

```math
\mathbf{G}_w - \mathbf{A}_{22} = (1-w)\,\mathbf{G} + w\,\mathbf{A}_{22} - \mathbf{A}_{22} = (1-w)\left(\mathbf{G} - \mathbf{A}_{22}\right)
```

blending shrinks the genomic correction of the non-genotyped animals by $(1-w)$:

```math
\text{Var}[\mathbf{a}_1] =
\left[\mathbf{A}_{11} + (1-w)\,\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\left(\mathbf{G} - \mathbf{A}_{22}\right)\mathbf{A}_{22}^{-1}\mathbf{A}_{21}\right]\sigma^2_a
```

In the same way, the whole matrix is $\mathbf{H}_w = (1-w)\,\mathbf{H} + w\,\mathbf{A}$: a weighted average of $\mathbf{H}$ built with the unblended $\mathbf{G}$ and the pedigree matrix $\mathbf{A}$. In $\mathbf{H}^{-1}$, $\mathbf{G}^{-1}$ is $\mathbf{G}_w^{-1}$.

Blending fixes the singularity of $\mathbf{G}$, but not the difference in scale between $\mathbf{G}$ and $\mathbf{A}_{22}$. Tuning $\mathbf{G}$ to $\mathbf{A}_{22}$ is a separate step and is not done here.

### 8.5 Variance of a

**Naive approach (incorrect).** Replace $\mathbf{A}_{22}$ by $\mathbf{G}$ and keep the other blocks of $\mathbf{A}$:

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

This is wrong. The non-genotyped animals are related to the genotyped ones through the pedigree, so when $\text{Var}[\mathbf{a}_2]$ changes from $\mathbf{A}_{22}\sigma^2_a$ to $\mathbf{G}\sigma^2_a$, then $\text{Var}[\mathbf{a}_1]$ and $\text{Cov}[\mathbf{a}_1, \mathbf{a}_2^\top]$ must change too. The correct blocks are derived below.

**Decomposition of $\mathbf{a}_1$.** From [8.1](#81-the-pedigree-propagates-genomic-information):

```math
\mathbb{E}[\mathbf{a}_1 \mid \mathbf{a}_2]
= \mathbf{A}_{12}\sigma^2_a\left(\mathbf{A}_{22}\sigma^2_a\right)^{-1}\mathbf{a}_2
= \mathbf{A}_{12}\sigma^2_a\mathbf{A}_{22}^{-1}\frac{1}{\sigma^2_a}\mathbf{a}_2
= \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{a}_2
= \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{Zg}
```

Then $\mathbf{a}_1$ is its conditional expectation plus a deviation $\boldsymbol{\varepsilon}$:

```math
\mathbf{a}_1 = \mathbb{E}[\mathbf{a}_1 \mid \mathbf{a}_2] + \boldsymbol{\varepsilon}
= \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{a}_2 + \boldsymbol{\varepsilon}
= \mathbf{D}\mathbf{a}_2 + \boldsymbol{\varepsilon}
\quad \Longrightarrow \quad
\boldsymbol{\varepsilon} = \mathbf{a}_1 - \mathbf{D}\mathbf{a}_2,
\qquad
\mathbf{D} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1} \in \mathbb{R}^{n_{a_1} \times n_{a_2}}
```

$\mathbf{D}\mathbf{a}_2$ is the part of $\mathbf{a}_1$ predicted from the genotyped animals through the pedigree. $\boldsymbol{\varepsilon}$ is the part that $\mathbf{a}_2$ cannot explain (e.g., the Mendelian sampling of the non-genotyped animals). Note that $\mathbf{D}^\top = (\mathbf{A}_{22}^{-1})^\top\mathbf{A}_{12}^\top = \mathbf{A}_{22}^{-1}\mathbf{A}_{21}$, because $\mathbf{A}_{22}^{-1}$ is symmetric and $\mathbf{A}_{12}^\top = \mathbf{A}_{21}$.

**Key assumption (Legarra et al., 2009).** The genomic information changes only the distribution of $\mathbf{a}_2$: $\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a$ instead of $\mathbf{A}_{22}\sigma^2_a$. The link between $\mathbf{a}_1$ and $\mathbf{a}_2$ ($\mathbf{D}$ and $\boldsymbol{\varepsilon}$) is still the one given by the pedigree. $\boldsymbol{\varepsilon}$ and its variance are therefore computed with the pedigree (co)variances, written $\text{Var}_{\mathbf{A}}$ and $\text{Cov}_{\mathbf{A}}$:

```math
\text{Var}_{\mathbf{A}}[\mathbf{a}_1] = \mathbf{A}_{11}\sigma^2_a, \quad
\text{Var}_{\mathbf{A}}[\mathbf{a}_2] = \mathbf{A}_{22}\sigma^2_a, \quad
\text{Cov}_{\mathbf{A}}\left[\mathbf{a}_1, \mathbf{a}_2^\top\right] = \mathbf{A}_{12}\sigma^2_a, \quad
\text{Cov}_{\mathbf{A}}\left[\mathbf{a}_2, \mathbf{a}_1^\top\right] = \mathbf{A}_{21}\sigma^2_a
```

#### Covariance between $\mathbf{D}\mathbf{a}_2$ and $\boldsymbol{\varepsilon}$

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

By the key assumption, $\boldsymbol{\varepsilon}$ stays independent of $\mathbf{a}_2$ when $\text{Var}[\mathbf{a}_2]$ changes to $\mathbf{G}\sigma^2_a$. Then

```math
\text{Cov}\left[\mathbf{D}\mathbf{a}_2, \boldsymbol{\varepsilon}^\top\right]
= \text{Cov}\left[\boldsymbol{\varepsilon}, \left(\mathbf{D}\mathbf{a}_2\right)^\top\right]^\top
= \mathbf{0} \in \mathbb{R}^{n_{a_1} \times n_{a_1}}
\qquad \Longrightarrow \qquad
\text{Var}[\mathbf{a}_1] = \text{Var}[\mathbf{D}\mathbf{a}_2 + \boldsymbol{\varepsilon}] = \text{Var}[\mathbf{D}\mathbf{a}_2] + \text{Var}[\boldsymbol{\varepsilon}]
```

#### Variance of $\mathbf{D}\mathbf{a}_2$

Here $\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a$ (genomic information):

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

#### Variance of $\boldsymbol{\varepsilon}$

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

#### Variance of $\mathbf{a}_1$

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

The second line uses $\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{22}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}$ to factor the last two terms. Interpretation:

- $\mathbf{A}_{11}$: pedigree relationships among non-genotyped animals.
- $\mathbf{G} - \mathbf{A}_{22}$: what the markers add to (or remove from) the pedigree relationships among genotyped animals.
- $\mathbf{A}_{12}\mathbf{A}_{22}^{-1}$: propagates this difference to the non-genotyped animals through their pedigree relationships with the genotyped animals.
- If $\mathbf{G} = \mathbf{A}_{22}$, then $\text{Var}[\mathbf{a}_1] = \mathbf{A}_{11}\sigma^2_a$ (pedigree only).

#### Variance of $\mathbf{a}_2$

```math
\text{Var}[\mathbf{a}_2] = \mathbf{G}\sigma^2_a
```

#### Covariance between $\mathbf{a}_1$ and $\mathbf{a}_2$

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

This is the transpose of $\text{Cov}_{\mathbf{A}}[\mathbf{a}_2, \boldsymbol{\varepsilon}^\top] = \mathbf{0}$, and by the key assumption it also holds with genomic information. Then

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

because $\mathbf{G}$ and $\mathbf{A}_{22}^{-1}$ are symmetric and $\mathbf{A}_{12}^\top = \mathbf{A}_{21}$.

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

$\mathbf{H}_{21} = \mathbf{H}_{12}^\top$, so $\mathbf{H}$ is symmetric. $\mathbf{H}$ is the relationship matrix that combines pedigree and genomic relationships (Legarra et al., 2009), with the animals in the order of $\mathbf{A}_{\text{ord}}$. In compact form, with $\mathbf{D} = \mathbf{A}_{12}\mathbf{A}_{22}^{-1}$:

```math
\mathbf{H} = \mathbf{A} + \begin{pmatrix}\mathbf{D} \\ \mathbf{I}\end{pmatrix}\left(\mathbf{G} - \mathbf{A}_{22}\right)\begin{pmatrix}\mathbf{D} \\ \mathbf{I}\end{pmatrix}^\top
```

Multiplying out gives the four blocks above.

### 8.6 Inverse of H

The MME need $\mathbf{H}^{-1}$, not $\mathbf{H}$. Aguilar et al. (2010) showed that

```math
\mathbf{H}^{-1} = \mathbf{A}^{-1} +
\begin{pmatrix}
\mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{G}^{-1} - \mathbf{A}_{22}^{-1}
\end{pmatrix}
```

Only the block of the genotyped animals changes. Three inverses are needed:

1. $\mathbf{A}^{-1}$: all animals, in the order of $\mathbf{A}_{\text{ord}}$.
2. $\mathbf{A}_{22}^{-1}$: the inverse of the block $\mathbf{A}_{22}$. This is **not** the 22 block of $\mathbf{A}^{-1}$. It is the plain pedigree block and is **not** blended.
3. $\mathbf{G}^{-1} = \mathbf{G}_w^{-1}$: the inverse of the blended matrix ([8.4](#84-blending-g)).

$\mathbf{A}$ and $\mathbf{A}_{22}$ are symmetric positive definite and are inverted from their Cholesky factors. Here $\mathbf{A}^{-1}$ is obtained by inverting the dense $\mathbf{A}$, which is fine for this pedigree. For large pedigrees, it is built directly from the pedigree (Henderson's rules), because it is sparse.

**Checks.**

1. $\mathbf{H}^{-1}$ is symmetric, with the animals in the order of $\mathbf{A}_{\text{ord}}$.
2. The rows of the non-genotyped animals are identical to those of $\mathbf{A}^{-1}$.
3. $\mathbf{H}$ is built from the compact form in [8.5](#85-variance-of-a) with $\mathbf{G} = \mathbf{G}_w$. The code checks that $\mathbf{H}_{22} = \mathbf{G}_w$ and that $\max\lvert\mathbf{H}^{-1}\mathbf{H} - \mathbf{I}\rvert < 10^{-8}$. $\mathbf{H}$ is kept, because REML uses it.

---

## 9. Mixed model equations

The model is $\mathbf{y} = \mathbf{Xb} + \mathbf{Wa} + \mathbf{e}$ with

```math
\mathbf{a} \sim N\left(\mathbf{0}, \mathbf{H}\sigma^2_a\right), \quad \mathbf{H}\sigma^2_a = \mathbf{Q} \in \mathbb{R}^{K \times K},
\qquad
\mathbf{e} \sim N\left(\mathbf{0}, \mathbf{I}\sigma^2_e\right), \quad \mathbf{I}\sigma^2_e = \mathbf{R} \in \mathbb{R}^{N \times N},
\qquad
\text{Cov}\left[\mathbf{a}, \mathbf{e}^\top\right] = \mathbf{0}
```

where $K$ is the number of animals in $\mathbf{a}$ (all animals in the pedigree) and $N$ the number of phenotypes. The MME (Henderson, 1984) are

```math
\begin{pmatrix}
\mathbf{X}^\top\mathbf{R}^{-1}\mathbf{X} & \mathbf{X}^\top\mathbf{R}^{-1}\mathbf{W} \\
\mathbf{W}^\top\mathbf{R}^{-1}\mathbf{X} & \mathbf{W}^\top\mathbf{R}^{-1}\mathbf{W} + \mathbf{H}^{-1}\frac{1}{\sigma^2_a}
\end{pmatrix}
\begin{pmatrix} \widehat{\mathbf{b}} \\ \widehat{\mathbf{a}} \end{pmatrix}
=
\begin{pmatrix} \mathbf{X}^\top\mathbf{R}^{-1}\mathbf{y} \\ \mathbf{W}^\top\mathbf{R}^{-1}\mathbf{y} \end{pmatrix}
```

with $\mathbf{H}^{-1}\frac{1}{\sigma^2_a} = \left(\mathbf{H}\sigma^2_a\right)^{-1} = \mathbf{Q}^{-1}$ and $\mathbf{R}^{-1} = \left(\mathbf{I}\sigma^2_e\right)^{-1} = \mathbf{I}\frac{1}{\sigma^2_e}$. Then

```math
\begin{pmatrix}
\mathbf{X}^\top\frac{1}{\sigma^2_e}\mathbf{X} & \mathbf{X}^\top\frac{1}{\sigma^2_e}\mathbf{W} \\
\mathbf{W}^\top\frac{1}{\sigma^2_e}\mathbf{X} & \mathbf{W}^\top\frac{1}{\sigma^2_e}\mathbf{W} + \mathbf{H}^{-1}\frac{1}{\sigma^2_a}
\end{pmatrix}
\begin{pmatrix} \widehat{\mathbf{b}} \\ \widehat{\mathbf{a}} \end{pmatrix}
=
\begin{pmatrix} \mathbf{X}^\top\frac{1}{\sigma^2_e}\mathbf{y} \\ \mathbf{W}^\top\frac{1}{\sigma^2_e}\mathbf{y} \end{pmatrix}
```

Multiplying both sides by $\sigma^2_e$:

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

$\mathbf{W} \in \mathbb{R}^{N \times K}$ has one 1 per row, in the column of the animal with the record. The columns of animals without a record ($A_{100}$, $A_{101}$) are zero.

---

## 10. Restricted maximum likelihood

The MME need $\alpha = \sigma^2_e/\sigma^2_a$, so the variance components $\boldsymbol{\theta} = (\sigma^2_a, \sigma^2_e)^\top$ must be estimated first. This is done by restricted (or residual) maximum likelihood, REML (Patterson and Thompson, 1971).

### 10.1 Why REML and not ML

Under the model, $\mathbf{y} \sim N(\mathbf{Xb}, \mathbf{V})$, with

```math
\mathbf{V} = \text{Var}[\mathbf{y}] = \text{Var}[\mathbf{Wa} + \mathbf{e}] = \mathbf{W}\mathbf{H}\mathbf{W}^\top\sigma^2_a + \mathbf{I}\sigma^2_e
```

because $\mathbf{a}$ and $\mathbf{e}$ are independent. $\mathbf{W}\mathbf{H}\mathbf{W}^\top$ is the block of $\mathbf{H}$ for the phenotyped animals, in the order of $\mathbf{y}$.

Maximum likelihood (ML) maximizes the likelihood of $\mathbf{y}$ as if $\mathbf{b}$ were known. It ignores that $\mathbf{b}$ is estimated from the same data, so the variances are biased downwards. For example, in $\mathbf{y} = \mathbf{Xb} + \mathbf{e}$, ML gives $\hat{\sigma}^2_e = \text{SSE}/N$ instead of $\text{SSE}/(N-p)$.

REML uses only the part of the data that is free of $\mathbf{b}$. Let $\mathbf{S}$ be any $N \times (N-p)$ matrix of full column rank with $\mathbf{S}^\top\mathbf{X} = \mathbf{0}$ (for example, $N - p$ linearly independent columns of $\mathbf{I} - \mathbf{X}(\mathbf{X}^\top\mathbf{X})^{-1}\mathbf{X}^\top$). Each element of $\mathbf{S}^\top\mathbf{y}$ is a linear combination of the records in which the fixed effects cancel out (an *error contrast*):

```math
\mathbf{S}^\top\mathbf{y} = \mathbf{S}^\top\mathbf{Xb} + \mathbf{S}^\top(\mathbf{Wa} + \mathbf{e}) = \mathbf{S}^\top(\mathbf{Wa} + \mathbf{e}) \sim N\left(\mathbf{0}, \mathbf{S}^\top\mathbf{V}\mathbf{S}\right)
```

REML is ML applied to $\mathbf{S}^\top\mathbf{y}$ instead of $\mathbf{y}$. Its distribution does not depend on $\mathbf{b}$, and the $p$ degrees of freedom used by the fixed effects are accounted for ($N - p$ contrasts instead of $N$ records). The result does not depend on the choice of $\mathbf{S}$, and $\mathbf{S}$ is never computed: the REML log-likelihood is written with $\mathbf{V}$ and $\mathbf{X}$ only.

### 10.2 REML log-likelihood

```math
L(\boldsymbol{\theta}) = -\frac{1}{2}\left[\log\left|\mathbf{V}\right| + \log\left|\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\right| + \mathbf{y}^\top\mathbf{P}\mathbf{y}\right] + \text{constant}
```

```math
\mathbf{P} = \mathbf{V}^{-1} - \mathbf{V}^{-1}\mathbf{X}\left(\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\right)^{-1}\mathbf{X}^\top\mathbf{V}^{-1}
```

$\mathbf{P}\mathbf{X} = \mathbf{0}$, so $\mathbf{P}$ removes the fixed effects. Also, $\mathbf{P}\mathbf{y} = \mathbf{V}^{-1}(\mathbf{y} - \mathbf{X}\hat{\mathbf{b}})$, where $\hat{\mathbf{b}}$ is the generalized least squares solution, $\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\hat{\mathbf{b}} = \mathbf{X}^\top\mathbf{V}^{-1}\mathbf{y}$. The three terms of $L$ are:

1. $\mathbf{y}^\top\mathbf{P}\mathbf{y} = (\mathbf{y} - \mathbf{X}\hat{\mathbf{b}})^\top\mathbf{V}^{-1}(\mathbf{y} - \mathbf{X}\hat{\mathbf{b}})$: a weighted sum of squares of residuals.
2. $\log\lvert\mathbf{V}\rvert$: depends only on the variance structure.
3. $\log\lvert\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\rvert$: a penalty for estimating $\mathbf{b}$. ML does not have this term.

### 10.3 Maximizing L

$L$ has no closed-form maximum, so $\boldsymbol{\theta}$ is found iteratively from the derivatives of $L$. $\mathbf{V}$ is linear in $\boldsymbol{\theta}$:

```math
\mathbf{V}_a = \frac{\partial\mathbf{V}}{\partial\sigma^2_a} = \mathbf{W}\mathbf{H}\mathbf{W}^\top,
\qquad
\mathbf{V}_e = \frac{\partial\mathbf{V}}{\partial\sigma^2_e} = \mathbf{I}
```

**Score** (first derivatives), for $i \in \lbrace a, e \rbrace$:

```math
\frac{\partial L}{\partial\theta_i} = -\frac{1}{2}\left[\text{tr}\left(\mathbf{P}\mathbf{V}_i\right) - \mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{y}\right]
```

At the maximum the score is zero: the quadratic form of the data, $\mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{y}$, equals its expected value, $\text{tr}(\mathbf{P}\mathbf{V}_i)$.

**Information matrices** ($2 \times 2$), element $(i, j)$:

| Information | Element $(i, j)$ | Method |
|:---|:---|:---|
| Observed, $-\partial^2 L/\partial\theta_i\partial\theta_j$ | $\mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j\mathbf{P}\mathbf{y} - \frac{1}{2}\text{tr}(\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j)$ | Newton-Raphson |
| Expected (Fisher) | $\frac{1}{2}\text{tr}(\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j)$ | Fisher scoring |
| Average (AI) | $\frac{1}{2}\mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{V}_j\mathbf{P}\mathbf{y}$ | AI-REML (Gilmour et al., 1995) |

The AI matrix is the average of the observed and the expected information. It has no trace terms, so it is cheap to compute. Each iteration updates

```math
\boldsymbol{\theta}^{(t+1)} = \boldsymbol{\theta}^{(t)} + \left(\mathbf{AI}^{(t)}\right)^{-1}\frac{\partial L}{\partial\boldsymbol{\theta}}\bigg|_{\boldsymbol{\theta}^{(t)}}
```

At convergence, $\text{Var}[\hat{\boldsymbol{\theta}}] \approx \mathbf{AI}^{-1}$, which gives the standard errors.

### 10.4 Computing L and its derivatives from the MME

When $N$ is large, $\mathbf{V}$ ($N \times N$) is too large to invert, and everything is obtained from the MME. Let $\mathbf{T} = \begin{pmatrix}\mathbf{X} & \mathbf{W}\end{pmatrix}$ and

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

With the MME solutions $\hat{\mathbf{b}}$, $\hat{\mathbf{a}}$ and the residuals $\hat{\mathbf{e}} = \mathbf{y} - \mathbf{X}\hat{\mathbf{b}} - \mathbf{W}\hat{\mathbf{a}}$:

```math
\mathbf{P}\mathbf{y} = \frac{\hat{\mathbf{e}}}{\sigma^2_e},
\qquad
\mathbf{y}^\top\mathbf{P}\mathbf{y} = \frac{\mathbf{y}^\top\hat{\mathbf{e}}}{\sigma^2_e}
```

**Log-likelihood** (e.g., Meyer, 1989), with $p$ the number of columns of $\mathbf{X}$ (full rank). $\log\lvert\mathbf{H}\rvert$ does not depend on $\boldsymbol{\theta}$; it is written explicitly so that this $L$ equals the one in 10.2, not only up to a constant:

```math
-2L = (N - p - K)\log\sigma^2_e + K\log\sigma^2_a + \log\left|\mathbf{H}\right| + \log\left|\mathbf{C}\right| + \frac{\mathbf{y}^\top\hat{\mathbf{e}}}{\sigma^2_e} + \text{constant}
```

**Score.** $\mathbf{C}^{22}$ is the block of $\mathbf{C}^{-1}$ for $\mathbf{a}$, with $\text{Var}[\mathbf{a} - \hat{\mathbf{a}}] = \mathbf{C}^{22}\sigma^2_e$ (prediction error variance, PEV):

```math
\frac{\partial L}{\partial\sigma^2_a} =
\frac{1}{2}\left[\frac{\hat{\mathbf{a}}^\top\mathbf{H}^{-1}\hat{\mathbf{a}}}{\sigma^4_a} - \frac{K}{\sigma^2_a} + \frac{\text{tr}\left(\mathbf{H}^{-1}\mathbf{C}^{22}\right)\sigma^2_e}{\sigma^4_a}\right]
```

```math
\frac{\partial L}{\partial\sigma^2_e} =
\frac{1}{2}\left[\frac{\hat{\mathbf{e}}^\top\hat{\mathbf{e}}}{\sigma^4_e} - \frac{N - p - K}{\sigma^2_e} - \frac{\text{tr}\left(\mathbf{H}^{-1}\mathbf{C}^{22}\right)}{\sigma^2_a}\right]
```

**AI matrix.** The vectors $\mathbf{V}_i\mathbf{P}\mathbf{y}$ are

```math
\mathbf{V}_a\mathbf{P}\mathbf{y} = \frac{\mathbf{W}\hat{\mathbf{a}}}{\sigma^2_a},
\qquad
\mathbf{V}_e\mathbf{P}\mathbf{y} = \frac{\hat{\mathbf{e}}}{\sigma^2_e}
```

With $\mathbf{F} = \begin{pmatrix}\mathbf{W}\hat{\mathbf{a}}/\sigma^2_a & \hat{\mathbf{e}}/\sigma^2_e\end{pmatrix} \in \mathbb{R}^{N \times 2}$, $\mathbf{AI} = \frac{1}{2}\mathbf{F}^\top\mathbf{P}\mathbf{F}$, and by the key identity

```math
\mathbf{AI} = \frac{1}{2\sigma^2_e}\left(\mathbf{F}^\top\mathbf{F} - \mathbf{F}^\top\mathbf{T}\mathbf{C}^{-1}\mathbf{T}^\top\mathbf{F}\right)
```

That is, the MME are solved with the two columns of $\mathbf{F}$ in place of $\mathbf{y}$ (right-hand side $\mathbf{T}^\top\mathbf{F}$).

### 10.5 EM-REML

Setting both scores to zero and rearranging gives (Mrode, 2014, Eqns 17.5 and 17.6)

```math
\sigma^2_e = \frac{\mathbf{y}^\top\hat{\mathbf{e}}}{N - p},
\qquad
\sigma^2_a = \frac{\hat{\mathbf{a}}^\top\mathbf{H}^{-1}\hat{\mathbf{a}} + \text{tr}\left(\mathbf{H}^{-1}\mathbf{C}^{22}\right)\sigma^2_e}{K}
```

Used as an update (right-hand sides computed with the current $\boldsymbol{\theta}$), this is the EM algorithm (Dempster et al., 1977). $\sigma^2_e$ comes from the residuals, and $\sigma^2_a$ from the predicted breeding values plus their PEV. EM estimates always stay positive, but convergence is slow. AI-REML converges in a few iterations, but an update can give a negative variance. A common strategy is to start with a few EM iterations, or to use an EM step when the AI step leaves the parameter space.

### 10.6 AI-REML algorithm

1. **Start values**, e.g., from a guess of heritability $h^2_0$ and the phenotypic variance $\sigma^2_y$: $\sigma^2_a = h^2_0\sigma^2_y$ and $\sigma^2_e = (1 - h^2_0)\sigma^2_y$.
2. **Compute what $L$, the score and AI need.** Either $\mathbf{V}$ and $\mathbf{P}$ (small $N$), or the MME solutions $\hat{\mathbf{b}}$, $\hat{\mathbf{a}}$, $\hat{\mathbf{e}}$ and $\mathbf{C}^{22}$ with $\alpha = \sigma^2_e/\sigma^2_a$ (large $N$).
3. **Compute** $L$, the score and the AI matrix.
4. **Update** $\boldsymbol{\theta}$. A variance cannot be negative. If the update makes it negative, it is set to a small positive value (the boundary), and it stays there while its score is negative. An EM step is an alternative.
5. **Repeat** steps 2–4 until convergence (changes in $L$ and $\boldsymbol{\theta}$ below a small tolerance).
6. **Report** $\hat{\sigma}^2_a$, $\hat{\sigma}^2_e$, their standard errors from $\mathbf{AI}^{-1}$, and $\hat{h}^2 = \hat{\sigma}^2_a/(\hat{\sigma}^2_a + \hat{\sigma}^2_e)$ with its standard error by the delta method.

**Note on Mrode (2014, Chapter 17).** Two formulas there seem to have misprints. $\mathbf{P}$ is printed without the first $\mathbf{V}^{-1}$ in the second term, and in Eqns 17.2 and 17.3 some signs and powers of the variances differ from the ones above. The expressions here are consistent with Eqns 17.5 and 17.6 and were checked numerically (Section 11.4).

---

## 11. AI-REML implementation

### 11.1 Which route

Here there are far fewer records $N$ than MME equations $K + p$ (the pedigree has many more animals than phenotypes). Factorizing $\mathbf{V}$ ($N \times N$) is then much cheaper than factorizing $\mathbf{C}$ ($(K+p) \times (K+p)$). $L$, the score and the AI matrix are therefore computed directly with the formulas of [10.3](#103-maximizing-l), one line of code per formula. The MME route of [10.4](#104-computing-l-and-its-derivatives-from-the-mme) gives the same numbers and is used as a check. With many records, the MME route with sparse matrices is the one to use.

REML needs only $\mathbf{y}$, $\mathbf{X}_{\text{ord}}$ and $\mathbf{W}\mathbf{H}\mathbf{W}^\top$: $\mathbf{V}$ depends on $\mathbf{H}$ only through the phenotyped animals. In the code, `Hyy = H[y_ids, y_ids]` is $\mathbf{W}\mathbf{H}\mathbf{W}^\top$, in the order of $\mathbf{y}$.

### 11.2 Steps 2 and 3: `reml_terms()`

For a given $\boldsymbol{\theta}$:

| Quantity | Formula | Computation |
|:---|:---|:---|
| $\mathbf{V}$ | $\mathbf{W}\mathbf{H}\mathbf{W}^\top\sigma^2_a + \mathbf{I}\sigma^2_e$ | `sa * Hyy + se * diag(N)` |
| $\mathbf{V}^{-1}$, $\log\lvert\mathbf{V}\rvert$ | Cholesky $\mathbf{V} = \mathbf{U}^\top\mathbf{U}$; $\log\lvert\mathbf{V}\rvert = 2\sum_i\log u_{ii}$ | `chol()`, `chol2inv()` |
| $\mathbf{P}$ | $\mathbf{V}^{-1} - \mathbf{V}^{-1}\mathbf{X}(\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}^\top\mathbf{V}^{-1}$ | `solve(XViX, t(ViX))` |
| $L$ | $-\frac{1}{2}\left[\log\lvert\mathbf{V}\rvert + \log\lvert\mathbf{X}^\top\mathbf{V}^{-1}\mathbf{X}\rvert + \mathbf{y}^\top\mathbf{P}\mathbf{y}\right]$ | `determinant()` for the second term |
| $\mathbf{V}_a\mathbf{P}\mathbf{y}$, $\mathbf{V}_e\mathbf{P}\mathbf{y}$ | $\mathbf{W}\mathbf{H}\mathbf{W}^\top\mathbf{P}\mathbf{y}$, $\mathbf{P}\mathbf{y}$ | |
| Score | $-\frac{1}{2}\left[\text{tr}(\mathbf{P}\mathbf{V}_i) - \mathbf{y}^\top\mathbf{P}\mathbf{V}_i\mathbf{P}\mathbf{y}\right]$ | $\text{tr}(\mathbf{P}\mathbf{V}_a) = \sum_{ij}p_{ij}(\mathbf{W}\mathbf{H}\mathbf{W}^\top)_{ij}$ (both symmetric), avoiding the matrix product |
| AI | $\frac{1}{2}\mathbf{F}^\top\mathbf{P}\mathbf{F}$, $\mathbf{F} = (\mathbf{V}_a\mathbf{P}\mathbf{y}, \mathbf{V}_e\mathbf{P}\mathbf{y})$ | `crossprod(Fm, P %*% Fm)` |

Here $p_{ij}$ is element $(i, j)$ of $\mathbf{P}$ (not an allele frequency).

### 11.3 Steps 1, 4, 5 and 6: `ai_reml()`

1. **Start values.** $\sigma^2_y$ is the residual variance of $\mathbf{y}$ after fitting only the fixed effects by ordinary least squares, $\sigma^2_y = \text{SSE}/(N - p)$. The start heritability is $h^2_0 = 0.3$.
2. **Update.** $\boldsymbol{\theta} \leftarrow \boldsymbol{\theta} + \mathbf{AI}^{-1}\,\partial L/\partial\boldsymbol{\theta}$. If the update makes a variance $\le 0$, it is set to $\epsilon = 10^{-8}\sigma^2_y$ (the boundary). While a variance is at the boundary and its score is negative ($L$ would only increase if the variance could be negative), it stays fixed. Only the other variance is updated, with its own element of $\mathbf{AI}$ (an active-set update).
3. **Convergence.** Stop when the largest change in $\boldsymbol{\theta}$, relative to $\sigma^2_y$, and the change in $L$ are both below $10^{-8}$.
4. **Estimates and standard errors.** $\text{Var}[\hat{\boldsymbol{\theta}}] \approx \mathbf{AI}^{-1}$. For $\hat{h}^2$, the delta method gives

```math
\text{Var}\left[\hat{h}^2\right] \approx \left(\nabla h^2\right)^\top\mathbf{AI}^{-1}\left(\nabla h^2\right),
\qquad
\nabla h^2 = \frac{\partial h^2}{\partial\boldsymbol{\theta}} = \frac{1}{\left(\sigma^2_a + \sigma^2_e\right)^2}\begin{pmatrix}\sigma^2_e \\ -\sigma^2_a\end{pmatrix}
```

If a variance is at the boundary, the large-sample theory does not apply. Its standard error, and that of $\hat{h}^2$, are not reported (`NA`), and the standard error of the other variance uses only its own element of $\mathbf{AI}$.

**Profile log-likelihood.** For each $h^2$ on a grid, $L$ is maximized over the total variance $\sigma^2_a + \sigma^2_e$:

```math
L_{\text{prof}}(h^2) = \max_{\sigma^2_a + \sigma^2_e}\; L\left(h^2\left(\sigma^2_a + \sigma^2_e\right),\; \left(1 - h^2\right)\left(\sigma^2_a + \sigma^2_e\right)\right)
```

Plotting $L_{\text{prof}}$ against $h^2$ shows the shape of the likelihood around $\hat{h}^2$, including when the estimate is at the boundary.

### 11.4 Checks

1. **The MME route gives the same $L$ and score.** At the start values, $L$ and the score are recomputed with the formulas of [10.4](#104-computing-l-and-its-derivatives-from-the-mme). This includes $\log\lvert\mathbf{H}\rvert$, so both values of $L$ are equal, not only up to a constant. The code builds $\mathbf{W}$, $\mathbf{T}$, $\mathbf{C}$, $\mathbf{C}^{-1}$, $\hat{\mathbf{a}}$, $\hat{\mathbf{e}}$, $\text{tr}(\mathbf{H}^{-1}\mathbf{C}^{22})$ and $\hat{\mathbf{a}}^\top\mathbf{H}^{-1}\hat{\mathbf{a}}$, and checks agreement within $10^{-6}$.
2. **A general-purpose optimizer finds the same maximum.** `optim()` (L-BFGS-B, without derivatives, with $\sigma^2_a, \sigma^2_e \ge 0$) maximizes $L$ directly. AI-REML must reach at least the same $L$.
3. **The code recovers known variances.** 100 data sets are simulated with the same $\mathbf{X}$ and $\mathbf{H}$, with $\sigma^2_a = 6$ and $\sigma^2_e = 9$ ($h^2 = 0.4$). If $\mathbf{z} \sim N(\mathbf{0}, \mathbf{I})$ and $\mathbf{U}_H$ is the Cholesky factor of $\mathbf{W}\mathbf{H}\mathbf{W}^\top$ ($\mathbf{W}\mathbf{H}\mathbf{W}^\top = \mathbf{U}_H^\top\mathbf{U}_H$), then $\mathbf{U}_H^\top\mathbf{z} \sim N(\mathbf{0}, \mathbf{W}\mathbf{H}\mathbf{W}^\top)$. So

```math
\mathbf{y}_{\text{sim}} = \mathbf{X}_{\text{ord}}\mathbf{1}_p + \sigma_a\,\mathbf{U}_H^\top\mathbf{z}_1 + \sigma_e\,\mathbf{z}_2,
\qquad \mathbf{z}_1, \mathbf{z}_2 \sim N(\mathbf{0}, \mathbf{I}_N)
```

The fixed effects do not change the REML estimates, so all are set to 1. The mean of the 100 estimates should be close to the true values (REML is not exactly unbiased in small samples).

---

## 12. Conclusion

There are few phenotypes. With the current data, the REML estimate of $\sigma^2_a$ for `delta` goes to the boundary, so $\hat{h}^2 \simeq 0$. At that point $\mathbf{V} = \mathbf{I}\sigma^2_e$, and $\hat{\sigma}^2_e$ is the residual variance of the model with only the fixed effects. The profile log-likelihood decreases as $h^2$ increases.

With this amount of data, $\sigma^2_a$ (and therefore $h^2$) cannot be estimated reliably. Subsequent predictions of additive genetic values, and their propagation to SNP effects, would be unreliable and misleading.

At the same time, the checks in [11.4](#114-checks) show that the pipeline itself is correct. The two routes agree, a general-purpose optimizer reaches the same maximum, and the method recovers known variance components in simulated data with the same structure.

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
