# baGWAS
### An Inter-Breed GWAS Framework Integrating Individual-Level Genotype and Breed-Average Phenotypes

**baGWAS** is an R package for inter-breed genome-wide association analysis using individual-level genotype data and breed-average phenotypes through two complementary workflows:

- **`genoGWAS`** — breed-level genotype data-based GWAS using GEMMA linear mixed models.
- **`freqGWAS`** — breed-level allele frequency-based GWAS using group-level allele-frequency matrices, an IBS-derived breed relationship matrix, and EMMREML.

The package provides a single main entry point, `baGWAS()`, and returns a standardized GWAS result table together with reproducible output files, diagnostic plots, and a run log.

> **Status:** Development release 1.0.0. The current release is intended for Linux/macOS (Unix-like systems) and requires external command-line software for the GWAS workflows.

## Main function

```r
baGWAS.Res <- baGWAS(
  strategy = "",
  organism = "",
  pheno = pheno,
  covariate = NULL,
  trait = "",
  geno = ""
)
```

Choose one of the two strategies:

```r
strategy = "genoGWAS"
```

or

```r
strategy = "freqGWAS"
```

## Installation

baGWAS is designed to run on Ubuntu Linux and requires R packages and external command-line software.

## 1. Install external software on Ubuntu

Open the Ubuntu Terminal and install PLINK 1.9 and GEMMA:
```r
sudo apt update
sudo apt install plink1.9 gemma
```

### Software requirements by strategy:
- genoGWAS — requires PLINK for genotype-data preparation and GEMMA for linear mixed-model association analysis.
- freqGWAS — requires PLINK for breed-level allele-frequency and IBS calculations.

> For the current implementation, PLINK 1.9 is the recommended starting point.
- These external programs must be installed separately and accessible through the system `PATH` or the environment-variable configuration supported by baGWAS.

## 2. Install R package dependencies

The following R packages are required by baGWAS:

- `data.table`
- `EMMREML`
- `dplyr`
- `ggplot2`
- `ggmanh`

The first four packages are available from CRAN. The `ggmanh` package is distributed through Bioconductor.

### Install R dependencies

Install the CRAN packages:

```r
install.packages(
  c("data.table", "EMMREML", "dplyr", "ggplot2"),
  repos = "https://cloud.r-project.org"
)
```

Install `ggmanh` from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

BiocManager::install("ggmanh")
```

### 3. Install baGWAS from GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("fallahi-m/baGWAS")
```

Load the package:

```r
library(baGWAS)
```

## Quick start

### 1. genoGWAS (`breed-level genotype data-based GWAS`)

```r
library(baGWAS)

baGWAS.Res <- baGWAS(
  strategy = "genoGWAS",
  organism = "dog",
  pheno = pheno,
  covariate = covariate,
  trait = "weight",
  geno = "/data/geno"
)
```

### 2. freqGWAS (`breed-level allele frequency-based GWAS`)

For `freqGWAS`, `group_col` identifies the group/breed/population used to calculate allele frequencies and group-level relationships.

```r
baGWAS.Res <- baGWAS(
  strategy = "freqGWAS",
  organism = "dog",
  pheno = pheno,
  covariate = covariate,
  trait = "weight",
  geno = "/data/geno",
  group_col = "FID"
)
```

The `covariate` argument may be `NULL`. In that case, the model uses an intercept only.


## Choosing a Strategy

baGWAS provides two complementary strategies for genome-wide association analysis using breed-average phenotypes. Both address the same broad question: which genetic loci show variation across breeds that is associated with differences in breed-mean phenotypes, after accounting for population structure? They differ in genotype representation, unit of analysis, and the way genetic relatedness is modelled. Their results may therefore be correlated, but they should not be expected to be identical.

### Conceptual Basis

In this framework, phenotypes are available as breed means rather than as individually measured, matched phenotypes. In `genoGWAS`, the breed-mean phenotype is assigned to each genotyped individual belonging to that breed. Consequently, the phenotype has no within-breed variation by construction. Within-breed genotype differences therefore do not provide direct within-breed genotype–phenotype information. The association signal is driven primarily by differences in genotype or allele frequency between breeds with different mean phenotypes.

In `freqGWAS`, genotypes are summarized as breed-level allele frequencies and analysed alongside the corresponding breed-mean phenotypes. This makes the breed the explicit unit of analysis and aligns the data representation directly with the level at which the phenotype is available.

### At-a-Glance Comparison

| Feature | `genoGWAS` | `freqGWAS` |
|---|---|---|
| Genotype representation | Individual-level genotypes | Breed-level allele frequencies |
| Unit represented in the association data | Individuals, with breed-mean phenotypes repeated within each breed | Breeds/groups, with one phenotype observation per breed/group |
| Association software/model | GEMMA linear mixed model; likelihood-ratio test in the current workflow | EMMREML mixed-model workflow |
| Relatedness matrix | Individual-level genomic relationship matrix | Breed/group-level IBS-derived relationship matrix |
| Covariates | Individual-level covariates, when supported by the workflow | Covariates summarized at the breed/group level |
| Main source of association signal | Differences in genotype distribution between breeds with different mean phenotypes | Covariation between breed-level allele frequencies and breed-mean phenotypes |
| Computational considerations | Generally more demanding when the individual-level relationship matrix is large | Generally less demanding when the number of breeds is much smaller than the number of individuals |
| Main caution | Repeated breed-mean phenotypes can create a mismatch between the nominal individual count and the number of independent phenotypic units | Breed-level allele frequencies can be imprecise when estimated from small samples; this uncertainty may not be fully represented in the model |

*Implementation details, including available covariates, quality-control filters, and test statistics, depend on the package version and analysis settings.*

### Key Differences That Can Affect Results

#### 1. Breed Sample Size and Weighting

In `genoGWAS`, breeds represented by more genotyped individuals can contribute more individual records to model fitting. In `freqGWAS`, each breed/group is represented by a row in the breed-level analysis. However, one row per breed does not, by itself, guarantee equal statistical weight; the effective contribution of each breed depends on the model and its variance assumptions.

#### 2. Repeated Phenotypes and Effective Replication

In `genoGWAS`, individuals from the same breed share an identical assigned phenotype. The number of individuals is therefore not the number of independent phenotypic observations. The number of breeds and their structure are central to the available between-breed information.

A mixed model and relationship matrix may account for some dependence, but they do not automatically guarantee valid statistical calibration when breed-mean phenotypes are repeated across individuals. Test-statistic inflation should be assessed rather than assumed absent.

#### 3. Uncertainty in Allele-Frequency Estimates

In `freqGWAS`, the precision of a breed's allele-frequency estimate depends partly on the number of successfully genotyped individuals in that breed and on genotype quality. Small breed samples can produce noisy estimates. Unless the analysis explicitly models this estimation uncertainty, it may not be fully reflected in standard errors or test statistics.

#### 4. Population Structure

`genoGWAS` models relatedness at the individual level, whereas `freqGWAS` uses a relationship matrix defined among breeds/groups. These matrices describe population structure at different levels and are not interchangeable. Their different parameterizations can change how population structure and potential confounding are controlled.

#### 5. Quality Control

Genotype quality control and variant filtering should be appropriate to the data representation and performed before interpreting associations. Individual-level filtering and breed-level allele-frequency filtering are not necessarily equivalent. Review the actual filters used by the selected workflow and document the thresholds in reports or publications.

#### 6. Computational Cost

When there are many individuals but relatively few breeds, a breed-level relationship matrix can be substantially smaller than an individual-level matrix. In that setting, `freqGWAS` may require fewer computational resources. Actual runtime depends on the number of individuals, breeds, variants, and model settings.

### Decision Guide

| Study Situation | Suggested Approach |
|---|---|
| Only breed-level allele frequencies are available | `freqGWAS` |
| Individual-level genotypes are available and the goal is association with breed-average phenotypes | Consider both strategies, with careful validation of `genoGWAS` calibration |
| Concern about inflated test statistics from repeated breed-mean phenotypes | Start with `freqGWAS` as the more direct breed-level representation; assess calibration for either method |
| Many individuals but relatively few breeds, with limited computational resources | `freqGWAS` may be more computationally practical |
| Individual-level covariates or individual-level genotype QC are essential | Consider `genoGWAS`, provided its statistical assumptions are appropriate and checked |
| Breed sample sizes are highly unbalanced | Choose based on the intended contribution of breeds and the model's weighting; investigate sensitivity to sample-size imbalance |
| Some breeds have very small genotype sample sizes | Interpret both approaches cautiously; breed-level allele-frequency estimates may be noisy |
| Independent assessment of candidate associations is desired | Run both where feasible and compare results, while also evaluating statistical calibration |

### Recommendation

Because the phenotype is defined at the breed level, **`freqGWAS` is the more direct representation of the stated analysis unit**. `genoGWAS` can provide a complementary analysis by retaining individual-level genotype records, but its use with repeated breed-mean phenotypes requires particular attention to model assumptions and test-statistic calibration.

When individual-level genotypes are available, running both strategies can be informative. We recommend evaluating:

- **Genomic inflation and Q–Q plots:** Inspect the distribution of test statistics and possible inflation. A genomic inflation factor (\(\lambda_{GC}\)) is a diagnostic, not proof that a model is valid or invalid.
- **Agreement across strategies:** Compare the correlation of \(-\log_{10}(P)\) values and inspect whether the same genomic regions show similar patterns.
- **Overlap of candidate loci:** Report the overlap under a clearly stated significance threshold and multiple-testing procedure.
- **Robustness:** Investigate whether results change with breed sample-size imbalance, the number of breeds, covariate choices, and population-structure assumptions.
- **Calibration and error rates:** Where possible, use simulations or other appropriate validation to assess type I error, false-positive rates, and power under realistic population structure.

Signals supported by both strategies may be useful candidates for follow-up, but agreement alone does not establish that an association is statistically valid. Independent replication and biological interpretation remain important.

> **Important note:** `genoGWAS` and `freqGWAS` are not statistically equivalent and should not be expected to return identical *P*-values. Differences can arise from the unit of analysis, breed sample-size balance, allele-frequency estimation, relationship matrices, covariates, and model assumptions. The number of breeds is a key determinant of the independent between-breed phenotypic information available to either strategy.

## Input data

### 1. Genotype data

`geno` is a PLINK binary dataset prefix. The following files must exist:

```text
genotypes.bed
genotypes.bim
genotypes.fam
```

Example:

```r
geno = "/home/user/data/geno"
```

A path ending in `.bed`, `.bim`, or `.fam` is normalized to the dataset prefix.

### 2. Phenotype data

`pheno` must be a data frame containing:

```text
FID   IID   Sex   Breed   <Trait>

FID: group/family identifier assigned to each breed. Individuals belonging to the same breed must have the same FID.
IID: individual identifier.
Sex: sex of the individual.
Breed: breed or population name.
<Trait>: phenotype trait to be analyzed. The trait column must be numeric.

group_col — specifies the column used to define the breed/group structure. By default, group_col = "FID".
```

Example:

```r
pheno <- data.frame(
  FID = c("1", "1", "2"),
  IID = c("1", "2", "3"),
  Sex = c(1, 2, 1),
  Breed = c("Breed_A", "Breed_A", "Breed_B"),
  weight = c(42.1, 42.1, 47.3)
)
```


### 3. Covariates
`Covariates` must be a data frame containing:

```text
FID   IID   Sex   Breed   <Covariate>

FID: group/family identifier assigned to each breed. Individuals belonging to the same breed must have the same FID.
IID: individual identifier.
Sex: sex of the individual.
Breed: breed or population name.
<Covariate>: numeric covariate variable included as a fixed effect in the GWAS model. Multiple covariates can be provided as additional columns.
```

## Supported organisms

The current PLINK species mapping implemented by baGWAS is:

| `organism` | PLINK flag |
|---|---|
| `dog` | `--dog` |
| `sheep` | `--sheep` |
| `cattle` | `--cow` |
| `cow` | `--cow` |
| `horse` | `--horse` |
| `mouse` | `--mouse` |
| `rice` | `--rice` |

These entries describe the organism keywords implemented by the package. Users should confirm that the selected PLINK option is appropriate for their dataset and analysis.

## Output

For a run with `output_dir = "<output_dir>"`, the package creates:

```text
<output_dir>/
└── Result/
    ├── genoGWAS/
    │   ├── baGWAS.Res
    │   ├── manhattan.png
    │   ├── qqplot.png
    │   └── baGWAS.log
    └── freqGWAS/
        ├── baGWAS.Res
        ├── manhattan.png
        ├── qqplot.png
        ├── baGWAS.log
        └── SNPs_with_errors.txt   # only when SNP models fail
```

`baGWAS.Res` contains the standardized columns:

```text
CHR   SNP   BP   P
```



## Diagnostic plots

Each successful strategy writes:

- `manhattan.png` — Manhattan plot of the GWAS p-values.
- `qqplot.png` — Q-Q plot with confidence band and the genomic inflation statistic (`lambda`).

## log file

baGWAS records the strategy, organism, trait, sample and SNP counts, run timing, and diagnostic information in `baGWAS.log`.

Intermediate files are created in temporary working directories and are removed after the analysis finishes.

## Statistical scope

`genoGWAS` and `freqGWAS` are distinct analytical workflows and should not be interpreted as interchangeable methods. Before interpreting results, users should evaluate genotype quality control, phenotype definition, sample structure, covariates, group definitions, missing data, model assumptions, and multiple-testing procedures.

The package performs the computational workflow; it does not replace study-specific quality control, model diagnostics, or biological interpretation.

## Testing

The repository contains `testthat` tests for package-level validation and core R functionality. Full end-to-end GWAS execution additionally requires appropriate genotype data and the external PLINK/GEMMA executables.

Run local tests with:

```r
testthat::test_local()
```

Run an R package check with:

```bash
R CMD check .
```

## Citation

Please cite **baGWAS** and the external methods/software used by your analysis. The methodological and software references are listed in [`REFERENCES.md`](REFERENCES.md).

For reproducibility, record at least:

- baGWAS version
- R version
- PLINK version
- GEMMA version
- EMMREML version
- ggmanh version

You can record the R package versions with:

```r
sessionInfo()
packageVersion("baGWAS")
packageVersion("EMMREML")
packageVersion("ggmanh")
```

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening an issue or pull request. Please do not upload private genotype or phenotype data to GitHub.

## License

baGWAS is released under the MIT License. See [`LICENSE`](LICENSE).

## Repository maintenance

The `.github/workflows/R-CMD-check.yaml` workflow is configured to run automated R package checks on GitHub Actions. The project also includes issue templates, contribution guidance, citation metadata, and a security policy.

## GitHub installation example

Once the repository is public, users can install the development version with:

```r
remotes::install_github("YOUR_GITHUB_USERNAME/baGWAS")
```

Replace `YOUR_GITHUB_USERNAME` with the actual GitHub account that owns the repository.
