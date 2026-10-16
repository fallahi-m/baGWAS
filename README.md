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

baGWAS is designed for studies in which new genotyping is not feasible. It repurposes publicly available genotype data from purebred animals and combines them with breed-average phenotypes (e.g., from breed standards or the literature) to identify loci associated with phenotypic differences between breeds.

### Requirements

- **Purebred animals only:** Breed identity should represent a distinct genetic background.
- **Adequate sample size per breed:** As a practical guideline, aim for at least 10–25 genotyped animals per breed to improve the stability of breed-level allele-frequency estimates and relatedness estimates. This is a guideline, not a guarantee of statistical adequacy.
- **Inter-breed objective:** The analysis targets genetic associations with phenotypic differences between breeds, rather than within-breed variation.

Because each individual is assigned its breed-average phenotype, within-breed phenotypic variance is zero by construction. The association signal therefore primarily reflects between-breed differences in allele frequencies and breed-average phenotypes. Results should not be interpreted as direct evidence of within-breed genetic effects.

### Strategies

Both strategies test whether genetic variation covaries with breed-average phenotypes while accounting for population structure. They differ in the representation of genotype data and the unit of analysis.

| Feature | `genoGWAS` | `freqGWAS` |
|---|---|---|
| Input | Individual-level genotypes | Breed-level allele frequencies |
| Unit of analysis | Individuals assigned breed-average phenotypes | Breeds |
| Statistical model | Linear mixed model (GEMMA) | Mixed model (EMMREML) |
| Relatedness matrix | Individual-level genomic relationship matrix | Breed-level IBS-derived relationship matrix |
| Phenotypic replication | Individuals are the nominal observations, but independent phenotypic information is primarily determined by the number of breeds | Number of breeds |
| Computational cost | Generally higher | Generally lower when the number of breeds is much smaller than the number of individuals |

### Considerations

- **Pseudo-replication:** In `genoGWAS`, individuals within a breed share the same assigned phenotype. Consequently, the nominal number of observations can exceed the number of independent breed-level phenotypic observations. If this dependence is not adequately accounted for, test statistics may be inflated. A mixed model does not automatically guarantee correct statistical calibration for this data structure.

- **Weighting:** The contribution of each breed to the `genoGWAS` analysis may depend on its sample size and the model implementation. In `freqGWAS`, each breed is represented by one allele-frequency observation, but this does not necessarily mean that all breeds receive equal statistical weight. The actual weighting depends on the model and its variance assumptions.

- **Allele-frequency precision:** In `freqGWAS`, the precision of breed-level allele-frequency estimates depends on the number of genotyped animals and genotype quality within each breed. Sampling uncertainty may not be explicitly modelled, particularly for breeds with small sample sizes.

- **Population structure:** Both strategies can be affected by population structure, uneven breed representation, and correlations between genetic ancestry and breed-average phenotypes. Relatedness adjustment does not eliminate every potential source of confounding.

### Selection Criteria

- **Use `freqGWAS`** when breed-level allele frequencies are already available, computational resources are limited, or a breed-level analysis is preferred. This approach aligns the unit of analysis more directly with the breed-average phenotype, but it does not automatically guarantee better-calibrated P-values.

- **Use `genoGWAS`** when individual-level genotype data are available and individual-level genotype quality control or covariate handling is important. Because the phenotype is repeated within breeds, model assumptions and statistical calibration should be evaluated carefully.

- **When both data types are available, consider running both strategies.** Compare the genomic inflation factor (\(\lambda_{GC}\)), the correlation of \(-\log_{10}(P)\) values, and the overlap of significant SNPs. Where possible, also assess type I error and statistical power using simulations that reflect the study's breed structure and sample sizes.

Associations supported by both strategies may be useful candidates for follow-up. However, agreement between strategies alone does not establish statistical validity; appropriate calibration and independent replication remain important.

> **Note:** `genoGWAS` and `freqGWAS` are not statistically equivalent. Their P-values may be correlated but are not expected to be identical because the strategies differ in their genotype representation, unit of analysis, relatedness matrices, and model assumptions.
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
