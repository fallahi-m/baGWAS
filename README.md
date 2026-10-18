# 2baGWAS
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

`baGWAS` is designed for inter-breed genome-wide association studies when collecting individual-level genotype and phenotype data for every animal is impractical. It integrates publicly available genotype data from purebred animals with breed-average phenotypes obtained from breed standards, published studies, or other reliable sources.

The primary objective is to identify genetic markers associated with **phenotypic differences between breeds**, such as differences in body size, height, or other breed-level traits. `baGWAS` provides two complementary strategies: `genoGWAS`, which retains individual-level genotypes, and `freqGWAS`, which summarizes genotypes as breed-level allele frequencies.

The choice of strategy depends on the desired representation of genetic variation, the analytical objective, and the computational performance observed for the available dataset.

### Study Design and Data Requirements

- **Purebred animals:** Use genetically representative, high-quality genotype data from animals with reliable breed assignments. Mixed-breed animals should not be assigned breed-average phenotypes as if they represented a single purebred population.

- **Sample size per breed:** As a practical design recommendation, aim for **at least 20 unrelated, high-quality genotyped animals per breed**, with 25–30 or more preferred when feasible. This recommendation is not a universally validated minimum for GWAS. Published canine breed-mapping studies have used up to 25 animals per breed, demonstrating the feasibility of this sampling scale ([Hayward et al., 2016](https://doi.org/10.1038/ncomms10460); [Hayward et al., 2019](https://doi.org/10.1371/journal.pgen.1008003)). Larger within-breed samples can improve the precision of allele-frequency estimates and the characterization of genetic relationships, although the benefits depend on breed-specific diversity and relatedness.

- **Number of breeds:** Aim to include **at least 80 genetically informative breeds whenever feasible**. A minimum target of approximately 70 breeds may be a practical starting point when data availability is limited, but neither threshold guarantees adequate statistical power. Previous canine breed-mapping studies have successfully analysed breed-average phenotypes across 80, 158, and 175 breeds ([Boyko et al., 2010](https://doi.org/10.1371/journal.pbio.1000451); [Hayward et al., 2016](https://doi.org/10.1038/ncomms10460); [Hayward et al., 2019](https://doi.org/10.1371/journal.pgen.1008003)). These studies provide useful design precedents rather than evidence of a universal minimum. Statistical power also depends on the range of breed-average phenotypes, allele frequencies, genetic architecture, effect sizes, population structure, and multiple-testing correction.

- **Phenotype diversity and consistency:** Include breeds spanning a broad range of phenotypic values whenever possible. Apply consistent phenotype definitions, measurement criteria, and data sources across breeds. Morphological traits such as adult height and body size may be easier to standardize than complex behavioural traits, which can be influenced by training, environment, and other non-genetic factors.

**Important design consideration:** In `genoGWAS`, individual-level genotypes are retained, but all sampled animals belonging to the same breed receive the same breed-average phenotype. Increasing the number of animals within a breed improves genetic characterization but does not create additional distinct breed-level phenotype observations. Consequently, both the number of animals per breed and the number and diversity of breeds should be considered when designing a study.

### Comparison of Strategies

Both strategies investigate associations between genetic variation and breed-average phenotypes while accounting for genetic relationships through their respective relationship matrices. Their main difference is how genetic variation is represented and analysed.

| Feature | `genoGWAS` | `freqGWAS` |
|---|---|---|
| Genotype representation | Individual-level genotypes | Breed-level allele frequencies |
| Phenotype representation | Breed-average phenotype assigned to each individual within the breed | Breed-average phenotype represented at the breed level |
| Unit of association analysis | Individual genotype records with breed-average phenotypes | Breeds |
| Statistical model | Linear mixed model implemented in GEMMA | Mixed-model analysis implemented using EMMREML |
| Relationship matrix | Individual-level genomic relationship matrix | Breed-level IBS-derived relationship matrix |
| Computational performance | Substantially faster execution in the developer's real-data benchmarks | Slower execution in the same benchmarks |
| Main analytical emphasis | Association testing that retains individual genotype information | Direct analysis of associations between breed-level allele frequencies and breed-average phenotypes |

*Computational performance reflects the developer's empirical tests on real datasets. Actual runtime may vary with dataset size, marker density, hardware, and software configuration.*

### How to Choose

- **Use `genoGWAS`** when retaining individual-level genotype information is important to the analytical design. This strategy evaluates individual genotype patterns in relation to assigned breed-average phenotypes while accounting for genetic relationships among animals. It is also the preferred starting point when computational speed is an important consideration, given the faster execution observed in the developer's real-data benchmarks.

- **Use `freqGWAS`** when the primary objective is to investigate associations between breed-level allele-frequency differences and breed-average phenotypic differences directly. This strategy summarizes genotypes within breeds before association testing and provides an explicit population-level representation of genetic variation.

- **Consider both strategies** when assessing the consistency of findings across analytical approaches. Compare effect directions, P-values, Manhattan and Q–Q plots, genomic inflation factors ($\lambda_{GC}$), and overlap among significant markers. Concordant findings provide complementary evidence of robustness but do not constitute independent replication because the strategies use the same underlying populations and breed-average phenotypic information.

### Statistical Considerations

- **Population structure and shared ancestry:** Genetic similarities among breeds can induce associations between genetic markers and breed-average phenotypes. Both strategies account for genetic relationships using their respective relationship matrices, but residual confounding may remain. Results should be interpreted in the context of the study population and the adequacy of the statistical model.

- **Allele-frequency precision:** In `freqGWAS`, breed-level allele frequencies are estimated from the available animals within each breed. Their precision depends on sample size, genotype quality, allele frequency, and the relatedness of sampled animals. Sampling uncertainty in the estimated allele frequencies may not be explicitly incorporated into the association model.

- **Statistical calibration:** Evaluate the reliability of association results using Manhattan and Q–Q plots, genomic inflation factors ($\lambda_{GC}$), and, where feasible, simulations under an appropriate null model. These diagnostics can identify potential calibration problems but cannot independently establish that a model is well calibrated.

- **Sensitivity analysis:** Examine the influence of individual breeds, differences in breed sample sizes, and alternative phenotype definitions where appropriate. Such analyses can help identify associations that depend heavily on particular breeds or assumptions.

> **Interpretation note:** `baGWAS` targets associations with **between-breed phenotypic differences**. Associations identified by either strategy should not be interpreted as direct evidence of within-breed genotype–phenotype effects or as proof of causality. Genetic differences associated with a trait may also reflect population history, shared ancestry, or other breed-associated factors. The strength of the conclusions depends on the number and diversity of breeds, genotype quality, the reliability of breed-average phenotypes, and the statistical performance of the selected model.
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
