# baGWAS
### An Inter-Breed GWAS Framework Integrating Individual-Level Genotype and Breed-Average Phenotypes

**baGWAS** is an R package for inter-breed genome-wide association analysis using individual-level genotype data and breed-average phenotypes through two complementary workflows:

- **`genoGWAS`** — retains individual-level genotypes and performs linear mixed-model association analysis using GEMMA.
- **`freqGWAS`** — summarizes genotypes as breed-level allele frequencies and performs mixed-model association analysis using EMMREML and an IBS-derived breed relationship matrix.

Both workflows return a standardized GWAS result table, diagnostic plots, and a run log.

> **Version:** 1.0.0 (development release)
> **Platform:** Unix-like systems; installation instructions below are for Ubuntu Linux.

## Installation

## 1. Install external software 

Open the Ubuntu Terminal and install PLINK 1.9 and GEMMA:
```r
sudo apt update
sudo apt install plink1.9 gemma
```

- genoGWAS — requires PLINK for genotype-data preparation and GEMMA for linear mixed-model association analysis.
- freqGWAS — requires PLINK for breed-level allele-frequency and IBS calculations.

> PLINK 1.9 is the recommended starting point for the current implementation. 
- These external programs must be installed separately and accessible through the system `PATH` or the supported environment-variable configuration.

## 2. Install R dependencies

The following R packages are required by baGWAS:

- `data.table`
- `EMMREML`
- `dplyr`
- `ggplot2`
- `ggmanh`


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

### 3. Install baGWAS 

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("fallahi-m/baGWAS")
```

Load the package and check external software availability:

```r
library(baGWAS)
baGWAS_check_tools()
```

## Quick start
The main function is baGWAS(). The geno argument specifies the PLINK binary dataset prefix, and trait identifies the numeric phenotype column.

### 1. genoGWAS 

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

### 2. freqGWAS 

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

The covariate argument can be NULL, in which case the model uses an intercept only. The optional threads argument controls parallel processing; by default, the package uses all available CPU cores except one.

## Input data

###  Genotype data

`geno` is a PLINK binary dataset prefix. The following files must exist:

```text
geno.bed
geno.bim
geno.fam
```

For example, use geno = "/home/user/data/geno" rather than a path to an individual file.

###  Phenotype data

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


###  Covariates
`Covariates` must be  NULL or a data frame containing:

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

## Choosing a Strategy

The two strategies differ in how genetic variation is represented and analysed.
| Feature | `genoGWAS` | `freqGWAS` |
|---|---|---|
| Genotype representation | Individual-level genotypes | Breed-level allele frequencies |
| Phenotype representation | Breed-average phenotype assigned to each individual within the breed | Breed-average phenotype represented at the breed level |
| Unit of association analysis | Individual genotype records with breed-average phenotypes | Breeds |
| Statistical model | Linear mixed model implemented in GEMMA | Mixed-model analysis implemented using EMMREML |
| Relationship matrix | Individual-level genomic relationship matrix | Breed-level IBS-derived relationship matrix |
| Computational performance | Substantially faster execution in the developer's real-data benchmarks | Slower execution in the same benchmarks |
| Main analytical emphasis | Association testing that retains individual genotype information | Direct analysis of associations between breed-level allele frequencies and breed-average phenotypes |

*Runtime comparisons reflect the developer's empirical tests; actual performance depends on dataset size, marker density, hardware, and software configuration..*


- **Use `genoGWAS`** when retaining individual-level genotype information is important to the analytical design. This strategy evaluates individual genotype patterns in relation to assigned breed-average phenotypes while accounting for genetic relationships among animals. It is also the preferred starting point when computational speed is an important consideration, given the faster execution observed in the developer's real-data benchmarks.

- **Use `freqGWAS`** when the primary objective is to investigate associations between breed-level allele-frequency differences and breed-average phenotypic differences directly. This strategy summarizes genotypes within breeds before association testing and provides an explicit population-level representation of genetic variation.

- **Consider both strategies** when assessing the consistency of findings across analytical approaches. Compare effect directions, P-values, Manhattan and Q–Q plots, genomic inflation factors ($\lambda_{GC}$), and overlap among significant markers. Concordant findings provide complementary evidence of robustness but do not constitute independent replication because the strategies use the same underlying populations and breed-average phenotypic information.


### Study Design Recommendations

- **Purebred animals:** Use genetically representative, high-quality genotype data from animals with reliable breed assignments. Mixed-breed animals should not be assigned breed-average phenotypes as if they represented a single purebred population.

- **Animals per breed:** As a practical design recommendation, aim for **at least 20 unrelated, high-quality genotyped animals per breed**, with 25–30 or more preferred when feasible. This recommendation is not a universally validated minimum for GWAS. Published canine breed-mapping studies have used up to 25 animals per breed, demonstrating the feasibility of this sampling scale ([Hayward et al., 2016](https://doi.org/10.1038/ncomms10460); [Hayward et al., 2019](https://doi.org/10.1371/journal.pgen.1008003)). Larger within-breed samples can improve the precision of allele-frequency estimates and the characterization of genetic relationships, although the benefits depend on breed-specific diversity and relatedness.

- **Number of breeds:** Aim to include **at least 80 genetically informative breeds whenever feasible**. Approximately 70 or more may be a practical starting point when data availability is limited, but neither threshold guarantees adequate statistical power. The cited studies demonstrate the feasibility of analyses involving 80–175 breeds; they do not establish a universal minimum ([Boyko et al., 2010](https://doi.org/10.1371/journal.pbio.1000451); [Hayward et al., 2016](https://doi.org/10.1038/ncomms10460); [Hayward et al., 2019](https://doi.org/10.1371/journal.pgen.1008003)). 

- **Phenotype coverage:** Include breeds spanning a broad range of phenotypic values and use consistent phenotype definitions and reliable sources across breeds.
- The number and diversity of breeds are particularly important because animals within a breed share the same breed-average phenotype in genoGWAS. Increasing the number of animals per breed improves genetic characterization but does not add distinct breed-level phenotype observations.


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

- **baGWAS.Res** — standardized GWAS results containing CHR, SNP, BP, and P.

- **manhattan.png** — Manhattan plot of association P-values.

- **qqplot.png** — Q–Q plot with a confidence band and genomic inflation statistic.

- **baGWAS.log** — records the analysis strategy, organism, trait, sample and SNP counts, runtime, and diagnostic information.

 The result table is also returned by baGWAS() and assigned to GWAS.res in the examples above. By default, output_dir is the current working directory, and the package creates the Result subdirectory there. Intermediate files are removed after execution.

## Statistical Considerations and Interpretation

- **Population structure:** Shared ancestry and genetic similarities among breeds can produce associations between markers and breed-average phenotypes. Relationship matrices help account for genetic structure, but residual confounding may remain.

- **Allele-frequency precision:** Estimates in freqGWAS depend on sample size, genotype quality, allele frequency, and relatedness among sampled animals. Sampling uncertainty may not be explicitly incorporated into the association model.

- **Model diagnostics:** Interpret P-values alongside Manhattan and Q–Q plots, genomic inflation factors, and sensitivity analyses. Where feasible, simulations under an appropriate null model can help assess calibration.

> **Interpretation note:** `baGWAS` targets genetic associations with phenotypic differences at the **breed level**, rather than within breeds. Identified markers are not direct evidence of within-breed genotype–phenotype effects or causality. Associations may also reflect population history, shared ancestry, or other breed-associated factors. Findings should be interpreted in light of the study design, phenotype reliability, genotype quality, and adequacy of the statistical model..

## Testing

Run the package tests from the repository root:

```r
testthat::test_local()
```

Run a package check from the terminal:

```bash
R CMD check .
```
Full end-to-end GWAS execution additionally requires suitable genotype data and the external executables required by the selected strategy.

## Citation

Please cite **baGWAS** and the methods and software used in your analysis. Relevant references are listed in [`REFERENCES.md`](REFERENCES.md).

For reproducibility, record the versions of baGWAS, R, PLINK, GEMMA, EMMREML, and ggmanh. R package versions can be recorded using:

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
