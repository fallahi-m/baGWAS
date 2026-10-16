
\name{baGWAS}
\alias{baGWAS}
\title{Breed-Average Genome-Wide Association Analysis}
\description{
Performs inter-breed genome-wide association analysis by integrating
individual-level genotype data with breed-average phenotypes.
Two complementary strategies are available: genotype-data-based
GWAS using GEMMA and allele-frequency-based GWAS using EMMREML.
The function prepares analysis inputs, runs the selected workflow,
generates association results and diagnostic plots, and records
execution information in a log file.
}
\usage{
baGWAS(strategy = "", organism = "", pheno,
       covariate = NULL, trait = "", geno = "",
       group_col = "FID", output_dir = ".")
}
\arguments{
  \item{strategy}{
  Character string specifying the analysis strategy.
  Use \code{"genoGWAS"} for breed-level genotype-data-based GWAS,
  or \code{"freqGWAS"} for breed-level allele-frequency-based GWAS.
  }

  \item{organism}{
  Character string identifying the organism. Supported keywords
  in the current implementation include \code{"dog"},
  \code{"sheep"}, \code{"cattle"}, \code{"cow"},
  \code{"horse"}, \code{"mouse"}, and \code{"rice"}.
  The selected keyword must be compatible with the PLINK options
  implemented by the package.
  }

  \item{pheno}{
  Data frame containing individual identifiers, group or breed
  information, and the numeric phenotype to be analyzed.
  It must contain \code{FID}, \code{IID}, and the column specified
  by \code{trait}. The grouping information must be consistent
  with the selected analysis design.
  }

  \item{covariate}{
  Optional data frame of covariates. Use \code{NULL} for an
  intercept-only model. When provided, the data must contain
  identifiers compatible with the phenotype data and the
  covariate columns required by the selected workflow.
  }

  \item{trait}{
  Character string giving the name of the numeric phenotype
  column in \code{pheno}.
  }

  \item{geno}{
  Character string specifying the PLINK binary genotype dataset
  prefix. The corresponding \code{.bed}, \code{.bim}, and
  \code{.fam} files must exist. A path ending in one of these
  extensions may be normalized to the dataset prefix.
  }

  \item{group_col}{
  Character string naming the column used to define groups,
  breeds, or populations for the frequency-based workflow.
  The default shown here is \code{"FID"}.
  }

  \item{output_dir}{
  Character string specifying the directory in which result
  files are written. The default shown here is the current
  working directory.
  }
}
\details{
The \code{genoGWAS} strategy uses individual-level genotype data
prepared with PLINK and performs linear mixed-model association
analysis with GEMMA.

The \code{freqGWAS} strategy uses breed- or group-level allele
frequencies, an identity-by-state (IBS)-derived group relationship
matrix, and EMMREML for association analysis.

These strategies represent different analytical workflows and
their results should not be treated as interchangeable.

External command-line programs must be installed separately and
available through the system PATH or the environment-variable
configuration supported by the package. The required tools depend
on the selected strategy.

Users should evaluate genotype quality control, phenotype
definition, sample structure, group definitions, covariates,
model assumptions, and multiple-testing procedures before
interpreting association results.
}
\value{
Returns a standardized GWAS result data frame. The documented
result columns include \code{CHR}, \code{SNP}, \code{BP}, and
\code{P}. The function also writes result files, diagnostic plots,
and a run log to the output directory, according to the selected
strategy.
}
\section{Output files}{
The analysis produces a result file named \code{baGWAS.Res},
a Manhattan plot named \code{manhattan.png}, a Q-Q plot named
\code{qqplot.png}, and a log file named \code{baGWAS.log}.
The frequency-based workflow may also produce
\code{SNPs_with_errors.txt} when SNP models fail.
Actual output locations depend on the package implementation
and the selected output directory.
}
\examples{
\dontrun{
# Genotype-data-based GWAS
result <- baGWAS(
  strategy = "genoGWAS",
  organism = "dog",
  pheno = pheno,
  covariate = NULL,
  trait = "weight",
  geno = "/data/geno"
)

# Allele-frequency-based GWAS
result <- baGWAS(
  strategy = "freqGWAS",
  organism = "dog",
  pheno = pheno,
  covariate = NULL,
  trait = "weight",
  geno = "/data/geno",
  group_col = "FID"
)
}
}
\seealso{
\code{\link{citation}},
\url{https://github.com/fallahi-m/baGWAS}
}
\keyword{models}
\keyword{documentation}
