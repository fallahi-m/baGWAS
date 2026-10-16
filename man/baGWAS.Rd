
\name{baGWAS}
\alias{baGWAS}
\title{breed-average Genome-Wide Association Analysis}
\description{
Performs inter-breed genome-wide association analysis using genetic
data and breed-average phenotype values. The package implements two
complementary analytical strategies, \code{genoGWAS} and
\code{freqGWAS}, both of which use breed-average phenotype values
to investigate genetic associations across breeds.

The two strategies differ in how genetic variation is represented
and analyzed:

\enumerate{
  \item \code{genoGWAS}: Uses individual-level genotype data, with
  breed-average phenotype values assigned to genotyped individuals
  according to their breed. Association analysis is performed using
  linear mixed models implemented in GEMMA.

  \item \code{freqGWAS}: Uses breed-level allele frequencies and
  breed-average phenotype values, together with a breed-level
  relationship matrix derived from identity-by-state (IBS)
  information. Association analysis is performed using EMMREML.
}
}
\usage{
baGWAS(strategy = "", organism = "", pheno,
       covariate = NULL, trait = "", geno = "",
       group_col = "FID", output_dir = ".")
}
\arguments{
  \item{strategy}{
  Character string specifying the analytical strategy.
  Accepted values are \code{"genoGWAS"} and \code{"freqGWAS"}.
  Both strategies use breed-average phenotype values, but differ
  in their representation of genetic variation and association
  analysis methods.
  }

  \item{organism}{
  Character string specifying the organism. Supported keywords in
  the current implementation include \code{"dog"}, \code{"sheep"},
  \code{"cattle"}, \code{"cow"}, \code{"horse"}, \code{"mouse"},
  and \code{"rice"}. The selected keyword determines the
  organism-specific options passed to PLINK.
  }

  \item{pheno}{
  Data frame containing the identifiers, breed or group information,
  and numeric phenotype values required by the selected strategy.
  The column specified by \code{trait} must be present. Phenotype
  values represent breed-level averages rather than necessarily
  measured individual-level phenotypes. Required identifier and
  grouping columns depend on the selected strategy and the
  implementation's input requirements.
  }

  \item{covariate}{
  Optional data frame containing covariates for inclusion in the
  association model. If \code{NULL}, no additional covariates are
  supplied. Covariate identifiers and values must be compatible
  with the input data and the selected strategy.
  }

  \item{trait}{
  Character string specifying the name of the numeric phenotype
  column in \code{pheno}.
  }

  \item{geno}{
  Character string specifying the prefix of a PLINK binary genotype
  dataset. The corresponding \code{.bed}, \code{.bim}, and
  \code{.fam} files must be available. Provide the dataset prefix
  in the format expected by the selected strategy.
  }

  \item{group_col}{
  Character string specifying the column used to define breeds,
  populations, or other groups for the frequency-based strategy.
  Defaults to \code{"FID"}.
  }

  \item{output_dir}{
  Character string specifying the directory in which output files
  are written. Defaults to the current working directory.
  }
}
\details{
\strong{Breed-average phenotype framework}

The package implements inter-breed GWAS in which phenotypic
information is represented by breed-level average values. Depending
on the selected strategy, these values are assigned to genotyped
individuals or analyzed alongside breed-level genetic summaries.
Consequently, phenotype values do not necessarily represent
measurements obtained from the same individuals for whom genotype
data are available.

\strong{1. genoGWAS}

The \code{genoGWAS} strategy uses individual-level genotype data
and assigns breed-average phenotype values to genotyped individuals
according to their breed. PLINK is used for genotype-data
preparation, and GEMMA is used to perform association analysis
with linear mixed models.

\strong{2. freqGWAS}

The \code{freqGWAS} strategy represents genetic variation using
breed-level allele frequencies. Breed-average phenotype values
are analyzed alongside these genetic summaries. A relationship
matrix derived from identity-by-state (IBS) information represents
genetic relationships among breeds or groups, and EMMREML is used
for association analysis.

\strong{External software and dependencies}

PLINK and GEMMA are external command-line programs that must be
installed separately on the operating system and made accessible
to the package. GEMMA is required for \code{genoGWAS}. PLINK is
used for genotype-data preparation and genetic summaries as
required by the selected strategy. The \code{freqGWAS} strategy
also uses the EMMREML R package.

\strong{Interpretation of results}

Breed-average phenotypes are not equivalent to independently
measured individual-level phenotypes. In \code{genoGWAS},
individuals from the same breed may receive the same phenotype
value. Therefore, the number of genotyped individuals should not
automatically be interpreted as the number of independent
phenotype observations.

In \code{freqGWAS}, the analysis uses breed-level genetic summaries
and corresponding breed-average phenotype values. Users should
consider breed structure, genetic relatedness, population
stratification, within-breed phenotype variability, model
assumptions, and multiple testing when interpreting association
results.

The two strategies use different genetic data representations
and statistical models. Their results should not be assumed to
be directly equivalent or interchangeable.
}
\value{
Returns a standardized GWAS result data frame. Documented result
columns include \code{CHR}, \code{SNP}, \code{BP}, and \code{P}.
The function may also generate result files, diagnostic plots,
and a run log according to the selected strategy and package
implementation.
}
\section{Output files}{
The analysis may generate the following files in the output directory:

\itemize{
  \item \code{baGWAS.Res}: association results.
  \item \code{manhattan.png}: Manhattan plot of association results.
  \item \code{qqplot.png}: quantile-quantile (Q-Q) plot.
  \item \code{baGWAS.log}: run log.
  \item \code{SNPs_with_errors.txt}: additional file that may be
  generated by \code{freqGWAS} when SNP models fail.
}

The files generated depend on the selected strategy and the
analysis outcome.
}
\examples{
\dontrun{
# 1. Individual-level genotypes with breed-average phenotypes
result_geno <- baGWAS(
  strategy = "genoGWAS",
  organism = "dog",
  pheno = pheno,
  covariate = NULL,
  trait = "weight",
  geno = "/data/geno"
)

# 2. Breed-level allele frequencies with breed-average phenotypes
result_freq <- baGWAS(
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
