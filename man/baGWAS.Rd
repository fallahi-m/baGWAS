\name{baGWAS}
\alias{baGWAS}
\title{Inter-Breed Genome-Wide Association Analysis with Breed-Average Phenotypes}
\description{
Runs an inter-breed genome-wide association study by integrating genotype data
with breed-level average phenotypes. The function provides two complementary
strategies, \code{genoGWAS} and \code{freqGWAS}, which differ in how genetic
variation is represented and analysed.
}
\usage{
baGWAS(
  strategy = c("genoGWAS", "freqGWAS"),
  organism,
  pheno,
  covariate = NULL,
  trait,
  geno,
  output_dir = getwd(),
  group_col = "FID",
  threads = NULL
)
}
\arguments{
  \item{strategy}{
  Character string selecting the analysis strategy: \code{"genoGWAS"} or
  \code{"freqGWAS"}. The default selects \code{"genoGWAS"}.
  }
  \item{organism}{
  Character string identifying the organism. Supported values in the current
  implementation are \code{"dog"}, \code{"sheep"}, \code{"cattle"},
  \code{"cow"}, \code{"horse"}, \code{"mouse"}, and \code{"rice"}.
  The value determines the species-specific option passed to PLINK.
  }
  \item{pheno}{
  Data frame containing at least \code{FID}, \code{IID}, and the numeric
  phenotype column named by \code{trait}. The identifiers must correspond to
  the PLINK genotype dataset. For \code{genoGWAS}, the phenotype assigned to
  each animal should represent the breed-average phenotype for its breed. For
  \code{freqGWAS}, the column specified by \code{group_col} must also be
  present, and phenotype values are summarized at the group level. An optional
  lowercase column named \code{sex} can be supplied for sex assignment during
  genotype preparation in \code{genoGWAS}.
  }
  \item{covariate}{
  Optional data frame containing covariates for the association model. If
  \code{NULL}, the model uses an intercept only. Otherwise, the data frame
  must contain \code{FID} and \code{IID}, along with the numeric covariate
  columns. For \code{freqGWAS}, it must also contain the column specified by
  \code{group_col}. Identifiers and group assignments must be consistent with
  the phenotype and genotype data. Missing numeric covariate values are
  median-imputed, and covariates with zero variance are excluded by the current
  implementation.
  }
  \item{trait}{
  Character string giving the name of the numeric phenotype column in
  \code{pheno}.
  }
  \item{geno}{
  Character string specifying the prefix of a PLINK binary genotype dataset.
  The corresponding \code{.bed}, \code{.bim}, and \code{.fam} files must
  exist. Supply the prefix without a file extension. A path ending in one of
  these extensions is also normalized to the dataset prefix by the current
  implementation.
  }
  \item{output_dir}{
  Character string specifying the directory in which the package creates the
  \code{Result} subdirectory. Defaults to the current working directory.
  }
  \item{group_col}{
  Character string naming the column in \code{pheno} that defines breeds,
  populations, or other groups for \code{freqGWAS}. Defaults to \code{"FID"}.
  When \code{covariate} is supplied for \code{freqGWAS}, that data frame must
  also contain this grouping column.
  }
  \item{threads}{
  Optional number of CPU threads to use. If \code{NULL}, the package selects
  a default based on the available logical CPU cores, generally using all
  cores except one. The value is converted to an integer and limited to at
  least one thread.
  }
}
\details{
\strong{Overview}

\code{baGWAS} is intended for studies of genetic associations with phenotypic
 differences at the breed level. Breed-average phenotypes may be obtained from
 breed standards, published studies, or other reliable sources, while genotype
 data may come from publicly available datasets. Phenotype values are therefore
 not necessarily measurements made on the same individuals from which genotype
 data were obtained.

\strong{Platform}

The package is intended for Unix-like operating systems, including Linux and
macOS. Installation examples in the README target Ubuntu Linux. Compatibility
with external command-line tools may depend on the operating system and the
available software builds.

\strong{The \code{genoGWAS} strategy}

\code{genoGWAS} retains individual-level genotypes and assigns the phenotype
value supplied for each animal to its genotype record. PLINK prepares the
genotype data, and GEMMA performs association testing using a linear mixed model
and an individual-level genomic relationship matrix. Animals from the same
breed may share the same breed-average phenotype value.

\strong{The \code{freqGWAS} strategy}

\code{freqGWAS} summarizes genotype data as group-level allele frequencies.
Phenotype values are summarized within the groups defined by \code{group_col}.
The current implementation requires at least five groups with usable phenotype
and genotype data. Genetic relationships among groups are represented using a
relationship matrix derived from identity-by-state (IBS) information. PLINK
calculates group-level allele frequencies and IBS information, and EMMREML is
used for mixed-model association analysis.

\strong{External software and dependencies}

PLINK and GEMMA are command-line programs that must be installed separately
and made accessible to the package. GEMMA is required for \code{genoGWAS};
PLINK is used by both strategies. PLINK 1.9 is recommended, particularly for
\code{freqGWAS}, which relies on PLINK 1.x commands for group-level frequency
and IBS calculations. The package searches for supported executable names on
the system \code{PATH}. The environment variables \code{BA_GWAS_PLINK},
\code{BA_GWAS_PLINK1}, and \code{BA_GWAS_GEMMA} can be used to specify
executable locations or command names.

\strong{Interpretation}

The number of genotyped animals in \code{genoGWAS} should not automatically be
interpreted as the number of independent phenotype observations, because
animals from the same breed can share a breed-average phenotype. In
\code{freqGWAS}, association testing uses group-level genetic summaries and
group-level phenotypes. Both approaches target associations with phenotypic
differences at the breed level; their results are not direct evidence of
within-breed genotype--phenotype effects or causality. Population structure,
shared ancestry, phenotype quality, genotype quality, model assumptions, and
multiple testing should be considered when interpreting findings.

}
\value{
Returns a standardized GWAS result data frame. Documented result columns are
\code{CHR}, \code{SNP}, \code{BP}, and \code{P}. The result may also carry
attributes describing the strategy, trait, organism, output directory,
genomic inflation statistic, number of SNPs, and runtime, where available.
}
\section{Output files}{
Files are written under
\code{<output_dir>/Result/<strategy>/}. The strategy-specific directory may
contain:
\itemize{
  \item \code{baGWAS.Res}: tab-delimited association results.
  \item \code{manhattan.png}: Manhattan plot of association P-values.
  \item \code{qqplot.png}: Q--Q plot and genomic inflation diagnostic.
  \item \code{baGWAS.log}: run information, summary counts, timing, and
  diagnostic information.
  \item \code{SNPs_with_errors.txt}: file written by \code{freqGWAS} if one
  or more SNP models fail.
}
The files generated depend on the selected strategy and analysis outcome.
}
\section{External tool check}{
Use \code{baGWAS_check_tools()} to check whether supported PLINK and GEMMA
executables can be found and run on the system.
}
\examples{
\dontrun{
# Individual-level genotype analysis with breed-average phenotypes
result_geno <- baGWAS(
  strategy = "genoGWAS",
  organism = "dog",
  pheno = pheno,
  covariate = NULL,
  trait = "weight",
  geno = "/data/geno",
  output_dir = "baGWAS_output"
)

# Group-level allele-frequency analysis
result_freq <- baGWAS(
  strategy = "freqGWAS",
  organism = "dog",
  pheno = pheno,
  covariate = NULL,
  trait = "weight",
  geno = "/data/geno",
  group_col = "FID",
  output_dir = "baGWAS_output"
)
}
}
\seealso{
\code{\link{baGWAS_check_tools}},
\url{https://github.com/fallahi-m/baGWAS}
}
\keyword{models}
\keyword{utilities}
