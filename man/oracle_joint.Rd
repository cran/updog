% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{oracle_joint}
\alias{oracle_joint}
\title{The joint probability of the genotype and the genotype estimate
of an oracle estimator.}
\usage{
oracle_joint(n, ploidy, seq, bias, od, dist)
}
\arguments{
\item{n}{The read-depth.}

\item{ploidy}{The ploidy of the individual.}

\item{seq}{The sequencing error rate.}

\item{bias}{The allele-bias.}

\item{od}{The overdispersion parameter.}

\item{dist}{The distribution of the alleles.}
}
\value{
A matrix. Element (i, j) is the joint probability of estimating
    the genotype to be i+1 when the true genotype is j+1. That is, the
    estimated genotype indexes the rows and the true genotype indexes
    the columns. This is when
    using an oracle estimator.
}
\description{
This returns the joint distribution of the true genotypes and an oracle
estimator given perfect knowledge of the data generating process. This is a useful
approximation when you have a lot of individuals.
}
\details{
To come up with \code{dist}, you need some additional assumptions.
For example, if the population is in Hardy-Weinberg equilibrium and
the allele frequency is \code{alpha} then you could calculate
\code{dist} using the R code: \code{dbinom(x = 0:ploidy, size = ploidy, prob = alpha)}.
Alternatively, if you know the genotypes of the individual's two parents are, say,
\code{ref_count1} and \code{ref_count2}, then you could use the \code{\link[updog]{get_q_array}}
function from the updog package: \code{get_q_array(ploidy)[ref_count1 + 1, ref_count2 + 1, ]}.

See the Examples to see how to reconcile the output of \code{oracle_joint}
with \code{\link{oracle_mis}} and \code{\link{oracle_mis_vec}}.
}
\examples{
## Hardy-Weinberg population with allele-frequency of 0.75.
## Moderate bias and moderate overdispersion.
ploidy <- 4
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
                   bias = 0.7, od = 0.01, dist = dist)
jd

## Get same output as oracle_mis this way:
1 - sum(diag(jd))
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)

## Get same output as oracle_mis_vec this way:
1 - diag(sweep(x = jd, MARGIN = 2, STATS = colSums(jd), FUN = "/"))
oracle_mis_vec(n = 100, ploidy = ploidy, seq = 0.001,
               bias = 0.7, od = 0.01, dist = dist)

}
\references{
\itemize{
  \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. \doi{10.1534/genetics.118.301468}.}
}
}
\seealso{
\describe{
  \item{\code{\link{oracle_plot}}}{For visualizing the joint distribution output from \code{oracle_joint}.}
  \item{\code{\link{oracle_mis_from_joint}}}{For obtaining the same results as \code{\link{oracle_mis}}
      directly from the output of \code{oracle_joint}.}
  \item{\code{\link{oracle_mis_vec_from_joint}}}{For obtaining the same results as \code{\link{oracle_mis_vec}}
      directly from the output of \code{oracle_joint}.}
  \item{\code{\link{oracle_cor_from_joint}}}{For obtaining the same results as \code{\link{oracle_cor}}
      directly from the output of \code{oracle_joint}.}
}
}
\author{
David Gerard
}
