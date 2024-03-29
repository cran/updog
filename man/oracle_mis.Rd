% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{oracle_mis}
\alias{oracle_mis}
\title{Calculate oracle misclassification error rate.}
\usage{
oracle_mis(n, ploidy, seq, bias, od, dist)
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
A double. The oracle misclassification error rate.
}
\description{
Given perfect knowledge of the data generating parameters,
\code{oracle_mis} calculates the misclassification error
rate, where the error rate is taken over both the data generation
process and the allele-distribution.
This is an ideal level of the misclassification error rate and
any real method will have a larger rate than this. This is a useful
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
}
\examples{
## Hardy-Weinberg population with allele-frequency of 0.75.
## Moderate bias and moderate overdispersion.
## See how oracle misclassification error rates change as we
## increase the ploidy.
ploidy <- 2
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)

ploidy <- 4
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)

ploidy <- 6
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)

}
\references{
\itemize{
  \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. \doi{10.1534/genetics.118.301468}.}
}
}
\author{
David Gerard
}
