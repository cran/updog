% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ashdog.R
\name{flexdog_full}
\alias{flexdog_full}
\title{Flexible genotyping for polyploids from next-generation sequencing data.}
\usage{
flexdog_full(
  refvec,
  sizevec,
  ploidy,
  model = c("norm", "hw", "bb", "s1", "s1pp", "f1", "f1pp", "flex", "uniform", "custom"),
  verbose = TRUE,
  mean_bias = 0,
  var_bias = 0.7^2,
  mean_seq = -4.7,
  var_seq = 1,
  mean_od = -5.5,
  var_od = 0.5^2,
  seq = 0.005,
  bias = 1,
  od = 0.001,
  update_bias = TRUE,
  update_seq = TRUE,
  update_od = TRUE,
  itermax = 200,
  tol = 10^-4,
  fs1_alpha = 10^-3,
  p1ref = NULL,
  p1size = NULL,
  p2ref = NULL,
  p2size = NULL,
  snpname = NULL,
  prior_vec = NULL,
  seq_upper = 0.05
)
}
\arguments{
\item{refvec}{A vector of counts of reads of the reference allele.}

\item{sizevec}{A vector of total counts.}

\item{ploidy}{The ploidy of the species. Assumed to be the same for each
individual.}

\item{model}{What form should the prior (genotype distribution) take?
See Details for possible values.}

\item{verbose}{Should we output more (\code{TRUE}) or less
(\code{FALSE})?}

\item{mean_bias}{The prior mean of the log-bias.}

\item{var_bias}{The prior variance of the log-bias.}

\item{mean_seq}{The prior mean of the logit of the sequencing
error rate.}

\item{var_seq}{The prior variance of the logit of the sequencing
error rate.}

\item{mean_od}{The prior mean of the logit of the overdispersion parameter.}

\item{var_od}{The prior variance of the logit of the overdispersion
parameter.}

\item{seq}{The starting value of the sequencing error rate.}

\item{bias}{The starting value of the bias.}

\item{od}{The starting value of the overdispersion parameter.}

\item{update_bias}{A logical. Should we update \code{bias}
(\code{TRUE}), or not (\code{FALSE})?}

\item{update_seq}{A logical. Should we update \code{seq}
(\code{TRUE}), or not (\code{FALSE})?}

\item{update_od}{A logical. Should we update \code{od}
(\code{TRUE}), or not (\code{FALSE})?}

\item{itermax}{The total number of EM iterations to
run.}

\item{tol}{The tolerance stopping criterion. The EM algorithm will stop
if the difference in the log-likelihoods between two consecutive
iterations is less than \code{tol}.}

\item{fs1_alpha}{The value at which to fix
the mixing proportion for the uniform component
 when \code{model = "f1"}, \code{model = "s1"},
 \code{model = "f1pp"}, or \code{model = "s1pp"}.
 I would recommend some small value such as \code{10^-3}.}

\item{p1ref}{The reference counts for the first parent if
\code{model = "f1"} or \code{model = "f1pp"}, or for
the only parent if \code{model = "s1"} or \code{model = "s1pp"}.}

\item{p1size}{The total counts for the first parent if
\code{model = "f1"} or \code{model = "f1pp"},
or for the only parent if \code{model = "s1"} or \code{model = "s1pp"}.}

\item{p2ref}{The reference counts for the second parent if
\code{model = "f1"} or \code{model = "f1pp"}.}

\item{p2size}{The total counts for the second parent if
\code{model = "f1"} or \code{model = "f1pp"}.}

\item{snpname}{A string. The name of the SNP under consideration.
This is just returned in the \code{input} list for your reference.}

\item{prior_vec}{The pre-specified genotype distribution. Only used if
\code{model = "custom"} and must otherwise be \code{NULL}. If specified,
then it should be a vector of length \code{ploidy + 1} with
non-negative elements that sum to 1.}

\item{seq_upper}{The upper bound on the possible sequencing error rate.
Default is 0.05, but you should adjust this if you have prior knowledge
on the sequencing error rate of the sequencing technology.}
}
\value{
An object of class \code{flexdog}, which consists
    of a list with some or all of the following elements:
\describe{
  \item{\code{bias}}{The estimated bias parameter.}
  \item{\code{seq}}{The estimated sequencing error rate.}
  \item{\code{od}}{The estimated overdispersion parameter.}
  \item{\code{num_iter}}{The number of EM iterations ran. You should
      be wary if this equals \code{itermax}.}
  \item{\code{llike}}{The maximum marginal log-likelihood.}
  \item{\code{postmat}}{A matrix of posterior probabilities of each
      genotype for each individual. The rows index the individuals
      and the columns index the allele dosage.}
  \item{\code{genologlike}}{A matrix of genotype \emph{log}-likelihoods of each
      genotype for each individual. The rows index the individuals
      and the columns index the allele dosage.}
  \item{\code{gene_dist}}{The estimated genotype distribution. The
      \code{i}th element is the proportion of individuals with
      genotype \code{i-1}.}
  \item{\code{par}}{A list of the final estimates of the parameters
      of the genotype distribution. The elements included in \code{par}
      depends on the value of \code{model}:
      \describe{
      \item{\code{model = "norm"}:}{
        \describe{
          \item{\code{mu}:}{The normal mean.}
          \item{\code{sigma}:}{The normal standard deviation (not variance).}
        }
      }
      \item{\code{model = "hw"}:}{
        \describe{
          \item{\code{alpha}:}{The major allele frequency.}
        }
      }
      \item{\code{model = "bb"}:}{
        \describe{
          \item{\code{alpha}:}{The major allele frequency.}
          \item{\code{tau}:}{The overdispersion parameter. See the
              description of \code{rho} in the Details of
              \code{\link{betabinom}()}.}
        }
      }
      \item{\code{model = "s1"}:}{
        \describe{
          \item{\code{pgeno}:}{The allele dosage of the parent.}
          \item{\code{alpha}:}{The mixture proportion of the discrete
              uniform (included and fixed at a small value mostly for
              numerical stability reasons). See the description
              of \code{fs1_alpha} in \code{\link{flexdog_full}()}.}
        }
      }
      \item{\code{model = "f1"}:}{
        \describe{
          \item{\code{p1geno}:}{The allele dosage of the first parent.}
          \item{\code{p2geno}:}{The allele dosage of the second parent.}
          \item{\code{alpha}:}{The mixture proportion of the discrete
              uniform (included and fixed at a small value mostly for
              numerical stability reasons). See the description
              of \code{fs1_alpha} in \code{\link{flexdog_full}()}.}
        }
      }
      \item{\code{model = "s1pp"}:}{
        \describe{
          \item{\code{ell1}:}{The estimated dosage of the parent.}
          \item{\code{tau1}:}{The estimated double reduction parameter
              of the parent. Available if \code{ell1} is \code{1}, \code{2},
              or \code{3}. Identified if \code{ell1} is \code{1} or \code{3}.}
          \item{\code{gamma1}:}{The estimated preferential pairing parameter.
              Available if \code{ell1} is \code{2}. However, it is not
              returned in an identified form.}
          \item{\code{alpha}:}{The mixture proportion of the discrete
              uniform (included and fixed at a small value mostly for
              numerical stability reasons). See the description
              of \code{fs1_alpha} in \code{\link{flexdog_full}()}.}
        }
      }
      \item{\code{model = "f1pp"}:}{
        \describe{
          \item{\code{ell1}:}{The estimated dosage of parent 1.}
          \item{\code{ell2}:}{The estimated dosage of parent 2.}
          \item{\code{tau1}:}{The estimated double reduction parameter
              of parent 1. Available if \code{ell1} is \code{1}, \code{2},
              or \code{3}. Identified if \code{ell1} is \code{1} or \code{3}.}
          \item{\code{tau2}:}{The estimated double reduction parameter
              of parent 2. Available if \code{ell2} is \code{1}, \code{2},
              or \code{3}. Identified if \code{ell2} is \code{1} or \code{3}.}
          \item{\code{gamma1}:}{The estimated preferential pairing parameter
              of parent 1. Available if \code{ell1} is \code{2}. However,
              it is not returned in an identified form.}
          \item{\code{gamma2}:}{The estimated preferential pairing parameter
              of parent 2. Available if \code{ell2} is \code{2}. However,
              it is not returned in an identified form.}
          \item{\code{alpha}:}{The mixture proportion of the discrete
              uniform (included and fixed at a small value mostly for
              numerical stability reasons). See the description
              of \code{fs1_alpha} in \code{\link{flexdog_full}()}.}
        }
      }
      \item{\code{model = "flex"}:}{\code{par} is an empty list.}
      \item{\code{model = "uniform"}:}{\code{par} is an empty list.}
      \item{\code{model = "custom"}:}{\code{par} is an empty list.}
      }}
  \item{\code{geno}}{The posterior mode genotype. These are your
      genotype estimates.}
  \item{\code{maxpostprob}}{The maximum posterior probability. This
      is equivalent to the posterior probability of correctly
      genotyping each individual.}
  \item{\code{postmean}}{The posterior mean genotype. In downstream
      association studies, you might want to consider using these
      estimates.}
  \item{\code{input$refvec}}{The value of \code{refvec} provided by
      the user.}
  \item{\code{input$sizevec}}{The value of \code{sizevec} provided
      by the user.}
  \item{\code{input$ploidy}}{The value of \code{ploidy} provided
      by the user.}
  \item{\code{input$model}}{The value of \code{model} provided by
      the user.}
  \item{\code{input$p1ref}}{The value of \code{p1ref} provided by the user.}
  \item{\code{input$p1size}}{The value of \code{p1size} provided by the user.}
  \item{\code{input$p2ref}}{The value of \code{p2ref} provided by the user.}
  \item{\code{input$p2size}}{The value of \code{p2size} provided by the user.}
  \item{\code{input$snpname}}{The value of \code{snpname} provided by the user.}
  \item{\code{prop_mis}}{The posterior proportion of individuals
      genotyped incorrectly.}
}
}
\description{
Genotype polyploid individuals from next generation
sequencing (NGS) data while assuming the genotype distribution is one of
several forms. \code{\link{flexdog_full}()} does this while accounting for allele bias,
overdispersion, and sequencing error. This function has more
options than \code{\link{flexdog}} and is only meant for expert users.
The method is described in detail in Gerard et. al. (2018) and
Gerard and Ferrão (2020).
}
\details{
Possible values of the genotype distribution (values of \code{model}) are:
\describe{
  \item{\code{"norm"}}{A distribution whose genotype frequencies are proportional
      to the density value of a normal with some mean and some standard deviation.
      Unlike the \code{"bb"} and \code{"hw"} options, this will allow for
      distributions both more and less dispersed than a binomial.
      This seems to be the most robust to violations in modeling assumptions, and so is the
      default. This prior class was developed in Gerard and Ferrão (2020).}
  \item{\code{"hw"}}{A binomial distribution that results from assuming that
      the population is in Hardy-Weinberg equilibrium (HWE). This actually does
      pretty well even when there are minor to moderate deviations from HWE.
      Though it does not perform as well as the `"norm"` option when there
      are severe deviations from HWE.}
  \item{\code{"bb"}}{A beta-binomial distribution. This is an overdispersed
      version of \code{"hw"} and can be derived from a special case of the Balding-Nichols model.}
  \item{\code{"s1"}}{This prior assumes the individuals are
      all full-siblings resulting
      from one generation of selfing. I.e. there is only
      one parent. This model assumes
      a particular type of meiotic behavior: polysomic
      inheritance with
      bivalent, non-preferential pairing.}
  \item{\code{"f1"}}{This prior assumes the individuals are all
      full-siblings resulting
      from one generation of a bi-parental cross.
      This model assumes
      a particular type of meiotic behavior: polysomic
      inheritance with
      bivalent, non-preferential pairing.}
  \item{\code{"f1pp"}}{This prior allows for double reduction
      and preferential pairing in an F1 population of tretraploids.}
  \item{\code{"s1pp"}}{This prior allows for double reduction
      and preferential pairing in an S1 population of tretraploids.}
  \item{\code{"flex"}}{Generically any categorical distribution. Theoretically,
      this works well if you have a lot of individuals. In practice, it seems to
      be much less robust to violations in modeling assumptions.}
  \item{\code{"uniform"}}{A discrete uniform distribution. This should never
      be used in practice.}
  \item{\code{"custom"}}{A pre-specified prior distribution. You specify
      it using the \code{prior_vec} argument. You should almost never
      use this option in practice.}
}

You might think a good default is \code{model = "uniform"} because it is
somehow an "uninformative prior." But it is very informative and tends to
work horribly in practice. The intuition is that it will estimate
the allele bias and sequencing error rates so that the estimated genotypes
are approximately uniform (since we are assuming that they are approximately
uniform). This will usually result in unintuitive genotyping since most
populations don't have a uniform genotype distribution.
I include it as an option only for completeness. Please don't use it.

The value of \code{prop_mis} is a very intuitive measure for
the quality of the SNP. \code{prop_mis} is the posterior
proportion of individuals mis-genotyped. So if you want only SNPS
that accurately genotype, say, 95\% of the individuals, you could
discard all SNPs with a \code{prop_mis} over \code{0.05}.

The value of \code{maxpostprob} is a very intuitive measure
for the quality of the genotype estimate of an individual.
This is the posterior probability of correctly genotyping
the individual when using \code{geno} (the posterior mode)
as the genotype estimate. So if you want to correctly genotype,
say, 95\% of individuals, you could discard all individuals
with a \code{maxpostprob} of under \code{0.95}. However, if you are
just going to impute missing genotypes later, you might consider
not discarding any individuals as \code{flexdog}'s genotype estimates will
probably be more accurate than other more naive approaches, such
as imputing using the grand mean.

In most datasets I've examined, allelic bias is a major issue. However,
you may fit the model assuming no allelic bias by setting
\code{update_bias = FALSE} and \code{bias_init = 1}.

Prior to using \code{flexdog}, during the read-mapping step,
you could try to get rid of allelic bias by
using WASP (\doi{10.1101/011221}). If you are successful
in removing the allelic bias (because its only source was the read-mapping
step), then setting \code{update_bias = FALSE} and \code{bias_init = 1}
would be reasonable. You can visually inspect SNPs for bias by
using \code{\link{plot_geno}()}.

\code{flexdog()}, like most methods, is invariant to which allele you
label as the "reference" and which you label as the "alternative".
That is, if you set \code{refvec} with the number of alternative
read-counts, then the resulting genotype estimates
will be the estimated allele dosage of the alternative allele.
}
\examples{
## A natural population. We will assume a
## normal prior since there are so few
## individuals.
data("uitdewilligen")
ploidy  <- 4
refvec  <- uitdewilligen$refmat[, 1]
sizevec <- uitdewilligen$sizemat[, 1]
fout    <- flexdog_full(refvec  = refvec,
                        sizevec = sizevec,
                        ploidy  = ploidy,
                        model   = "norm")
plot(fout)

}
\references{
\itemize{
  \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. \doi{10.1534/genetics.118.301468}.}
  \item{Gerard, David, and Luís Felipe Ventorim Ferrão. "Priors for genotyping polyploids." Bioinformatics 36, no. 6 (2020): 1795-1800. \doi{10.1093/bioinformatics/btz852}.}
}
}
\seealso{
Run \code{browseVignettes(package = "updog")} in R for example usage.
    Other useful functions include:
\describe{
    \item{\code{\link{multidog}}}{For running \code{flexdog} on multiple
          SNPs in parallel.}
    \item{\code{\link{flexdog}}}{For a more user-friendly version of
          \code{flexdog_full}.}
    \item{\code{\link{rgeno}}}{For simulating genotypes under the allowable
          prior models in \code{flexdog}.}
    \item{\code{\link{rflexdog}}}{For simulating read-counts under the
          assumed likelihood model in \code{flexdog}.}
    \item{\code{\link{plot.flexdog}}}{For plotting the output of
          \code{flexdog}.}
    \item{\code{\link{oracle_mis}}}{For calculating the oracle genotyping
          error rates. This is useful for read-depth calculations
          \emph{before} collecting data. After you have data, using
          the value of \code{prop_mis} is better.}
    \item{\code{\link{oracle_cor}}}{For calculating the correlation
          between the true genotypes and an oracle estimator
          (useful for read-depth calculations \emph{before}
          collecting data).}
}
}
\author{
David Gerard
}
