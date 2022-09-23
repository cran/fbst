\name{bdm}
\alias{bdm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bdm
}
\description{
Calculates the Bayesian discrepancy measure for a precise null hypothesis.
}
\usage{
bdm(posteriorDensityDraws, nullHypothesisValue=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{posteriorDensityDraws}{
    Vector of (MCMC) posterior parameter draws.
}
  \item{nullHypothesisValue}{
    Parameter value of the precise null hypothesis. Defaults to zero.
}
}
\details{
The BDM is calculated as \eqn{\delta_H(x):=2\cdot P(\theta \in I_H(x)|x)} where \eqn{I_H(x):=(m,\theta_0)} if \eqn{m<\theta_0}, \eqn{I_H(x):=\{m\}} if \eqn{m=\theta_0} and \eqn{I_H(x):=(\theta_0,m)} if \eqn{m>\theta_0}, where \eqn{m} denotes the posterior median of the parameter \eqn{\theta}, and the null hypothesis specifies \eqn{H_0:\theta=\theta_0}.
}
\value{
  Returns the value \eqn{\delta_H(x)} of the BDM.
}
\references{
  For details, see: https://arxiv.org/abs/2105.13716
}
\author{
  Riko Kelter
}
\examples{
set.seed(57)
grp1=rnorm(50,0,1.5)
grp2=rnorm(50,0.8,3.2)

p = as.vector(BayesFactor::ttestBF(x=grp1,y=grp2, 
  posterior = TRUE, iterations = 3000, 
  rscale = "medium")[,4])

bdm(p,0)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
