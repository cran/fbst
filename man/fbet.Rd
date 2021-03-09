\name{fbet}
\alias{fbet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fbet
}
\description{
The function computes the Full Bayesian Evidence Test (FBST) and the Bayesian evidence value (the generalized e-value which obtains the e-value of the FBST as a special case), which is the Bayesian evidence against an interval null hypothesis. The function assumes posterior MCMC draws and constructs a posterior density based on a kernel density estimator subsequently. The Bayesian evidence interval is computed using a linear search based on the evidence-threshold and the calculation of the Bayesian evidence value is performed using numerical integration.
}
\usage{
fbet(posteriorDensityDraws, interval, nu=1, FUN=NULL, par=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{posteriorDensityDraws}{
    Vector of MCMC posterior parameter draws
  }
  \item{interval}{
    Vector of two numerical values containing the boundaries of the interval null hypothesis to be tested
  }
  \item{nu}{
    Numerical value which provides the evidence-threshold based on which the Bayesian evidence interval is calculated
  }
  \item{FUN}{
    Reference function
  }
  \item{par}{
    Additional parameters of the reference function 
  }
}
\details{
If no reference function is specified, a flat reference function \eqn{r(\theta)=1} is used as default reference function.
}
\value{
  Returns an object of class fbet.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
For a details, see: https://arxiv.org/abs/2001.10577.
}
\author{
  Riko Kelter
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(57)
grp1=rnorm(50,0,1.5)
grp2=rnorm(50,0.3,3.2)

p = as.vector(BayesFactor::ttestBF(x=grp1,y=grp2, 
  posterior = TRUE, iterations = 3000, 
  rscale = "medium")[,4])

# flat reference function, nu = 0
res = fbet(p, interval = c(-0.1,0.1), nu=0, FUN=NULL, par=NULL)
summary(res)
plot(res)

# flat reference function, nu = 1
res = fbet(p, interval = c(-0.1,0.1), nu=1, FUN=NULL, par=NULL)
summary(res)
plot(res)

# medium Cauchy C(0,1) reference function, nu = 1
res_med = fbet(p = p, interval = c(-0.1,0.1), nu = 1, 
FUN = dcauchy, par = list(location = 0, scale=sqrt(2)))
summary(res_med)
plot(res_med)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.