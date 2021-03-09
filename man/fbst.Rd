\name{fbst}
\alias{fbst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fbst
}
\description{
The function computes the Full Bayesian Significance Test (FBST) and the e-value, which is the Bayesian evidence against a sharp null hypothesis. The function assumes posterior MCMC draws and constructs a posterior density based on a kernel density estimator subsequently.
}
\usage{
fbst(posteriorDensityDraws, nullHypothesisValue, FUN, par, 
dimensionTheta, dimensionNullset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{posteriorDensityDraws}{
    Vector of MCMC posterior parameter draws
  }
  \item{nullHypothesisValue}{
    Parameter value belonging to the sharp null hypothesis
  }
  \item{FUN}{
    Reference function
  }
  \item{par}{
    Additional parameters of the reference function 
  }
  \item{dimensionTheta}{
    Dimension of the parameter space
  }
  \item{dimensionNullset}{
    Dimension of the null set corresponding to the null hypothesis
  }
}
\details{
If no reference function is specified, a flat reference function \eqn{r(\theta)=1} is used as default reference function.
}
\value{
  Returns an object of class fbst.
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
grp2=rnorm(50,0.8,3.2)

p = as.vector(BayesFactor::ttestBF(x=grp1,y=grp2, 
  posterior = TRUE, iterations = 3000, 
  rscale = "medium")[,4])

# flat reference function
res = fbst(posteriorDensityDraws = p, nullHypothesisValue = 0, 
dimensionTheta = 2, dimensionNullset = 1)
summary(res)
plot(res)

# medium Cauchy C(0,1) reference function
res_med = fbst(posteriorDensityDraws = p, nullHypothesisValue = 0, dimensionTheta = 2, 
dimensionNullset = 1, FUN = dcauchy, par = list(location = 0, scale = sqrt(2)/2))
summary(res_med)
plot(res_med)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.