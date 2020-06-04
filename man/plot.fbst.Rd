\name{plot.fbst}
\alias{plot.fbst}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.fbst
}
\description{
Plots the results of a Full Bayesian Significance Test.
}
\usage{
\method{plot}{fbst}(x, ..., leftBoundary = -100, rightBoundary = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An Object of class \code{"fbst"}.
  }
  \item{...}{
    Additional parameters, see \code{"plot(x, ...)"}.
  }
  \item{leftBoundary}{
    x-coordinate for the left boundary to which is used for visualising the evidence in support of the null hypothesis. Defaults to -100.
  }
  \item{rightBoundary}{
    x-coordinate for the right boundary to which is used for visualising the evidence in support of the null hypothesis. Defaults to 100.
  }
}
\details{
 Plots the surprise function, the supremum of the surprise function restricted to the null set (blue point) and visualises the Bayesian e-value against the sharp null hypothesis as the blue shaded area under the surprise function. The Bayesian e-value in favour of the sharp null hypothesis is visualised as the red shaded area under the surprise function.
}
\value{
  Returns a plot.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
For a details, see: https://arxiv.org/abs/2001.10577 and https://arxiv.org/pdf/2001.10577.pdf.
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
dimensionTheta = 3, dimensionNullset = 2)
plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
