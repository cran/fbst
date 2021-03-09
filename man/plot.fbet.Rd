\name{plot.fbet}
\alias{plot.fbet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.fbet
}
\description{
Plots the results of a Full Bayesian Evidence Test.
}
\usage{
\method{plot}{fbet}(x, ..., leftBoundary = -100, rightBoundary = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An Object of class \code{"fbet"}.
  }
  \item{...}{
    Additional parameters, see \code{"plot(x, ...)"}.
  }
  \item{leftBoundary}{
    x-coordinate for the left boundary to which is used for visualising the results. Defaults to -100.
  }
  \item{rightBoundary}{
    x-coordinate for the right boundary to which is used for visualising the results. Defaults to 100.
  }
}
\details{
 Plots the resulting surprise function, the interval null hypothesis (dotted blue lines), the resulting Bayesian evidence interval (solid blue lines), the evidence-threshold \eqn{\nu} (dotted black line) and the resulting Bayesian evidence values. The Bayesian evidence value in favour of the interval null hypothesis is visualized as the blue area, and the Bayesian evidence value in favour of the alternative hypothesis is visualized as the red area.
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
grp2=rnorm(50,0.3,3.2)

p = as.vector(BayesFactor::ttestBF(x=grp1,y=grp2, 
  posterior = TRUE, iterations = 3000, 
  rscale = "medium")[,4])

# flat reference function
res = fbet(p, interval = c(-0.1,0.1), nu=1, FUN=NULL, par=NULL)
summary(res)
plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
