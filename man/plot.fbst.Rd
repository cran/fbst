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
\method{plot}{fbst}(x, ..., leftBoundary = -100, rightBoundary = 100, type = "contour", parNames = NULL, 
xlimleft = NULL, xlimright = NULL, xlabString = "Parameter", ylabString = NULL)
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
  \item{type}{
    Relevant only if dim=2. Defaults to "contour" which provides a contour plot of the posterior, with a magenta point that shows the supremum over the null set. Alternatively, "persp" provides a 3-dimensional perspective plot of the posterior.
  }
  \item{parNames}{
    Vector of two entries which specifies the names for the parameters. Only relevant if dimensionTheta=2.
  }
  \item{xlimleft}{
    The left value for the x-axis range for the plot. Defaults to the minimum value provided in the posterior draws stored in the FBST object.  
  }
  \item{xlimright}{
    The right value for the x-axis range for the plot. Defaults to the maximum value provided in the posterior draws stored in the FBST object.
  }
  \item{xlabString}{
    String for the x-axis label. Defaults to "Parameter".
  }
  \item{ylabString}{
    String for the y-axis label. Default to "density".
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
dimensionTheta = 2, dimensionNullset = 1)
plot(res)
plot(res, xlimleft = -1.5, xlimright = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
