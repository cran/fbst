\name{show.fbst}
\alias{show.fbst}
\alias{show}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
show.fbst
}
\description{
Prints the main results of a Full Bayesian Significance Test to the console.
}
\usage{
\method{show}{fbst}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An Object of class \code{"fbst"}.
  }
}
\details{
Shows the main results of a Full Bayesian Significance Test stored in an object of class \code{fbst}.
}
\value{
  Prints the results onto the console.
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
show(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
