\name{names.fbet}
\alias{names.fbet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
names.fbet
}
\description{
Plots the names of the objects stored in the \code{data} object of a Full Bayesian Evidence Test.
}
\usage{
\method{names}{fbet}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An Object of class \code{"fbet"}.
  }
}
\details{
Plots the names of the objects stored in the \code{data} object of a Full Bayesian Evidence Test.
}
\value{
  Returns a list of names.
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
res = fbet(p, interval = c(-0.1,0.1), nu=1, FUN=NULL, par=NULL)
names(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
