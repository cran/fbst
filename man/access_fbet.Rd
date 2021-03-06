\name{access_fbet-method}
\alias{$,fbet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\docType{methods}
\title{
 Returns an object from an object of class \code{fbet}.
}
\description{
  Returns an object from an object of class \code{fbet}
}    
%- maybe also 'usage' for other objects documented here.
\details{
 -
}
\value{
  -
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
\examples{
set.seed(57)
grp1=rnorm(50,0,1.5)
grp2=rnorm(50,0.8,3.2)

p = as.vector(BayesFactor::ttestBF(x=grp1,y=grp2, 
  posterior = TRUE, iterations = 3000, 
  rscale = "medium")[,4])

# flat reference function
res = fbet(p, interval = c(-0.1,0.1), nu=1, FUN=NULL, par=NULL)

# Return the Bayesian evidence value for the interval null hypothesis
res$eValueH0
}
