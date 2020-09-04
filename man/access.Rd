\name{access-method}
\alias{$,fbst-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\docType{methods}
\title{
 Returns an object from an object of class \code{fbst}.
}
\description{
  Returns an object from an object of class \code{fbst}
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
res = fbst(posteriorDensityDraws = p, nullHypothesisValue = 0, 
dimensionTheta = 3, dimensionNullset = 2)

# Return the e-value from an fbst object
res$eValue
}
