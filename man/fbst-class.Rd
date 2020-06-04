\name{fbst-class}
\Rdversion{1.0}
\docType{class}
\alias{fbst-class}
%%\alias{testMethod,testClass-method}
\title{Class \code{"fbst-class"}}
\description{Class for modelling the results of a Full Bayesian Significance Test}
\section{Objects from the Class}{Store the results of a FBST}
\section{Slots}{
\describe{\item{\code{posteriorDensityDraws}:}{Object of class \code{"numeric"} holding the posterior MCMC parameter draws. ~~ }
          \item{\code{postEffSizeSorted}:}{Object of class \code{"numeric"} storing the sorted posterior MCMC parameter draws. ~~ }           \item{\code{densZero}:}{Object of class \code{"numeric"} storing the surprise function value at the sharp null hypothesis parameter value. ~~ }
          \item{\code{postDensValues}:}{Object of class \code{"numeric"} storing the posterior density values. ~~ }
          \item{\code{indices}:}{Object of class \code{"numeric"} storing indices for deciding which values are located inside the tangential set. ~~ }
          \item{\code{nullHypothesisValue}:}{Object of class \code{"numeric"} storing the sharp null hypothesis parameter value. ~~ }
          \item{\code{prior}:}{Object of class \code{"character"} holding the name of the reference function used. ~~ }
          \item{\code{dimensionTheta}:}{Object of class \code{"numeric"} holding the dimension of the parameter space. ~~ }
          \item{\code{dimensionNullset}:}{Object of class \code{"numeric"} holding the dimension of the null set corresponding to the null hypothesis. ~~ }
          \item{\code{eValue}:}{Object of class \code{"numeric"} holding the Bayesian evidence against the sharp null hypothesis, the e-value. ~~ }
          \item{\code{pValue}:}{Object of class \code{"numeric"} holding the p-value associated with the Bayesian e-value in favour of the sharp null hypothesis. ~~ }
          \item{\code{sev_H_0}:}{Object of class \code{"numeric"} holding the standardized e-value as a replacement of the frequentist p-value. ~~ }}}
\keyword{classes}

