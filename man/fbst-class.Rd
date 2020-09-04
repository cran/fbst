\name{fbst-class}
\Rdversion{1.0}
\docType{class}
\alias{fbst-class}
%%\alias{testMethod,testClass-method}
\title{Class \code{"fbst-class"}}
\description{Class for modelling the results of a Full Bayesian Significance Test}
\section{Objects from the Class}{Store the results of a FBST}
\section{Slots}{
\describe{\item{\code{data}:}{Object of class \code{"list"} holding the results of the Full Bayesian Significance Test. \code{posteriorDensityDraws} holds the posterior MCMC parameter draws, \code{postEffSizeSorted} stores the sorted posterior MCMC parameter draws, \code{densZero} stores the surprise function value at the sharp null hypothesis parameter value, \code{postDensValues} stores the posterior density values, \code{indices} stores the indices for deciding which values are located inside the tangential set, \code{nullHypothesisValue} stores the sharp null hypothesis parameter value, \code{referenceFunction} holds the name of the reference function used, \code{dimensionTheta} holds the dimension of the parameter space, \code{dimensionNullset} holds the dimension of the null set corresponding to the null hypothesis, \code{eValue} holds the Bayesian evidence against the sharp null hypothesis, the e-value, \code{pValue} holds the p-value associated with the Bayesian e-value in favour of the sharp null hypothesis, \code{sev_H_0} holds the standardized e-value as a replacement of the frequentist p-value.}}}
\keyword{classes}

