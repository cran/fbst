\name{fbet-class}
\Rdversion{1.0}
\docType{class}
\alias{fbet-class}
%%\alias{testMethod,testClass-method}
\title{Class \code{"fbet-class"}}
\description{Class for modelling the results of a Full Bayesian Evidence Test}
\section{Objects from the Class}{Store the results of a FBET}
\section{Slots}{
\describe{\item{\code{data}:}{Object of class \code{"list"} holding the results of the Full Bayesian Evidence Test. \code{posteriorDensityDraws} holds the posterior MCMC parameter draws, \code{posteriorDensityDrawsSorted} stores the sorted posterior MCMC parameter draws, \code{postDensValues} stores the posterior density values, \code{indices} stores the indices for deciding which values pass the evidence-threshold \eqn{\nu}, \code{interval} stores the boundaries of the interval null hypothesis, \code{referenceFunction} stores the name of the reference function used, \code{nu} specifies the evidence-threshold used for computation of the Bayesian evidence interval, \code{evidenceInterval} holds the endpoints of the resulting Bayesian evidence interval, \code{eValueH0} holds the Bayesian evidence value in favour of the interval null hypothesis, \code{eValueH1} holds the Bayesian evidence value in favour of the alternative hypothesis (or equivalently, against the interval null hypothesis)}}}
\keyword{classes}