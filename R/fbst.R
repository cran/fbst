# fbst - The Full Bayesian Significance Test and the e-value for testing 
# a sharp hypothesis against its alternative
#     Copyright (C) 2020  Riko Kelter
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Conduct the Full Bayesian Significance Test and return an fbst class object.
#'
#' @method 
#' @export
#' @examples
fbst <- function(posteriorDensityDraws, nullHypothesisValue=0, FUN=NULL, par=NULL, dimensionTheta, dimensionNullset){
  # Sort posterior draws ascending
  postEffSizeSorted = sort(posteriorDensityDraws, decreasing = FALSE)
  # Construct density of posterior effect size
  postDens <- approxfun(density(postEffSizeSorted), rule = 2) # rule = 2 means: use closest data extreme for NA values produced
    
  # Calculate posterior density at the null set, that is at null hypothesis H_0: delta = 0
  if(!is.null(FUN)) {
    parNull = c(list(x=nullHypothesisValue),par)
    densZero = postDens(nullHypothesisValue)/do.call(FUN, parNull) 
    }
  if(is.null(FUN)){ 
    densZero = postDens(nullHypothesisValue)/1 
  }
  
  # Calculate posterior density values for all posterior effect size samples
  if(!is.null(FUN)){
    par[[1]] = rep(par[[1]],length(posteriorDensityDraws))
    par[[2]] = rep(par[[2]],length(posteriorDensityDraws))
    par = c(list(x=postEffSizeSorted),par)

    postDensValues = postDens(postEffSizeSorted)/do.call(FUN, par) # user defined reference function
  } else {
    postDensValues = postDens(postEffSizeSorted)/1 # flat reference function
  }
  
  # Compute p-values for Bayesian evidence value in support of $H_0$
  m_0 = postDens(nullHypothesisValue)
  M_0 = bayestestR::map_estimate(postEffSizeSorted)[1]
  d_0 = abs((m_0-M_0)^2)
  p_value_ev_H_0 = pchisq(d_0, df=dimensionTheta, lower.tail = TRUE)
  
  # Get indices of posterior draws belonging to the tangential set (flat prior)
  indices = which(postDensValues > densZero)
  
  # Calculate posterior probability mass of tangential set to H_0: delta = 0 (this is the Bayesian evidence value against H, \bar{ev})
  if (length(indices)>0){ # is there at least one posterior density value in the tangential set?
    barEv = integrate(postDens, lower = postEffSizeSorted[min(indices)], upper = postEffSizeSorted[max(indices)])$value
  } else { # tangential set is empty, Lebesgue-integral over it is zero then
    barEv = 0
  }
  
  # Compute standardized e-value sev(H_0)
  sev_H_0 = 1-pchisq(qchisq(barEv, df=dimensionTheta),df=dimensionTheta-dimensionNullset)
  
  if (is.null(FUN)){
    priorString = "Flat"
  } else {
    priorString = "User-defined"
  }
  # return fbst object
  res = new("fbst",posteriorDensityDraws=posteriorDensityDraws,
                   postDensValues=postDensValues,
                   postEffSizeSorted=postEffSizeSorted,
                   densZero=densZero,
                   indices=indices,
                   nullHypothesisValue=nullHypothesisValue,
                   prior=priorString,
                   dimensionTheta=dimensionTheta,
                   dimensionNullset=dimensionNullset,
                   eValue=barEv,
                   pValue=p_value_ev_H_0,
                   sev_H_0 = sev_H_0)
  res
}



#' fbst class
#'
#' Stores the results of a Full Bayesian Significance Test
#'
#' @slot posteriorDensityDraws A numeric (vector) of posterior MCMC parameter draws.
#' @slot postEffSizeSorted A numeric (vector) of sorted posterior MCMC parameter draws.
#' @slot densZero A numeric storing the surprise function value at the sharp null hypothesis parameter value.
#' @slot postDensValues A numeric (vector) of posterior density values.
#' @slot indices A numeric (vector) storing indices for deciding which values are located inside the tangential set.
#' @slot nullHypothesisValue A numeric storing the sharp null hypothesis parameter value.
#' @slot prior A character holding the name of the reference function used.
#' @slot dimensionTheta A numeric holding the dimension of the parameter space.
#' @slot dimensionNullset A numeric holding the dimension of the null set.
#' @slot eValue A numeric holding the Bayesian evidence against the sharp null hypothesis, the e-value.
#' @slot pValue A numeric holding the p-value associated with the Bayesian e-value in favour of the sharp null hypothesis.
#' @slot sev_H_0 A numeric holding the standardized e-value as a replacement of the frequentist p-value
#' @name fbst-class
#' @rdname fbst-class
#' @export
setClass(Class="fbst",
         slots =c(
           posteriorDensityDraws="numeric",
           postEffSizeSorted="numeric",
           densZero="numeric",
           postDensValues="numeric",
           indices="numeric",
           nullHypothesisValue="numeric",
           prior="character",
           dimensionTheta="numeric",
           dimensionNullset="numeric",
           eValue="numeric",
           pValue="numeric",
           sev_H_0="numeric"
         ),
         prototype = NULL,         
         validity = function(object) return(TRUE)
)


#' plot object of class fbst
#' @usage \\method{plot}{fbst}(x, ...)
#' @export
plot.fbst <- function(x, ..., leftBoundary= -100, rightBoundary = 100){
  postDens <- approxfun(x=x@postEffSizeSorted,y=x@postDensValues, rule = 2)
  # prior-posterior plot
  plot(x@postEffSizeSorted,x@postDensValues,ty="l",lty=1,xlim=c(min(x@postEffSizeSorted),max(x@postEffSizeSorted)),
       main="",ylab="surprise function density",xlab="Parameter")
  
  # tangential area
  from.z <- x@postEffSizeSorted[min(x@indices)]
  to.z <- x@postEffSizeSorted[max(x@indices)]
  
  S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
  S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
  polygon(S.x,S.y, col = rgb(red = 0, green = 0.8, blue = 1, alpha = 1))
  
  # left tail (null) area
  from.z <- leftBoundary
  to.z <- x@postEffSizeSorted[min(x@indices)]
  
  # plot(postEffSizeSorted,postDensValues,ty="l", main = "", xlab = expression(paste("p(", delta, "| x)")), ylab = "Density")
  S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
  S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
  polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
  
  # right tail (null) area
  from.z <- x@postEffSizeSorted[max(x@indices)]
  to.z <- rightBoundary
  S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
  S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
  polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
  
  # null density value and separating line for tangential set
  points(x=0,y=x@densZero,col="blue",pch=19)
  abline(h=x@densZero,lty=2,lwd=1,col="blue")
}

#' Print summary of an object of class fbst
#' @usage \\method{summary}{fbst}(object, ...)
#' @export
summary.fbst <- function(object, ...){
  cat("Full Bayesian Significance Test for testing a sharp hypothesis against its alternative:\n")
  cat("Prior:", object@prior, "\n")
  cat("Testing Hypothesis H_0:Parameter=", object@nullHypothesisValue, "against its alternative H_1\n")
  cat("Bayesian e-value against H_0:", object@eValue, "\n")
  cat("p-value associated with the Bayesian e-value in favour of the null hypothesis:", object@pValue, "\n")
  cat("Standardized e-value:", object@sev_H_0, "\n")
}