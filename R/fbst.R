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
fbst <- function(posteriorDensityDraws, nullHypothesisValue=0, FUN=NULL, par=NULL, dimensionTheta, dimensionNullset, dim=1, gridSize = 1000){
  if(dim == 2){
    #  Check whether a matrix of posterior draws has been submitted
    if(class(posteriorDensityDraws)[1]!="matrix"){
      stop("Please provide a matrix with two columns for the argument posteriorDensityDraws.")
    }

    # STEP 1 - Multivariate kernel density estimation
    H = diag(c(1, 1)) # bandwith matrix
    kde <- ks::kde(x = posteriorDensityDraws, H = H, gridsize = c(gridSize, gridSize)) # grid of size 2000x2000 in R^2 for creating the kernel density

    # STEP 2 - Determine supremum of multivariate kernel density estimate over the null hypothesis set
    index_nullHypothesisValue = min(which(round(kde$eval.points[[1]],2) %in% c(nullHypothesisValue)))
    index_par2MaxTangentialSet = which.max(kde$estimate[index_nullHypothesisValue,])

    kde_nullHypothesisValue_par2MaxTangentialSet = kde$estimate[index_nullHypothesisValue,index_par2MaxTangentialSet]

    # STEP 3 - Evaluate the multivariate kernel density estimate at the posterior MCMC draws
    kde_eval <- ks::kde(x = posteriorDensityDraws, H = H, gridsize = c(gridSize, gridSize),
                        eval.points = posteriorDensityDraws) # grid of size gridSizexgridSize in R^2 for creating the kernel density
    # note that we specify eval.points as the posteriorDensityDraws, this allows later to use the kernel density estimates
    # at the values of the posterior draws we obtained
    # thus, the call below simply evaluates the kde above at the posteriorDensityDraws values


    # STEP 4 - Compute the evidence against the null hypothesis via Monte Carlo integration
    # iterate over posterior draws
    barEv = 0
    for(m in 1:nrow(kde_eval$eval.points)){
      if(kde_eval$estimate[m] > kde_nullHypothesisValue_par2MaxTangentialSet){
        # barEv = barEv + kde_eval$estimate[m]
        barEv = barEv + 1
      }
    }
    barEv = barEv / nrow(kde_eval$eval.points)

    # set reference function string
    refString = "Flat"

    # Compute standardized e-value sev(H_0)
    sev_H_0 = 1-pchisq(qchisq(barEv, df=dimensionTheta),df=dimensionTheta-dimensionNullset)
    
    # return fbst object
    res = new("fbst", data = list(posteriorDensityDraws=posteriorDensityDraws,
                                  postEffSizeSorted=sort(posteriorDensityDraws),
                                  densZero=kde_nullHypothesisValue_par2MaxTangentialSet,
                                  postDensValues=kde_eval$estimate,
                                  indices=c(index_nullHypothesisValue,index_par2MaxTangentialSet),
                                  nullHypothesisValue=nullHypothesisValue,
                                  referenceFunction=refString,
                                  dimensionTheta=dimensionTheta,
                                  dimensionNullset=dimensionNullset,
                                  eValue=barEv,
                                  pValue=NULL,
                                  sev_H_0 = sev_H_0))
  } else {
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
    # WARNING: The solution given in the 2008 paper is wrong. ev can be approximated as the upper tail of the 
    # cumulative chi_k^2 distribution function starting from ||m-M||^2
    p_value_ev_H_0 = 1-pchisq(d_0, df=dimensionNullset, lower.tail = TRUE)
    
    # Get indices of posterior draws belonging to the tangential set (flat prior)
    indices = which(postDensValues > densZero)
    
    # Calculate posterior probability mass of tangential set to H_0: delta = 0 (this is the Bayesian evidence value against H, \bar{ev})
    if (length(indices)>0){ # is there at least one posterior density value in the tangential set?
      barEv = cubature::cubintegrate(postDens, lower = postEffSizeSorted[min(indices)], upper = postEffSizeSorted[max(indices)])$integral
    } else { # tangential set is empty, Lebesgue-integral over it is zero then
      barEv = 0
    }
    
    if(barEv > 1){ # if integrate function makes some rounding errors set e-value against H_0 to 1, to avoid the quantile functions below to produce errors like NaN
      barEv =1
    }
    if(barEv < 0){
      barEv = 0
    }
    
    # Compute standardized e-value sev(H_0)
    sev_H_0 = 1-pchisq(qchisq(barEv, df=dimensionTheta),df=dimensionTheta-dimensionNullset)
    
    if (is.null(FUN)){
      refString = "Flat"
    } else {
      refString = "User-defined"
    }
    # return fbst object
    res = new("fbst", data = list(posteriorDensityDraws=posteriorDensityDraws,
                                  postEffSizeSorted=postEffSizeSorted,
                                  densZero=densZero,
                                  postDensValues=postDensValues,
                                  indices=indices,
                                  nullHypothesisValue=nullHypothesisValue,
                                  referenceFunction=refString,
                                  dimensionTheta=dimensionTheta,
                                  dimensionNullset=dimensionNullset,
                                  eValue=barEv,
                                  pValue=p_value_ev_H_0,
                                  sev_H_0 = sev_H_0))
  }
  return(res)
}


#' fbst class
#'
#' Stores the results of a Full Bayesian Significance Test
#'
#' @slot data A named list for storing the user-accessible data of an fbst object
#' posteriorDensityDraws A numeric (vector) of posterior MCMC parameter draws.
#' postEffSizeSorted A numeric (vector) of sorted posterior MCMC parameter draws.
#' densZero A numeric storing the surprise function value at the sharp null hypothesis parameter value.
#' postDensValues A numeric (vector) of posterior density values.
#' indices A numeric (vector) storing indices for deciding which values are located inside the tangential set.
#' nullHypothesisValue A numeric storing the sharp null hypothesis parameter value.
#' referenceFunction A character holding the name of the reference function used.
#' dimensionTheta A numeric holding the dimension of the parameter space.
#' dimensionNullset A numeric holding the dimension of the null set.
#' eValue A numeric holding the Bayesian evidence against the sharp null hypothesis, the e-value.
#' pValue A numeric holding the p-value associated with the Bayesian e-value in favour of the sharp null hypothesis.
#' sev_H_0 A numeric holding the standardized e-value as a replacement of the frequentist p-value
#' @name fbst-class
#' @rdname fbst-class
#' @export
setClass("fbst", representation(data="list"),
         prototype = NULL,
         validity = function(object) return(TRUE)
)



#' Conduct the Full Bayesian Evidence Test and return an fbet class object.
#'
#' @method 
#' @export
#' @examples
fbet <- function(posteriorDensityDraws=NULL, interval, nu=1, FUN=NULL, par=NULL, posterior=NULL, par_posterior=NULL){
  if(is.null(posterior) && is.null(posteriorDensityDraws)){
    stop("You must either provide the posteriorDensityDraws argument or the posterior argument!") 
  }
  if((!is.null(posterior)) && (!is.null(posteriorDensityDraws))){
    stop("You must either provide the posteriorDensityDraws argument or the posterior argument, not both!") 
  }
  
  if(is.null(posterior)){ # posterior draws provided as input
    # Sort posterior draws ascending
    postEffSizeSorted = sort(posteriorDensityDraws, decreasing = FALSE)
    # Construct density of posterior effect size
    postDens <- approxfun(density(postEffSizeSorted), rule = 2) # rule = 2 means: use closest data extreme for NA values produced
    
    # Calculate surprise function values for all posterior effect size samples
    if(!is.null(FUN)){ # reference function provided
      for(i in 1:length(par)){
        par[[i]] = rep(par[[i]],length(posteriorDensityDraws))
      }
      par = c(list(x=postEffSizeSorted),par)
      
      postDensValues = postDens(postEffSizeSorted)/do.call(FUN, par) # user defined reference function
    } else { # no reference function provided, then use flat improper prior as reference function
      postDensValues = postDens(postEffSizeSorted)/1 # flat reference function
    }
    
    if(nu==0 && is.null(FUN)){ # if flat ref function and nu == 0 the posterior is the surprise function
      evidenceInterval = c(-Inf,Inf)
      bayesianEvidenceValueIntervalNullHyp = cubintegrate(postDens, lower = interval[1], upper = interval[2])$integral
      bayesianEvidenceValueAlternative = 1-bayesianEvidenceValueIntervalNullHyp
      indices = "All"
      # rounding errors
      if(bayesianEvidenceValueIntervalNullHyp > 1){
        bayesianEvidenceValueIntervalNullHyp = 1
      }
      if(bayesianEvidenceValueAlternative > 1){
        bayesianEvidenceValueAlternative = 1
      }
    } else {
      # Compute indices of posterior parameter draws which belong to the expanded tangential set \tilde{T}(\nu)
      indices = which(postDensValues >= nu)
      
      # Get smallest and largest values as the boundaries of the evidence interval
      if(length(indices)>0){
        evidenceInterval = c(postEffSizeSorted[min(indices)],postEffSizeSorted[max(indices)])
      } else { # expanded tangential set is empty
        evidenceInterval = "Empty set"
      }
      
      if(!inherits(evidenceInterval, "character")){ # evidence interval is not the empty set? then compute Bayesian evidence values
        # Calculate the Bayesian evidence value for the interval null hypothesis
        if (length(indices)>0 && max(evidenceInterval[1],interval[1])<min(evidenceInterval[2],interval[2])){ # is there at least one posterior density value in the expanded tangential set?
          # integrate over intersection between the evidence interval and rope, indices [1] and [2] denote the lower and upper bounds of the evidence interval and rope
          bayesianEvidenceValueIntervalNullHyp = cubintegrate(postDens, lower = max(evidenceInterval[1],interval[1]), upper = min(evidenceInterval[2],interval[2]))$integral
        } else { # expanded tangential set is empty, Lebesgue-integral over an empty set is zero then
          bayesianEvidenceValueIntervalNullHyp = 0
        }
        # Calculate the Bayesian evidence value for the alternative hypothesis
        bayesianEvidenceValueAlternative = 0
        if(evidenceInterval[1]<interval[1] && evidenceInterval[2]<interval[1]){
          # intersection of Bayesian evidence interval and alternative hypothesis (left part)
          from.z <- evidenceInterval[1]
          to.z <- evidenceInterval[2]
          
          bayesianEvidenceValueAlternative = bayesianEvidenceValueAlternative + cubintegrate(postDens, lower = from.z, upper = to.z)$integral
        }
        
        if(evidenceInterval[1]>interval[2] && evidenceInterval[2]>interval[2]){
          # intersection of Bayesian evidence interval and alternative hypothesis (right part)
          from.z <- evidenceInterval[1]
          to.z <- evidenceInterval[2]
          
          bayesianEvidenceValueAlternative = bayesianEvidenceValueAlternative + cubintegrate(postDens, lower = from.z, upper = to.z)$integral
        }
        
        if(evidenceInterval[1]<interval[1] && evidenceInterval[2]>interval[1]){
          # intersection of Bayesian evidence interval and alternative hypothesis (right part)
          from.z <- evidenceInterval[1]
          to.z <- interval[1]
          
          bayesianEvidenceValueAlternative = bayesianEvidenceValueAlternative + cubintegrate(postDens, lower = from.z, upper = to.z)$integral
        }
        
        if(evidenceInterval[1]<interval[2] && evidenceInterval[2]>interval[2]){
          # intersection of Bayesian evidence interval and alternative hypothesis (right part)
          from.z <- interval[2]
          to.z <- evidenceInterval[2]
          
          bayesianEvidenceValueAlternative = bayesianEvidenceValueAlternative + cubintegrate(postDens, lower = from.z, upper = to.z)$integral
        }
        
        # rounding errors
        if(bayesianEvidenceValueIntervalNullHyp > 1){
          bayesianEvidenceValueIntervalNullHyp = 1
        }
        if(bayesianEvidenceValueAlternative > 1){
          bayesianEvidenceValueAlternative = 1
        }
      } else { # if evidence interval is empty, set Bayesian evidence values accordingly
        bayesianEvidenceValueIntervalNullHyp = 0
        bayesianEvidenceValueAlternative = 0
      }
    }
    if (is.null(FUN)){
      refString = "Flat"
    } else {
      refString = "User-defined"
    }
    
    res = new("fbet", data=list(
      posteriorDensityDraws = posteriorDensityDraws,
      posteriorDensityDrawsSorted = postEffSizeSorted,
      postDensValues = postDensValues,
      indices = indices,
      interval = interval,
      referenceFunction = refString,
      nu = nu,
      evidenceInterval = evidenceInterval, 
      eValueH0 = bayesianEvidenceValueIntervalNullHyp,
      eValueH1 = bayesianEvidenceValueAlternative))
  }
  if(is.null(posteriorDensityDraws)){ # calculation of FBET based on analytic posterior
    integrand <- function(x){
      par_posterior = c(x=x,par_posterior)
      numerator = do.call(posterior,par_posterior)
      par = c(x=x,par)
      denominator = do.call(FUN,par)
      if(numerator/denominator >= nu){
        res = numerator
      }
      else {
        res = 0
      }
      res
    }
    
    # integrate over null Hypothesis Interval
    EvH0 = cubintegrate(integrand,lower=interval[1],upper=interval[2])
    #cat("Evidence for the hypothesis specified inside interval argument:\n")
    res = EvH0$integral
  }
  res
}




#' fbet class
#'
#' Stores the results of a Full Bayesian Evidence Test
#'
#' @slot data A named list for storing the user-accessible data of an fbet object
#' posteriorDensityDraws A numeric (vector) of posterior MCMC parameter draws.
#' posteriorDensityDrawsSorted A numeric (vector) of sorted posterior MCMC parameter draws.
#' postDensValues A numeric (vector) of posterior density values.
#' indices A numeric (vector) storing indices for deciding which values are located inside the tangential set.
#' interval A numeric (vector) storing the boundaries of the interval null hypothesis (ROPE).
#' referenceFunction A character holding the name of the reference function used.
#' nu Evidence-threshold used for computation of the evidence interval
#' evidenceInterval A numeric (vector) storing the lower and upper boundaries of the resulting Bayesian evidence interval
#' eValueH0 Bayesian evidence value for the interval null hypothesis
#' eValueH1 Bayesian evidence value for the alternative hypothesis
#' @name fbet-class
#' @rdname fbet-class
#' @export
setClass("fbet", representation(data="list"),
         prototype = NULL,
         validity = function(object) return(TRUE)
)




#' plot object of class fbst
#' @usage \\method{plot}{fbst}(x, ...)
#' @export
plot.fbst <- function(x, ..., leftBoundary= -100, rightBoundary = 100, type = "contour", parNames = NULL, xlimleft = NULL, xlimright = NULL, xlabString = "Parameter", ylabString = NULL){
  if(inherits(x$posteriorDensityDraws, "numeric")){
    postDens <- approxfun(x=x$postEffSizeSorted,y=x$postDensValues, rule = 2)
    # prior-posterior plot
    if(is.null(xlimleft)){
      xlimleft = min(x$postEffSizeSorted)
    }
    if(is.null(xlimright)){
      xlimright = max(x$postEffSizeSorted)
    }
    if(is.null(ylabString)){
      ylabString = "density"
    }
    plot(x$postEffSizeSorted,x$postDensValues,ty="l",lty=1,xlim=c(xlimleft,xlimright),
         main="",ylab=ylabString,xlab=xlabString)
    
    # tangential area
    from.z <- x$postEffSizeSorted[min(x$indices)]
    to.z <- x$postEffSizeSorted[max(x$indices)]
    
    S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
    S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
    polygon(S.x,S.y, col = rgb(red = 0, green = 0.8, blue = 1, alpha = 1))
    
    # left tail (null) area
    from.z <- leftBoundary
    to.z <- x$postEffSizeSorted[min(x$indices)]
    
    # plot(postEffSizeSorted,postDensValues,ty="l", main = "", xlab = expression(paste("p(", delta, "| x)")), ylab = "Density")
    S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
    S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
    polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
    
    # right tail (null) area
    from.z <- x$postEffSizeSorted[max(x$indices)]
    to.z <- rightBoundary
    S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
    S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
    polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
    
    # null density value and separating line for tangential set
    points(x=0,y=x$densZero,col="blue",pch=19)
    abline(h=x$densZero,lty=2,lwd=1,col="blue")
  }
  if(inherits(x$posteriorDensityDraws, "matrix")){
    # Contour plot
    if(type == "contour"){
      H = diag(c(1, 1))
      kde <- ks::kde(x = x$posteriorDensityDraws, H = H, gridsize = c(1000, 1000))
      # contour plot
      graphics::image(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate,
            col = viridis::viridis(20), xlab = parNames[1], ylab = parNames[2])
      points(kde$x)
      index_nullHypothesisValue = min(which(round(kde$eval.points[[1]],2) %in% c(x$nullHypothesisValue)))
      index_par2MaxTangentialSet = which.max(kde$estimate[index_nullHypothesisValue,])
      kde_nullHypothesisValue_par2MaxTangentialSet = kde$estimate[index_nullHypothesisValue,index_par2MaxTangentialSet]
      points(kde$eval.points[[1]][index_nullHypothesisValue],kde$eval.points[[2]][index_par2MaxTangentialSet],col="magenta",pch=18)
    }
    if(type == "persp"){
      H = diag(c(1, 1))
      kde <- ks::kde(x = x$posteriorDensityDraws, H = H, gridsize = c(250, 250))
      plot(kde, display = "persp", col.fun = viridis::viridis, xlab = parNames[1], ylab = parNames[2])
    }
  }
}



#' plot object of class fbet
#' @usage \\method{plot}{fbet}(x, ...)
#' @export
plot.fbet <- function(x, ..., leftBoundary= -100, rightBoundary = 100, type = "posterior", legendposition = "topleft", main = ""){
  if(type=="surprise"){
    postDens <- approxfun(x=x$posteriorDensityDrawsSorted,y=x$postDensValues, rule = 2)
    # prior-posterior plot
    plot(x$posteriorDensityDrawsSorted,x$postDensValues,ty="l",lty=1,xlim=c(min(x$posteriorDensityDrawsSorted),max(x$posteriorDensityDrawsSorted)),
         main=main,ylab="density",xlab="Parameter")
  }
  if(type=="posterior"){
    postDens <- approxfun(density(x$posteriorDensityDraws))
    plot(density(x$posteriorDensityDraws),ty="l",lty=1,xlim=c(min(x$posteriorDensityDrawsSorted),max(x$posteriorDensityDrawsSorted)),
         main=main,ylab="density",xlab="Parameter")
  }
  
  # Visualize interval hypothesis
  abline(v=x$interval[1],lty=2,lwd=1,col="blue")
  abline(v=x$interval[2],lty=2,lwd=1,col="blue")
  
  # Visualize Bayesian evidence interval
  abline(v=x$evidenceInterval[1],lty=1,lwd=1,col="blue")
  abline(v=x$evidenceInterval[2],lty=1,lwd=1,col="blue")
  
  # intersection of Bayesian evidence interval and interval null hypothesis
  if (max(x$evidenceInterval[1],x$interval[1])<min(x$evidenceInterval[2],x$interval[2])){
    # intersection between the evidence interval and rope, indices [1] and [2] denote the lower and upper bounds of the evidence interval and rope
    from.z <- max(x$evidenceInterval[1],x$interval[1])
    to.z <- min(x$evidenceInterval[2],x$interval[2])
    
    S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
    S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
    polygon(S.x,S.y, col = rgb(red = 0, green = 0.8, blue = 1, alpha = 1))
  }

  # intersection of Bayesian evidence interval and alternative hypothesis
  if(x$evidenceInterval[1]<x$interval[1] && x$evidenceInterval[2]<x$interval[1]){
    # intersection of Bayesian evidence interval and alternative hypothesis (left part)
    from.z <- min(x$evidenceInterval[1],x$interval[1])
    to.z <- min(x$evidenceInterval[2],x$interval[1])
    
    if(from.z < to.z){
      S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
      S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
      polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
    }
  }
  
  if(is.infinite(x$evidenceInterval[1]) && is.infinite(x$evidenceInterval[2])){
    S.x  <- c(-1000, seq(-1000, x$interval[1], by = 0.001), x$interval[1])
    S.y  <- c(0, postDens(seq(-1000, x$interval[1], 0.001)), 0)
    polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
  }
  
  if(x$evidenceInterval[1]>x$interval[2] && x$evidenceInterval[2]>x$interval[2]){
  # intersection of Bayesian evidence interval and alternative hypothesis (right part)
    from.z <- x$evidenceInterval[1]
    to.z <- x$evidenceInterval[2]
    
    if(from.z < to.z){
      S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
      S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
      polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
    }
  }
  
  if(x$evidenceInterval[1]<x$interval[1] && x$evidenceInterval[2]>x$interval[1]){
    # intersection of Bayesian evidence interval and alternative hypothesis (right part)
    from.z <- x$evidenceInterval[1]
    to.z <- x$interval[1]
    if(is.infinite(from.z)){
      from.z = -100
    }
    
    if(from.z < to.z){
      S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
      S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
      polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
    }
  }
  
  if(x$evidenceInterval[1]<x$interval[2] && x$evidenceInterval[2]>x$interval[2]){
    # intersection of Bayesian evidence interval and alternative hypothesis (right part)
    from.z <- x$interval[2]
    to.z <- x$evidenceInterval[2]
    if(is.infinite(to.z)){
      to.z = 100
    }
    
    if(from.z < to.z){
      S.x  <- c(from.z, seq(from.z, to.z, by = 0.0001), to.z)
      S.y  <- c(0, postDens(seq(from.z, to.z, 0.0001)), 0)
      polygon(S.x,S.y, col = rgb(red = 1, green = 0, blue = 0, alpha = 1))
    }
  }
  
  if(type=="surprise"){
    # separating line for evidence-threshold nu used for the Bayesian evidence interval
    abline(h=x$nu,lty=2,lwd=1,col="black")
    
    graphics::legend(legendposition, legend=c("information function", "Null hypothesis", "Evidence Interval", "Evidence threshold"),
                     col=c("black", "blue", "blue", "black"), lty=c(1,2,1,2), cex=1)
  }
  if(type=="posterior"){
    graphics::legend(legendposition, legend=c("posterior density", "Null hypothesis", "Evidence Interval"),
                     col=c("black", "blue", "blue"), lty=c(1,2,1), cex=1)
  }
  
}




#' Print summary of an object of class fbst
#' @usage \\method{summary}{fbst}(object, ...)
#' @export
summary.fbst <- function(object, ...){
  cat("Full Bayesian Significance Test for testing a sharp hypothesis against its alternative:\n")
  cat("Reference function:", object$referenceFunction, "\n")
  cat("Hypothesis H_0:Parameter=", object$nullHypothesisValue, "against its alternative H_1\n")
  cat("Bayesian e-value against H_0:", object$eValue, "\n")
  
  cat("Standardized e-value:", object$sev_H_0, "\n")
}

#' Print summary of an object of class fbet
#' @usage \\method{summary}{fbet}(object, ...)
#' @export
summary.fbet <- function(object, ...){
  cat("Full Bayesian Evidence Test:\n")
  cat("Reference function:", object$referenceFunction, "\n")
  cat("Testing Hypothesis H_0:=[",object$interval[1], ",",object$interval[2],"]\n")
  cat("Bayesian e-value in favour of H_0:", object$eValueH0, "\n")
  cat("Bayesian e-value against H_0:", object$eValueH1, "\n")
}



#' Show method for an object of class fbst
#' @usage \\method{show}{fbst}(object)
#' @export
show.fbst <- function(object){
  cat("FBST for testing H_0:Parameter =", object$nullHypothesisValue, "against its alternative H_1\n")
  cat("Reference function:", object$referenceFunction, "\n")
  cat("Bayesian e-value against H_0:", object$eValue, "\n")
}

#' Show method for an object of class fbet
#' @usage \\method{show}{fbet}(object)
#' @export
show.fbet <- function(object){
  cat("FBET for testing H_0:Parameter in", object$interval, "against its alternative H_1\n")
  cat("Reference function:", object$referenceFunction, "\n")
  cat("Bayesian e-value in favour of H_0:", object$eValueH0, "\n")
}



#' Access results stored in the data slot of an object of class fbst
#' @usage \\method{$}{fbst}(object, ...)
#' @export
setMethod('$', signature="fbst",
          definition=function(x, name) {
            return_value = x@data[[name]]
            names(return_value) = name
            return(return_value)}
)

#' Access results stored in the data slot of an object of class fbet
#' @usage \\method{$}{fbet}(object, ...)
#' @export
setMethod('$', signature="fbet",
          definition=function(x, name) {
            return_value = x@data[[name]]
            names(return_value) = name
            return(return_value)}
)




#' Get names of the data slot of an object of class fbst
#' @usage \\method{names}{fbst}(object, ...)
#' @export
names.fbst <- function(x){
  names(x@data)
}

#' Get names of the data slot of an object of class fbet
#' @usage \\method{names}{fbet}(object, ...)
#' @export
names.fbet <- function(x){
  names(x@data)
}

#' Compute the Bayesian discrepancy measure (BDM)
#' @method 
#' @export
#' @examples
bdm <- function(posteriorDensityDraws, nullHypothesisValue=0){
  # Step 1: Compute median of posterior draws
  m=stats::median(posteriorDensityDraws)
  # Step 2: Define discrepancy interval
  if(m<nullHypothesisValue){
    I_H=c(m,nullHypothesisValue)
  }
  if(m==nullHypothesisValue){
    I_H=m
  }
  if(m>nullHypothesisValue){
    I_H=c(nullHypothesisValue,m)
  }
  # Step 3: Build kernel density estimate of posterior
  # Sort posterior draws ascending
  postEffSizeSorted = sort(posteriorDensityDraws, decreasing = FALSE)
  # Construct density of posterior effect size
  postDens <- approxfun(density(postEffSizeSorted), rule = 2) # rule = 2 means: use closest data extreme for NA values produced
  
  # Step 4: Integrate over discrepancy interval and double the value, which is the BDM
  bdm = 2 * cubature::cubintegrate(postDens, lower = I_H[1], upper = I_H[2])$integral
  
  # Step 5: return BDM
  return(bdm)
}

