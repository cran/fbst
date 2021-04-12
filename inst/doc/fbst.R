## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
sleep

## -----------------------------------------------------------------------------
library(BayesFactor)
grp1 = sleep[1:10,]$extra
grp2 = sleep[11:20,]$extra
ttestBF(x=grp1,y=grp2, rscale="medium")

## -----------------------------------------------------------------------------
posteriorDraws = ttestBF(x=grp1,y=grp2, rscale="medium", posterior = TRUE, iterations = 100000)[,4]
library(fbst)
result = fbst(posteriorDensityDraws = posteriorDraws, nullHypothesisValue = 0, dimensionTheta = 2, dimensionNullset = 1)
summary(result)

## ----fig.align='center',dpi=300,out.width="80%"-------------------------------
plot(result)

## -----------------------------------------------------------------------------
result2 = fbst(posteriorDensityDraws = posteriorDraws, nullHypothesisValue = 0, dimensionTheta = 2, dimensionNullset = 1, FUN = dcauchy, par=list(location = 0, scale = sqrt(2)/2))
summary(result2)

## ----fig.align='center',dpi=300,out.width="80%"-------------------------------
plot(result2)

## -----------------------------------------------------------------------------
result$sev_H_0
result2$sev_H_0

