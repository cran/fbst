## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(fbst)

## -----------------------------------------------------------------------------
sleep

## -----------------------------------------------------------------------------
library(BayesFactor)
grp1 = sleep[1:10,]$extra
grp2 = sleep[11:20,]$extra
ttestBF(x=grp1,y=grp2, rscale="medium")

## -----------------------------------------------------------------------------
set.seed(42)
posteriorDraws = ttestBF(x=grp1,y=grp2, rscale="medium", posterior = TRUE, 
                         iterations = 100000)[,4]

## -----------------------------------------------------------------------------
set.seed(42)
result = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.1,0.1), nu=0)
summary(result)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result)

## -----------------------------------------------------------------------------
set.seed(42)
result2 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.2,0.2), nu=0)
summary(result2)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result2)

## -----------------------------------------------------------------------------
set.seed(42)
result3 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.1,0.1), 
               nu=1, FUN = dcauchy, par=list(location = 0, scale = sqrt(2)/2))
summary(result3)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result3)

## -----------------------------------------------------------------------------
result4 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.2,0.2), 
               nu=1, FUN = dcauchy, par=list(location = 0, scale = sqrt(2)/2))
summary(result4)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result4)

## -----------------------------------------------------------------------------
set.seed(42)
result5 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.2,0.2), nu=2, 
               FUN = dcauchy, par=list(location = 0, scale = sqrt(2)/2))
summary(result5)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result5)

