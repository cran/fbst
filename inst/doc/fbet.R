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
require(BayesFactor)
grp1 = sleep[1:10,]$extra
grp2 = sleep[11:20,]$extra
diff = grp1 - grp2
ttestBF(x = grp1, y = grp2, rscale = "medium", paired = TRUE)

## -----------------------------------------------------------------------------
set.seed(42)
posteriorDraws = ttestBF(x = grp1, y = grp2, rscale = "medium", paired = TRUE,
                         posterior = TRUE, iterations = 100000)[,3]

## -----------------------------------------------------------------------------
set.seed(42)
result = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.2,0.2), nu=0)
summary(result)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result, type = "posterior")

## -----------------------------------------------------------------------------
set.seed(42)
result2 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.25,0.25), nu=0)
summary(result2)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result2, type = "posterior")

## -----------------------------------------------------------------------------
set.seed(42)
result3 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.5,0.5), 
               nu=1, FUN = dcauchy, par=list(location = 0, scale = sqrt(2)/2))
summary(result3)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result3, type="surprise")

## -----------------------------------------------------------------------------
result4 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.2,0.2), 
               nu=1, FUN = dcauchy, par=list(location = 0, scale = sqrt(2)/2))
summary(result4)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result4, type="surprise")

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result4)

## -----------------------------------------------------------------------------
set.seed(42)
result5 = fbet(posteriorDensityDraws = posteriorDraws, interval = c(-0.2,0.2), nu=2, 
               FUN = dcauchy, par=list(location = 0, scale = sqrt(2)/2))
summary(result5)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result5)

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(result5, type = "surprise")

