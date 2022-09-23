## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(fbst)

## ---- warning=FALSE-----------------------------------------------------------
library(rstanarm)
set.seed(42)
mu_true = 0.25
sigma_true = 1
stanData = data.frame(y=rnorm(50,mu_true,sigma_true))
fit1 <- stan_glm(y ~ NULL, data = stanData,
                  family = gaussian(link = "identity"),
                  prior_intercept = normal(0,10),
                  prior_aux = exponential(1),
                  seed = 12345)

## -----------------------------------------------------------------------------
print(fit1)

## -----------------------------------------------------------------------------
posteriorDrawsMatrix = as.matrix(fit1)
head(posteriorDrawsMatrix)

## -----------------------------------------------------------------------------
resFlat = fbst(posteriorDrawsMatrix, nullHypothesisValue=0, dimensionTheta=2, dimensionNullset=1, dim = 2, gridSize = 1000)
summary(resFlat)

## -----------------------------------------------------------------------------
resFlat$eValue

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "500", out.height = "350"----
plot(resFlat, type = "contour", parNames=c("mu","sigma"))

## ----fig.align='center', dpi=300, fig.width = 7, fig.height = 5, out.width = "600", out.height = "450"----
plot(resFlat, type="persp", parNames=c("mu","sigma"))

## -----------------------------------------------------------------------------
resFlat$eValue

## -----------------------------------------------------------------------------
resFlat$sev_H_0

