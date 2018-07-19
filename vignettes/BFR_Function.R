## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(GFR)

## ------------------------------------------------------------------------
data('Wheat_GFR')

## ------------------------------------------------------------------------
ETA1 <- ETAGenerate(Wheat_GFR, basisType = 'Fourier.Basis', Bands = Wheat_Bands, Wavelengths = Wheat_Wavelengths, priorType = 'BayesA', method = 'Alternative2', nBasis = 21)

## ------------------------------------------------------------------------
FM1 <- BFR(Wheat_GFR, ETA = ETA1, nIter = 1500, burnIn = 1000)

## ------------------------------------------------------------------------
CrossV <- list(Type = 'KFold', nFolds = 3)
PM1 <- BFR(ETA = ETA1, nIter = 1500, burnIn = 1000, CrossValidation = CrossV, set_seed = 10)
summary(PM1)
plot(PM1)
boxplot(PM1)

## ------------------------------------------------------------------------
ETA3 <- list(Env = list(X = model.matrix(~0+as.factor(Wheat_GFR$Env)), model = 'FIXED'),
             Line = list(X = model.matrix(~0+as.factor(Wheat_GFR$Line)), model = 'BRR'),
             Bands = list(X = Bspline.Basis(Wheat_Bands, Wheat_Wavelengths, nBasis = 23), model = 'BayesA'))

## ------------------------------------------------------------------------
FM2 <- BFR(Wheat_GFR, ETA = ETA3, nIter = 15000, burnIn = 10000, CrossValidation = CrossV, set_seed = 10)

## ------------------------------------------------------------------------
CrossV <- list(Type = 'KFold', nFolds = 5)

PM2 <- BFR(Wheat_GFR, ETA = ETA3, nIter = 15000, burnIn = 10000, CrossValidation = CrossV, set_seed = 10)
plot(PM2, select = 'MSEP')
boxplot(PM2, select = 'MSEP')

## ------------------------------------------------------------------------
cleanDat(T)

