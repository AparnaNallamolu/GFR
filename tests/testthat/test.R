library(testthat)
library(BGFRA)

context('Bspline Basis')

test_that('Bspline.Basis Test', {
  data('wheat_BGFRA')
  pos <- c(1:100)
  basis <- 21

  Bands_Drought <- Bands[pos,]
  Wavelengths <- Wavelengths
  Bspline <- Bspline.Basis(Bands = Bands_Drought, Wavelengths = Wavelengths, n.basis = basis, interaction = NULL)

  expect_is(Bspline, 'matrix')
  expect_equal(dim(Bspline)[1], dim(Bands_Drought)[1])
  expect_equal(dim(Bspline)[2], basis)
})

context('Fourier Basis')

test_that('Fourier.Basis Test', {
  data('wheat_BGFRA')
  pos <- c(1:100)
  basis <- 21

  Bands_Drought <- Bands[pos,]
  Wavelengths <- Wavelengths
  Fourier <- Fourier.Basis(Bands = Bands_Drought, Wavelengths = Wavelengths, n.basis = basis, interaction = NULL)
  expect_is(Fourier, 'matrix')

  expect_is(Fourier, 'matrix')
  expect_equal(dim(Fourier)[1], dim(Bands_Drought)[1])
  expect_equal(dim(Fourier)[2], basis)
})

context('CrossValidation Tests')
test_that('IBCF function', {
  data('wheat_BGFRA')

  Folds <- 5
  CV <- crossvalidation(Wheat, set_seed = 123)
  CV2 <- crossvalidation(Wheat, Folds, set_seed = 123)

  expect_equal(CV, CV2)
  expect_output(str(CV), 'List of 3')
  expect_false(any(is.na(CV$cv)))
  expect_equal(length(CV$cv), Folds)
  expect_equal(length(CV$ng), Folds)
  expect_is(CV$ng, 'numeric')
  expect_is(CV$n_CL, 'integer')
})

context('BGFRA Gaussian Tests')
test_that('BGFRA Fitting model Tests',{
  data('wheat_BGFRA')
  data <- Wheat

  ETA <- list(Env = list(X = model.matrix(~0+as.factor(data$Env)), model = "FIXED"),
               Line = list(X = model.matrix(~0+as.factor(data$Line)), model = "BRR"),
               Bands = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = NULL), model = "BRR"),
               BandsxEnv = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = data$Env), model = "BRR")
  )

  fm <- BGFRA(data, ETA = ETA, nIter = 1000, burnIn = 750, verbose = F)
  expect_output(str(fm), 'List of 21')

  expect_false(any(is.na(fm$response)))
  expect_false(any(is.na(fm$predictions)))
  expect_length(fm$response, length(fm$predictions))
})


test_that('BGFRA Predictive model Tests', {
  data('wheat_BGFRA')
  data <- Wheat

  ETA <- list(Env = list(X = model.matrix(~0+as.factor(data$Env)), model = "FIXED"),
              Line = list(X = model.matrix(~0+as.factor(data$Line)), model = "BRR"),
              Bands = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = NULL), model = "BRR"),
              BandsxEnv = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = data$Env), model = "BRR")
  )

  pm <- BGFRA(data, ETA = ETA, folds = 5, nIter = 1000, burnIn = 750, set_seed = 10, verbose = F)
  ## Error in data$Predictions[Pos_NA] <- predicted[Pos_NA] :
  ## replacement has length zero
  expect_output(str(pm), 'List of 4')
  expect_is(pm, 'BGFRACV')

  expect_output(str(pm$predictions_Summary), '18 obs. of  4 variables')
  expect_false(any(is.na(pm$predictions_Summary)))
  expect_is(pm$predictions_Summary[1, 1], 'character')
  #expect_is(pm$predictions_Summary[1, 2], 'numeric')

  expect_output(str(pm$cv), 'List of 5')
  expect_false(any(is.na(pm$cv)))

  expect_false(any(is.na(pm$response)))
  expect_false(any(is.na(pm$predictions)))
  expect_length(pm$response, length(pm$predictions))
})

context('BGFRA Binary Tests')
test_that('BGFRA Fitting model Tests', {
  data('wheat_BGFRA')
  data <- Wheat

  data$Response <- ifelse(data$Response > 3.5, 1, 0)

  ETA <- list(Env = list(X = model.matrix(~0+as.factor(data$Env)), model = "FIXED"),
              Line = list(X = model.matrix(~0+as.factor(data$Line)), model = "BRR"),
              Bands = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = NULL), model = "BRR"),
              BandsxEnv = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = data$Env), model = "BRR")
  )

  fm <- BGFRA(data, response_type = 'ordinal', ETA = ETA, nIter = 1000, burnIn = 750, verbose = F)
  expect_output(str(fm), 'List of 27')

  expect_false(any(is.na(fm$response)))
  expect_false(any(is.na(fm$predictions)))
  expect_length(fm$response, length(fm$predictions))
})


test_that('BGFRA Predictive model Tests',{
  data('wheat_BGFRA')
  data <- Wheat
  data$Response <- ifelse(data$Response > 3.5, 1, 0)

  ETA <- list(Env = list(X = model.matrix(~0+as.factor(data$Env)), model = "FIXED"),
              Line = list(X = model.matrix(~0+as.factor(data$Line)), model = "BRR"),
              Bands = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = NULL), model = "BRR"),
              BandsxEnv = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = data$Env), model = "BRR")
  )

  pm <- BGFRA(data, response_type = 'ordinal', ETA = ETA, folds = 5, nIter = 1000, burnIn = 750, verbose = F)

  expect_output(str(pm), 'List of 4')
  expect_is(pm, 'BGFRACV')

  expect_output(str(pm$predictions_Summary), '18 obs. of  4 variables')
  expect_false(any(is.na(pm$predictions_Summary)))
  expect_is(pm$predictions_Summary[1, 1], 'character')
  #expect_is(pm$predictions_Summary[1, 2], 'numeric')

  expect_output(str(pm$cv), 'List of 5')
  expect_false(any(is.na(pm$cv)))

  expect_false(any(is.na(pm$response)))
  expect_false(any(is.na(pm$predictions)))
  expect_length(pm$response, length(pm$predictions))
})

context('BGFRA Ordinal Tests')
test_that('BGFRA Fitting model Tests', {
  data('wheat_BGFRA')
  data <- Wheat
  y <- data$Response
  data$Response <-  ifelse(y < quantile(y,1/4),1, ifelse(y < quantile(y,3/4),2,3))

  ETA <- list(Env = list(X = model.matrix(~0+as.factor(data$Env)), model = "FIXED"),
              Line = list(X = model.matrix(~0+as.factor(data$Line)), model = "BRR"),
              Bands = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = NULL), model = "BRR"),
              BandsxEnv = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = data$Env), model = "BRR")
  )

  fm <- BGFRA(data, response_type = 'ordinal', ETA = ETA, nIter = 1000, burnIn = 750, verbose = F)
  expect_output(str(fm), 'List of 27')

  expect_false(any(is.na(fm$response)))
  expect_false(any(is.na(fm$predictions)))
  expect_length(fm$response, length(fm$predictions))
})


test_that('BGFRA Predictive model Tests',{
  data('wheat_BGFRA')
  data <- Wheat
  y <- data$Response
  data$Response <-  ifelse(y < quantile(y,1/4),1, ifelse(y < quantile(y,3/4),2,3))

  ETA <- list(Env = list(X = model.matrix(~0+as.factor(data$Env)), model = "FIXED"),
              Line = list(X = model.matrix(~0+as.factor(data$Line)), model = "BRR"),
              Bands = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = NULL), model = "BRR"),
              BandsxEnv = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = data$Env), model = "BRR")
  )

  pm <- BGFRA(data, response_type = 'ordinal', ETA = ETA, folds = 5, nIter = 1000, burnIn = 750, verbose = F)

  expect_output(str(pm), 'List of 4')
  expect_is(pm, 'BGFRACV')

  expect_output(str(pm$predictions_Summary), '18 obs. of  4 variables')
  expect_false(any(is.na(pm$predictions_Summary)))
  expect_is(pm$predictions_Summary[1, 1], 'character')
  #expect_is(pm$predictions_Summary[1, 2], 'numeric')

  expect_output(str(pm$cv), 'List of 5')
  expect_false(any(is.na(pm$cv)))

  expect_false(any(is.na(pm$response)))
  expect_false(any(is.na(pm$predictions)))
  expect_length(pm$response, length(pm$predictions))
})



path_delete <- dir()
path_delete <- path_delete[which(path_delete != 'test.R')]
file.remove(path_delete)




