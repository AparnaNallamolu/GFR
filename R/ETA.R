#' ETAGenerate
#'
#' Function to generate a Linear Predictor from a dataset to BFR.
#'
#' @param dataset TidyFormat
#' @param basisType Basis function, by default is Fourier.Basis also could be Bspline.Basis.
#' @param Bands Bands
#' @param Wavelengths Wavelenths
#' @param priorType Prior to assign, by default is 'FIXED', could be 'BRR', 'BayesA', 'BayesB', 'BayesC' or 'BL'
#' @param method Model to apply in bands, by default 'Alternative' will be used, also could be 'Simple', 'Complex'
#' @param nBasis Number of basis by default only use 1 basis.
#' @param datasetID column with the identifiers.
#'
#' @return
#' @export
#'
#' @examples
ETAGenerate <- function(dataset, datasetID = 'Line', priorType = 'FIXED', Bands = NULL, Wavelengths = NULL,
                        method = 'Alternative', basisType = 'Fourier.Basis', nBasis = 1, ...) {
  dataset <- validate.dataset(dataset, datasetID)
  Design <- check(dataset, Bands)
  switch(Design,
         'Single' = {
           ETA <- NULL
         }, 'SingleBands' = {
           ETA <- list(Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), model = priorType))
         }, 'Env' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Line = list(X = model.matrix(~0+as.factor(dataset$Line)), model = priorType),
                       LinexEnv = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Env)), model = priorType))
         }, 'EnvBands' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Line = list(X = model.matrix(~0+as.factor(dataset$Line)), model = priorType),
                       LinexEnv = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Env)), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxEnv = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Env, ...), model = priorType))
         }, 'Trait' = {
           ETA <- list(Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = model.matrix(~0+as.factor(dataset$Line)), model = priorType),
                       LinexTrait = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Trait)), model = priorType))
         }, 'TraitBands' = {
           ETA <- list(Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = model.matrix(~0+as.factor(dataset$Line)), model = priorType),
                       LinexTrait = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Trait)), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxTrait = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Trait, ...), model = priorType))
         }, 'Multi' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = model.matrix(~0+as.factor(dataset$Line)), model = priorType),
                       LinexEnv = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Env)), model = priorType),
                       LinexTrait = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Trait)), model = priorType),
                       EnvxTrait = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait)), model = priorType),
                       EnvxTraitxLine = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait):as.factor(dataset$Line)), model = priorType))
         }, 'MultiBands' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = model.matrix(~0+as.factor(dataset$Line)), model = priorType),
                       LinexEnv = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Env)), model = priorType),
                       LinexTrait = list(X = model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Trait)), model = priorType),
                       EnvxTrait = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait)), model = priorType),
                       EnvxTraitxLine = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait):as.factor(dataset$Line)), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxEnv = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Env, ...), model = priorType),
                       # BandsxTrait = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Trait, ...), model = priorType))
         }, stop('Error in dataset or Bands provided')
  )
  out <- list(ETA = ETA, #Linear predictor
              dataset = dataset, #Fixed Dataset
              ID = datasetID,
              Design = Design,
              Basis = basisType,
              Prior = priorType,
              Method = method
              ) #Type of model
  class(out) <- 'ETA'
  return(out)
}

validate.dataset <- function(dataset, datasetID, orderData=T) {
  if (is.null(dataset$Env)) {
    dataset$Env <- ''
  }

  if (is.null(dataset$Trait)) {
    dataset$Trait <- ''
  }

  if (is.null(dataset$Response)) {
    stop("No '$Response' provided in dataset")
  }

  try <- tryCatch({
    dataset[, datasetID]
  }, warning = function(w) {
    warning('Warning with the dataset')
  }, error = function(e) {
    stop("No identifier provided in dataset, use datesetID parameter to select the column of identifiers")
  })

  if (orderData) {
    dataset <- dataset[order(dataset$Trait, dataset$Env),]
  }
  dataset <- dataset[order(dataset$Trait, dataset$Env),]
  return(dataset)
}

check <- function(dataset, Bands = NULL) {
  nEnv <- length(unique(dataset$Env))
  nTrait <- length(unique(dataset$Trait))
  bandsExist <- !is.null(Bands)

  if (nEnv == 1 & nTrait == 1) {
    if (bandsExist) {
      return('SingleBands')
    } else {
      return('Single')
    }
  } else if (nEnv > 1 & nTrait == 1) {
    if (bandsExist) {
      return('EnvBands')
    } else {
      return('Env')
    }
  } else if (nEnv == 1 & nTrait > 1) {
    if (bandsExist) {
      return('TraitBands')
    } else {
      return('Trait')
    }
  } else {
    if (bandsExist) {
      return('MultiBands')
    } else {
      return('Multi')
    }
  }
}

bandsModel <- function(type, Bands, Wavelengths, basisType, nBasis, period = diff(range(c(Wavelengths))), ...) {
  switch(type,
         'Simple' = {
           return(data.matrix(Bands))
         }, 'Complex' = {
           if (basisType == 'Bspline.Basis') {
             Phi <- fda::eval.basis(c(Wavelengths), fda::create.bspline.basis(range(c(Wavelengths)), nbasis = nBasis, breaks = NULL, norder = 4))
           } else {
             Phi <- fda::eval.basis(c(Wavelengths), fda::create.fourier.basis(range(c(Wavelengths)), nbasis = nBasis, period = period))
           }
           return(data.matrix(Bands) %*% Phi %*% solve(t(Phi) %*% Phi) %*% t(Phi))
         }, 'Alternative' = {
           if (basisType == 'Bspline.Basis') {
             Phi <- fda::eval.basis(c(Wavelengths), fda::create.bspline.basis(range(c(Wavelengths)), nbasis = nBasis, breaks = NULL, norder = 4))
           } else {
             Phi <- fda::eval.basis(c(Wavelengths), fda::create.fourier.basis(range(c(Wavelengths)), nbasis = nBasis, period = period))
           }
           return(data.matrix(Bands) %*% Phi)
         },
         stop('Error: method ', type, ' no exist.')
  )
}
