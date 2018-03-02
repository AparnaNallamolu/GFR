#' ETAGenerate
#'
#' Function to generate a Linear Predictor from a dataset to GFR.
#'
#' @param dataset TidyFormat
#' @param datasetID column with the identifiers.
#' @param GenomicMatrix lalalalal
#' @param priorType Prior to assign, by default is 'FIXED', could be 'BRR', 'BayesA', 'BayesB', 'BayesC' or 'BL'
#' @param Bands Bands
#' @param Wavelengths Wavelenths
#' @param method Model to apply in bands, by default 'Alternative' will be used, also could be 'Simple', 'Complex'
#' @param basisType Basis function, by default is Fourier.Basis also could be Bspline.Basis.
#' @param nBasis Number of basis by default only use 1 basis.
#'
#'
#' @return
#' @export
#'
#' @examples
ETAGenerate <- function(dataset, datasetID = 'Line', GenomicMatrix = NULL, priorType = 'FIXED', Bands = NULL, Wavelengths = NULL,
                        method = 'Simple', basisType = 'Fourier.Basis', nBasis = 1, ...) {
  dataset <- validate.dataset(dataset, datasetID)
  Design <- checkDesign(dataset, Bands)
  priorType <- validate.prior(priorType)

  if (!is.null(GenomicMatrix)) {
    GenomicMatrix <- validate.GenomicMatrix(GenomicMatrix, dataset[, datasetID])
    L <- t(chol(GenomicMatrix))
    XL <- model.matrix(~0+as.factor(dataset[, datasetID])) %*% L
  } else {
    XL <- model.matrix(~0+as.factor(dataset[, datasetID]))
  }

  switch(Design,
         'Single' = {
           ETA <- NULL
         }, 'SingleBands' = {
           ETA <- list(Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), model = priorType))
         }, 'Env' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Line = list(X = XL, model = priorType),
                       LinexEnv = list(X = model.matrix(~0+XL:as.factor(dataset$Env)), model = priorType))
         }, 'EnvBands' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Line = list(X = XL, model = priorType),
                       LinexEnv = list(X = model.matrix(~0+XL:as.factor(dataset$Env)), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxEnv = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Env, ...), model = priorType))
         }, 'Trait' = {
           ETA <- list(Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = XL, model = priorType),
                       LinexTrait = list(X = model.matrix(~0+XL:as.factor(dataset$Trait)), model = priorType))
         }, 'TraitBands' = {
           ETA <- list(Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = XL, model = priorType),
                       LinexTrait = list(X = model.matrix(~0+XL:as.factor(dataset$Trait)), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxTrait = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Trait, ...), model = priorType))
         }, 'Multi' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = XL, model = priorType),
                       LinexEnv = list(X = model.matrix(~0+XL:as.factor(dataset$Env)), model = priorType),
                       LinexTrait = list(X = model.matrix(~0+XL:as.factor(dataset$Trait)), model = priorType),
                       EnvxTrait = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait)), model = priorType),
                       EnvxTraitxLine = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait):XL), model = priorType))
         }, 'MultiBands' = {
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = priorType),
                       Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = priorType),
                       Line = list(X = XL, model = priorType),
                       LinexEnv = list(X = model.matrix(~0+XL:as.factor(dataset$Env)), model = priorType),
                       LinexTrait = list(X = model.matrix(~0+XL:as.factor(dataset$Trait)), model = priorType),
                       EnvxTrait = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait)), model = priorType),
                       EnvxTraitxLine = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait):XL), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxEnv = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Env, ...), model = priorType),
                       # BandsxTrait = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Trait, ...), model = priorType))
         }, stop('Error in dataset or Bands provided, a bad design detected.')
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

validate.prior <- function(prior) {
  if (prior == 'RHKS') {prior <- 'BRR'}
  validPriors <- c('BRR', 'BayesA', 'BayesB', 'BayesC', 'BL')
  if (prior %in% validPriors) {
    return(prior)
  }

  stop('ERROR: The prior provided (', prior, ') is not available, check for misspelling or use a valid prior.')
}

validate.GenomicMatrix <- function(GenomicMatrix, Lines) {
  if (is.null(rownames(GenomicMatrix))) {
    stop("Row names of GenomicMatrix are not provided, use rownames() function to asign the names.")
  }

  GenomicDimension <- dim(GenomicMatrix)
  check1 <- GenomicDimension[1] == GenomicDimension[2]
  check2 <- GenomicDimension[1] == length(Lines)
  if (check1 && check2) {
    GenomicMatrix <- GenomicMatrix[Lines, ]
    return(GenomicMatrix)
  } else {
    stop('Error with the GenomicMatrix dimensions')
  }

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

checkDesign <- function(dataset, Bands = NULL) {
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
         }, 'Conventional' = {
           if (basisType == 'Bspline.Basis') {
             return(Bspline.Basis(Bands, Wavelengths, nBasis, ...))
           } else {
             return(Fourier.Basis(Bands, Wavelengths, nBasis, ...))
           }
         }, 'Alternative1' = {
           if (basisType == 'Bspline.Basis') {
             Phi <- fda::eval.basis(c(Wavelengths), fda::create.bspline.basis(range(c(Wavelengths)), nbasis = nBasis, breaks = NULL, norder = 4))
           } else {
             Phi <- fda::eval.basis(c(Wavelengths), fda::create.fourier.basis(range(c(Wavelengths)), nbasis = nBasis, period = period))
           }
           return(data.matrix(Bands) %*% Phi %*% solve(t(Phi) %*% Phi) %*% t(Phi))
         }, 'Alternative2' = {
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
