#' ETAGenerate
#'
#' Function to generate a Linear Predictor from a dataset to GFR.
#'
#' @param dataset TidyFormat
#' @param datasetID column with the identifiers.
#' @param GenomicMatrix lalalalal
#' @param REML Logical, By default is NULL, If is TRUE, the priorType is ommited else If is FALSE, the priorType is ommited and.
#' @param priorType Prior to assign, by default is 'FIXED', could be 'BRR', 'BayesA', 'BayesB', 'BayesC' or 'BL'
#' @param Bands Bands
#' @param Wavelengths Wavelenths
#' @param method Model to apply in bands, by default 'Alternative' will be used, also could be 'Simple', 'Complex'
#' @param basisType Basis function, by default is Fourier.Basis also could be Bspline.Basis.
#' @param nBasis Number of basis by default only use 1 basis.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' @export
#'
#' @examples
ETAGenerate <- function(dataset, datasetID = 'Line', GenomicMatrix = NULL, REML= NULL, priorType = 'FIXED', Bands = NULL, Wavelengths = NULL,
                        method = 'Simple', basisType = 'Fourier.Basis', nBasis = 1, ...) {
  dataset <- validate.dataset(dataset, datasetID)
  Design <- checkDesign(dataset, Bands, REML)
  priorType <- validate.prior(priorType)
  REML <- validateParadigm(REML)

  if (!is.null(GenomicMatrix)) {
    GenomicMatrix <- validate.GenomicMatrix(GenomicMatrix, dataset[, datasetID])
    L <- t(chol(GenomicMatrix))
    XL <- model.matrix(~0+as.factor(dataset[, datasetID])) %*% L
  } else {
    XL <- model.matrix(~0+as.factor(dataset[, datasetID]))
  }

  switch(Design,
         'Bayes-Single' = {
           ETA <- NULL
         }, 'Bayes-SingleBands' = {
           ETA <- list(Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), model = priorType))
         }, 'Bayes-Env' = {
           Env_prior <- ifelse(length(unique(dataset$Env)) < 10, "FIXED", priorType)
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = Env_prior),
                       Line = list(X = XL, model = priorType),
                       LinexEnv = list(X = model.matrix(~0+XL:as.factor(dataset$Env)), model = priorType))
         }, 'Bayes-EnvBands' = {
           Env_prior <- ifelse(length(unique(dataset$Env)) < 10, "FIXED", priorType)
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = Env_prior),
                       Line = list(X = XL, model = priorType),
                       LinexEnv = list(X = model.matrix(~0+XL:as.factor(dataset$Env)), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxEnv = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Env, ...), model = priorType))
         }, 'Bayes-Trait' = {
           Trait_prior <- ifelse(length(unique(dataset$Trait)) < 10, "FIXED", priorType)
           ETA <- list(Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = Trait_prior),
                       Line = list(X = XL, model = priorType),
                       LinexTrait = list(X = model.matrix(~0+XL:as.factor(dataset$Trait)), model = priorType))
         }, 'Bayes-TraitBands' = {
           Trait_prior <- ifelse(length(unique(dataset$Trait)) < 10, "FIXED", priorType)
           ETA <- list(Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = Trait_prior),
                       Line = list(X = XL, model = priorType),
                       LinexTrait = list(X = model.matrix(~0+XL:as.factor(dataset$Trait)), model = priorType),
                       Bands = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis), model = priorType, ...))
                       # BandsxTrait = list(X = bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, interaction = dataset$Trait, ...), model = priorType))
         }, 'Bayes-Multi' = {
           Env_prior <- ifelse(length(unique(dataset$Env)) < 10, "FIXED", priorType)
           Trait_prior <- ifelse(length(unique(dataset$Trait)) < 10, "FIXED", priorType)
           ETA <- list(Env = list(X = model.matrix(~0+as.factor(dataset$Env)), model = Env_prior),
                       Trait = list(X = model.matrix(~0+as.factor(dataset$Trait)), model = Trait_prior),
                       Line = list(X = XL, model = priorType),
                       LinexEnv = list(X = model.matrix(~0+XL:as.factor(dataset$Env)), model = priorType),
                       LinexTrait = list(X = model.matrix(~0+XL:as.factor(dataset$Trait)), model = priorType),
                       EnvxTrait = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait)), model = priorType),
                       EnvxTraitxLine = list(X = model.matrix(~0+as.factor(dataset$Env):as.factor(dataset$Trait):XL), model = priorType))
         }, 'Bayes-MultiBands' = {
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
         }, 'Frequentist-Single' = {
           ETA <- NULL
         }, 'Frequentist-SingleBands' = {
           ZBands <- model.matrix(~0+rownames(Bands))
           colnames(ZBands) <- rownames(Bands)
           Bands <- as.matrix(Bands)
           GBands <- Bands %*% t(Bands) / ncol(Bands)

           ETA <- list(Bands = list(Z = ZBands, K = GBands))
         }, 'Frequentist-Env' = {
           X <- model.matrix(~0+as.factor(dataset$Env))

           ZLines <- model.matrix(~0+as.factor(dataset$Line))
           ZLinEnv <- model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Env))

           if (!is.null(GenomicMatrix)) {
             GLinEnv <- kronecker(X = XL, Y = diag(length(unique(dataset$Env))))

             ETA <- list(X, Line = list(Z = ZLines, K = XL),
                         LinexEnv = list(Z = ZLinEnv, K = GLinEnv))
           } else {
            ETA <- list(X, Line = list(Z = ZLines),
                         LinexEnv = list(Z = ZLinEnv))
           }
         }, 'Frequentist-EnvBands' = {
           X <- model.matrix(~0+as.factor(dataset$Env))

           ZLines <- model.matrix(~0+as.factor(dataset$Line))
           ZLinEnv <- model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Env))
           ZBands <- model.matrix(~0+rownames(Wheat_Bands))
           colnames(ZBands) <- rownames(Bands)
           Bands <- as.matrix(Bands)
           GBands <- Bands %*% t(Bands) / ncol(Bands)

           if (!is.null(GenomicMatrix)) {
             GLinEnv <- kronecker(X = XL, Y = diag(length(unique(dataset$Env))))
             ETA <- list(X, Line = list(Z = ZLines, K = XL),
                         LinexEnv = list(Z = ZLinEnv, K = GLinEnv),
                         Bands = list(Z = ZBands, K = GBands))
           } else {
             ETA <- list(X, Line = list(Z = ZLines),
                         LinexEnv = list(Z = ZLinEnv),
                         Bands = list(Z = ZBands, K = GBands))
           }


         }, 'Frequentist-Trait' = {
           X <- model.matrix(~0+as.factor(dataset$Trait))

           ZLines <- model.matrix(~0+as.factor(dataset$Line))
           ZLinTrait <- model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Trait))

           if (!is.null(GenomicMatrix)) {
             GLinTrait <- kronecker(X = XL, Y = diag(length(unique(dataset$Trait))))
             ETA <- list(Line = list(Z = ZLines, K = XL),
                         LinexTrait = list(Z = ZLinTrait, K = GLinTrait))
           } else {
             ETA <- list(X, Line = list(Z = ZLines),
                         LinexTrait = list(Z = ZLinEnv))
           }

         }, 'Frequentist-TraitBands' = {
           X <- model.matrix(~0+as.factor(dataset$Trait))

           ZLines <- model.matrix(~0+as.factor(dataset$Line))
           ZLinTrait <- model.matrix(~0+as.factor(dataset$Line):as.factor(dataset$Trait))
           ZBands <- model.matrix(~0+rownames(Wheat_Bands))
           colnames(ZBands) <- rownames(Bands)
           Bands <- as.matrix(Bands)
           GBands <- Bands %*% t(Bands) / ncol(Bands)

           if (!is.null(GenomicMatrix)) {
             GLinTrait <- kronecker(X = XL, Y = diag(length(unique(dataset$Trait))))
             ETA <- list(X, Line = list(Z = ZLines, K = XL),
                         LinexTrait = list(Z = ZLinTrait, K = GLinTrait),
                         Bands = list(Z = ZBands, K = GBands))
           } else {
             ETA <- list(X, Line = list(Z = ZLines),
                         LinexTrait = list(Z = ZLinEnv),
                         Bands = list(Z = ZBands, K = GBands))
           }
         }, 'Frequentist-Multi' = {

         }, 'Frequentist-MultiBands' = {

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

validateParadigm <- function(REML) {
  changeInterpretation <- !is.null(REML)
  message("The statistical paradigm was changed from Bayes to Frequentist, so priorType parameter is ommited.")
}

checkDesign <- function(dataset, Bands = NULL, REML = NULL) {
  nEnv <- length(unique(dataset$Env))
  nTrait <- length(unique(dataset$Trait))
  bandsExist <- !is.null(Bands)
  changeInterpretation <- !is.null(REML)
  if (changeInterpretation) {
    return(checkDesign.Frequentist(nEnv, nTrait, bandsExist))
  } else {
    return(checkDesign.Bayes(nEnv, nTrait, bandsExist))
  }
}

checkDesign.Bayes <- function(nEnv, nTrait, bandsExist){
  if (nEnv == 1 & nTrait == 1) {
    if (bandsExist) {
      return('Bayes-SingleBands')
    } else {
      return('Bayes-Single')
    }
  } else if (nEnv > 1 & nTrait == 1) {
    if (bandsExist) {
      return('Bayes-EnvBands')
    } else {
      return('Bayes-Env')
    }
  } else if (nEnv == 1 & nTrait > 1) {
    if (bandsExist) {
      return('Bayes-TraitBands')
    } else {
      return('Bayes-Trait')
    }
  } else {
    if (bandsExist) {
      return('Bayes-MultiBands')
    } else {
      return('Bayes-Multi')
    }
  }
}

checkDesign.Frequentist <- function(nEnv, nTrait, bandsExist){
  if (nEnv == 1 & nTrait == 1) {
    if (bandsExist) {
      return('Frequentist-SingleBands')
    } else {
      return('Frequentist-Single')
    }
  } else if (nEnv > 1 & nTrait == 1) {
    if (bandsExist) {
      return('Frequentist-EnvBands')
    } else {
      return('Frequentist-Env')
    }
  } else if (nEnv == 1 & nTrait > 1) {
    if (bandsExist) {
      return('Frequentist-TraitBands')
    } else {
      return('Frequentist-Trait')
    }
  } else {
    if (bandsExist) {
      return('Frequentist-MultiBands')
    } else {
      return('Frequentist-Multi')
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
