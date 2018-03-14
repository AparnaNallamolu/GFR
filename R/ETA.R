#' ETAGenerate
#'
#' Function to generate a Linear Predictor from a dataset to GFR.
#'
#' @param dataset TidyFormat
#' @param datasetID column with the identifiers.
#' @param Multivariate By default, when de dataset includes more than one Environment and more than one Trait (MTME) the solution is adjusted by the method traditional, also is possible adjust the MTME by "SVD" <doi: >
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
ETAGenerate <- function(dataset, datasetID = 'Line', Multivariate = 'Traditional', GenomicMatrix = NULL, REML= NULL, priorType = 'BRR', Bands = NULL, Wavelengths = NULL,
                        method = 'Simple', basisType = 'Fourier.Basis', nBasis = 1, ...) {
  dataset <- validate.dataset(dataset, datasetID, TRUE, Multivariate)
  Design <- checkDesign(dataset, Bands, REML)
  priorType <- validate.prior(priorType, Multivariate)
  REML <- validateParadigm(REML)
  # Bands <- validateBands(Bands, Wavelengths)
  # Wavelengths <- validateWavelengths()

  Env_prior <- ifelse(length(unique(dataset$Env)) < 10, "FIXED", priorType)
  Trait_prior <- ifelse(length(unique(dataset$Trait)) < 10, "FIXED", priorType)

  if (!is.null(GenomicMatrix)) {
    L <- t(chol(GenomicMatrix))
    XL <- model.matrix(~0+as.factor(dataset[, datasetID])) %*% L
    if (Multivariate == 'SVD') {
      XL <- kronecker(X = GenomicMatrix, Y = diag(length(unique(dataset$Env))))
    }
    rownames(XL) <- dataset[, datasetID]
    XL <- validate.GenomicMatrix(XL, dataset[, datasetID])
    Line_prior <- 'RKHS'
    LineK <- TRUE
  } else {
    XL <- model.matrix(~0+as.factor(dataset[, datasetID]))
    Line_prior <- priorType
    LineK <- FALSE
  }

  switch(Design,
         'Bayes-Single' = {
          ETA <- NULL
         }, 'Bayes-SingleBands' = {
          ETA <- list(Bands = ETAList(bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), priorType, TRUE))
         }, 'Bayes-Env' = {
          ETA <- list(Env = ETAList(dataset$Env, Env_prior),
                      Line = ETAList(XL, Line_prior, TRUE, LineK),
                      LinexEnv = ETAList(XL, priorType, interaction1 = dataset$Env))
         }, 'Bayes-EnvBands' = {
          ETA <- list(Env = ETAList(dataset$Env, Env_prior),
                      Line = ETAList(XL, Line_prior, TRUE, LineK),
                      LinexEnv = ETAList(XL, priorType, interaction1 = dataset$Env),
                      Bands = ETAList(bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), priorType, TRUE))
         }, 'Bayes-Trait' = {
          ETA <- list(Trait = ETAList(dataset$Trait, Trait_prior),
                      Line = ETAList(XL, Line_prior, TRUE, LineK),
                      LinexTrait = ETAList(XL, priorType, interaction1 = dataset$Trait))
         }, 'Bayes-TraitBands' = {
          ETA <- list(Trait = ETAList(dataset$Trait, Trait_prior),
                      Line = ETAList(XL, Line_prior, TRUE, LineK),
                      LinexTrait = ETAList(XL, priorType, interaction1 = dataset$Trait),
                      Bands = ETAList(bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), priorType, TRUE))
         }, 'Bayes-Multi' = {
          if (Multivariate == 'Traditional') {
            ETA <- list(Env = ETAList(dataset$Env, Env_prior),
                        Trait = ETAList(dataset$Trait, Trait_prior),
                        Line = ETAList(XL, Line_prior, TRUE, LineK),
                        LinexEnv = ETAList(XL, priorType, interaction1 = dataset$Env),
                        LinexTrait = ETAList(XL, priorType, interaction1 = dataset$Trait),
                        EnvxTrait = ETAList(dataset$Env, priorType, interaction1 = dataset$Trait),
                        EnvxTraitxLine = ETAList(dataset$Env, priorType, interaction1 = dataset$Trait, interaction2 = dataset$Line))
          } else if (Multivariate == 'SVD') {
            ETA <- list(Env = ETAList(dataset$Env, Env_prior),
                        # Trait = list(X = model.matrix(~0+as.factor(data_Long$Trait)), model = Trait_prior),
                        Line = ETAList(XL, Line_prior, TRUE, LineK),
                        LinexEnv = ETAList(XL, priorType, interaction1 = dataset$Env, likeKernel = TRUE, withK = T))
                        # LinexTrait = list(X = model.matrix(~0+XL:as.factor(data_Long$Trait)), model = priorType),
                        # EnvxTrait = list(X = model.matrix(~0+as.factor(data_Long$Env):as.factor(data_Long$Trait)), model = priorType),
                        # EnvxTraitxLine = list(X = model.matrix(~0+as.factor(data_Long$Env):as.factor(data_Long$Trait):XL), model = priorType))

          }

         }, 'Bayes-MultiBands' = {
           if (Multivariate == 'Traditional') {
             ETA <- list(Env = ETAList(dataset$Env, Env_prior),
                         Trait = ETAList(dataset$Trait, Trait_prior),
                         Line = ETAList(XL, priorType, TRUE, LineK),
                         LinexEnv = ETAList(XL, priorType, interaction1 = dataset$Env),
                         LinexTrait = ETAList(XL, priorType, interaction1 = dataset$Trait),
                         EnvxTrait = ETAList(dataset$Env, priorType, interaction1 = dataset$Trait),
                         EnvxTraitxLine = ETAList(dataset$Env, priorType, interaction1 = dataset$Trait, interaction2 = dataset$Line),
                         Bands = ETAList(bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), priorType, TRUE))
           } else if (Multivariate == 'SVD') {
             ETA <- list(Env = ETAList(dataset$Env, Env_prior),
                         # Trait = list(X = model.matrix(~0+as.factor(data_Long$Trait)), model = Trait_prior),
                         Line = ETAList(XL, priorType, TRUE, LineK),
                         LinexEnv = ETAList(XL, priorType, interaction1 = dataset$Env, likeKernel = TRUE, withK = T),
                         Bands = ETAList(bandsModel(method, Bands, Wavelengths, basisType, nBasis = nBasis, ...), priorType, TRUE, likeKernel = TRUE, withK = T))
                         # LinexTrait = list(X = model.matrix(~0+XL:as.factor(data_Long$Trait)), model = priorType),
                         # EnvxTrait = list(X = model.matrix(~0+as.factor(data_Long$Env):as.factor(data_Long$Trait)), model = priorType),
                         # EnvxTraitxLine = list(X = model.matrix(~0+as.factor(data_Long$Env):as.factor(data_Long$Trait):XL), model = priorType))

           }
         }, 'Frequentist-Single' = {
           ETA <- NULL
         }, 'Frequentist-SingleBands' = {
           ETA <- list(Bands = ZList(rownames(Bands), Bands %*% t(Bands) / ncol(Bands)))
         }, 'Frequentist-Env' = {
           X <- model.matrix(~0+as.factor(dataset$Env))
           ETA <- list(X, Line = ZList(dataset$Line, K = XL, WithK = !is.null(GenomicMatrix)),
                       LinexEnv = ZList(dataset$Line, K = kronecker(X = XL, Y = diag(length(unique(dataset$Env)))), interaction1 = dataset$Env, WithK = !is.null(GenomicMatrix)))
         }, 'Frequentist-EnvBands' = {
           X <- model.matrix(~0+as.factor(dataset$Env))
           ETA <- list(X, Line = ZList(Z = dataset$Line, K = XL, WithK = !is.null(GenomicMatrix)),
                       LinexEnv = ZList(Z = dataset$Line, K = kronecker(X = XL, Y = diag(length(unique(dataset$Env)))), interaction1 = dataset$Env, WithK = !is.null(GenomicMatrix) ),
                       Bands = ZList(Z = rownames(Bands), K = Bands %*% t(Bands) / ncol(Bands)))
         }, 'Frequentist-Trait' = {
           X <- model.matrix(~0+as.factor(dataset$Trait))
           ETA <- list(X, Line = ZList(Z = dataset$Line, K = XL, WithK = !is.null(GenomicMatrix)),
                       LinexTrait = ZList(Z = dataset$Line, K = kronecker(X = XL, Y = diag(length(unique(dataset$Trait)))), interaction1 = dataset$Trait, WithK = !is.null(GenomicMatrix) ),
                       Bands = ZList(Z = rownames(Bands), K = Bands %*% t(Bands) / ncol(Bands)))
         }, 'Frequentist-TraitBands' = {
           X <- model.matrix(~0+as.factor(dataset$Trait))
           ETA <- list(X, Line = ZList(Z = dataset$Line, K = XL, WithK = !is.null(GenomicMatrix)),
                       LinexTrait = ZList(Z = dataset$Line, K = kronecker(X = XL, Y = diag(length(unique(dataset$Trait)))), interaction1 = dataset$Trait, WithK = !is.null(GenomicMatrix) ),
                       Bands = ZList(Z = rownames(Bands), K = Bands %*% t(Bands) / ncol(Bands)))
         }, 'Frequentist-Multi' = {

         }, 'Frequentist-MultiBands' = {

         }, Error('Error in dataset or Bands provided, a bad design detected.')
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

validate.prior <- function(prior, Multivariate) {
  validPriors <- c('BRR', 'BayesA', 'BayesB', 'BayesC', 'BL', 'RHKS', 'FIXED')

  if (prior %in% validPriors) {
    if (Multivariate == 'SVD') {
      return('RKHS')
    }
    return(prior)
  }
  Error('The prior provided (', prior, ') is not available, check for misspelling or use a valid prior.')
}

validate.GenomicMatrix <- function(GenomicMatrix, Lines) {
  if (is.null(rownames(GenomicMatrix))) {
    Error("Row names of GenomicMatrix are not provided, use rownames() function to asign the names.")
  }

  GenomicDimension <- dim(GenomicMatrix)
  check1 <- GenomicDimension[1] == GenomicDimension[2]
  check2 <- GenomicDimension[1] == length(Lines)
  # if (!check1) {
  #   GenomicMatrix <- GenomicMatrix %*% t(GenomicMatrix)
  #   check1 <- GenomicDimension[1] == GenomicDimension[2]
  # }

  if (check1 && check2) {
    GenomicMatrix <- GenomicMatrix[Lines, ]
    return(GenomicMatrix)
  } else {
    Error('Error with the GenomicMatrix dimensions')
  }

}

validate.dataset <- function(dataset, datasetID, orderData = T, Multivariate = 'Traditional') {
  if (is.null(dataset$Env)) {
    dataset$Env <- ''
  } else {
    dataset$Env <- as.factor(dataset$Env)
  }

  if (is.null(dataset$Trait)) {
    dataset$Trait <- ''
  } else {
    dataset$Trait <- as.factor(dataset$Trait)
  }

  if (is.null(dataset$Response)) {
    Error("No '$Response' provided in dataset")
  }

  dataset[, datasetID] <- tryCatch({
     as.factor(dataset[, datasetID])
  }, warning = function(w) {
    Warning('Warning with the dataset')
  }, error = function(e) {
    Error("No identifier provided in dataset, use datesetID parameter to select the column of identifiers")
  })

  if (Multivariate == 'SVD') {
    dataset <- tidyr::spread(dataset, 'Trait', 'Response')
  }

  if (orderData) {
    if (length(unique(dataset$Trait)) > 1) {
      dataset <- dataset[order(dataset$Trait, dataset$Env),]
    }
    dataset <- dataset[order(dataset$Env),]
  }

  return(dataset)
}

validateParadigm <- function(REML) {
  changeInterpretation <- !is.null(REML)
  if (changeInterpretation) {
    Message("The statistical paradigm was changed from Bayes to Frequentist, so priorType parameter is ommited.")
  }
  return(REML)
}

validateBands <- function(Bands, Wavelengths){
  Bands <- as.matrix(Bands)
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
         Error('Error: method ', type, ' no exist.')
  )
}

ETAList <- function(variable, priorType = 'BRR', no.model = FALSE, withK = FALSE, likeKernel = FALSE, interaction1 = NULL, interaction2 = NULL) {
  if (no.model) {
    out <- list(X = variable, model = priorType) ## Witout model.matrix function
  } else if (!is.null(interaction2)) {
    out <- list(X = model.matrix(~0+variable:interaction1:interaction2), model = priorType) # With double interaction case
  } else if (!is.null(interaction1)) {
    out <- list(X = model.matrix(~0+variable:interaction1), model = priorType) # With one interaction case
  } else {
    out <- list(X = model.matrix(~0+variable), model = priorType) # Without interaction case
  }

  if (withK) {
    if (likeKernel) {
      out[['X']] <- out[['X']] %*% t(out[['X']])# / ncol(out[['X']])
    }
    names(out)[1] <- 'K' #Rename X by K to kernel use.
  }
  return(out)
}

ZList <- function(variable, K = NULL, no.model = FALSE, interaction1 = NULL, interaction2 = NULL, withK = FALSE){
  if (!withK) {
    K = NULL
  }

  if (no.model) {
    out <- list(Z = variable, K = K) ## Witout model.matrix function
  } else if (!is.null(interaction2)) {
    Z.m <- model.matrix(~0+variable:interation1:interation2)
    out <- list(Z = Z.m, K = K) # With double interaction case
  } else if (!is.null(interaction1)) {
    Z.m <- model.matrix(~0+variable:interation1)
    out <- list(Z = model.matrix(~0+variable:interaction1), K = K) # With one interaction case
  } else {
    Z.m <- model.matrix(~0+variable)
    colnames(Z.m) <- rownames(variable)
    out <- list(Z = model.matrix(~0+variable), K = K) # Without interaction case
  }
  return(out)
}
