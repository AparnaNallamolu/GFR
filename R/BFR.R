#' GFR
#'
#' Bayesian Functional Regression
#'
#' BFR is an modificated version of BGLR which implements a Gibbs sampler for a Bayesian regression model, this new version allows to... More details about the package.
#'
#' @seealso \code{\link[BGLR]{BGLR}}
#'
#' @docType package
"_PACKAGE"


#' BFR
#'
#' Bayesian Functional Regression
#'
#' More details about this function
#'
#' @param data (\code{data.frame}) the data with the $n$ $Response, also needs $Line and $Env for Cross Validation defined on it (NAs allowed).
#' @param datasetID (\code{string}) The name of the column with the identifiers of each line.
#' @param Multivariate By default, when de dataset includes more than one Environment and more than one Trait (MTME) the solution is adjusted by the method traditional, also is possible adjust the MTME by "SVD" <doi: >
#' @param a (\code{numeric}, $n$) Only requiered for censored outcomes. It's a vector specifying lower bounds for censored observation. By default is null.
#' @param b (\code{numeric}, $n$) Only requiered for censored outcomes. It's a vector specifying upper bounds for censored observation. By default is null.
#' @param ETA (\code{list}) Two level list used to specify the regression function, also could be generate by ETAGenerate() function for easy-use.
#' @param nIter (\code{integer}) The number of iterations.
#' @param burnIn (\code{integer}) The number of burn-in.
#' @param thin (\code{integer}) The number of thinning.
#' @param saveAt (\code{character}) This may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.
#' @param S0 (\code{numeric}) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes.
#' @param df0 (\code{numeric}) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes.
#' @param R2 (\code{numeric},$0<R2<1$) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified.
#' @param weights (\code{numeric}, $n$) a vector of weights, may be NULL. The residual variance of each data-point is set to be proportional to the inverse of the squared-weight. Only used with Gaussian outcomes.
#' @param verbose (\code{logical}) By default is \code{TRUE} and shows a fitting model progress bar if Folds <=1 or cross validation progress bar if Folds > 2.
#' @param rmExistingFiles (\code{logical}) By default is \code{TRUE} and removes existing output files from previous runs.
#' @param groups (\code{factor}) A vector of the same lenght of \code{data$Response} that associates observations with groups, each group will have an associated variance component for the error term.
#' @param CrossValidation (\code{list}) Especified list to KFold Crossvalidation use list(Type = 'KFold', nFolds = 5), and to Random Partiton Cross validation use list(Type = 'RandomPartition', nPartitions = 5, pTesting = 0.20, Traits.testing = NULL)
#' @param set_seed (\code{integer}) A seed for replicable research.
#' @param dec (\code{integer}) Number of decimals to show on the predictions.
#'
#' @importFrom stats lm rnorm var vcov
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom tidyr gather
#' @importFrom foreach %dopar%
#'
#' @export

BFR <- function(data = NULL, datasetID = 'Line', a=NULL, b=NULL, ETA = NULL, nIter = 1500, burnIn = 500, thin = 5,
                saveAt = '', S0 = NULL, df0 = 5, R2 = 0.5, weights = NULL, progressBar = TRUE, verbose = TRUE,
                rmExistingFiles = TRUE, groups = NULL, testingSet = NULL, set_seed = NULL, digits = 4, parallelCores = 1,
                Multivariate = NULL, CrossValidation = NULL, dec = NULL, response_type = NULL) {

  ########   DEPRECATED PARAMETERS  #########
  if(!is.null(CrossValidation)){Message('CrossValidation is deprecated, use testingSet instead, see more using help(BFR).', verbose)}
  if(!is.null(dec)){Message('dec is deprecated, use digits instead, see more using help(BFR).', verbose)}
  if(!is.null(Multivariate)){Message('Multivariate is deprecated, use BSVD() function instead, see more using help(BFR).', verbose)}
  if(!is.null(response_type)){Message('response_type is deprecated, use factor() or as.factor() in $Response column in data instead for ordinal analysis, or numeric to gaussian analysis, see more help(BFR).', verbose)}

  nCores <- detectCores()
  if(nCores < parallelCores){
    Message(paste0('The number of cores available is less than the specified in the parallelCores parameter we will use only ', nCores,'.'), verbose)
    parallelCores <- nCores
  }
  ########        INIT VALUES       #########
  response_type <- ifelse(is.factor(data$Response), 'ordinal', 'gaussian')
  data$Response <- as.numeric(data$Response)

  if (inherits(ETA, 'ETA')) {
    data <- ETA$data
    design <- ETA$Design
    datasetID <- ETA$ID
    ETA <- ETA$ETA
  } else {
    data <- validate.dataset(data, datasetID, orderData = F)
    design <- 'Handmade'
  }

  if (is.null(testingSet)) {
    out <- BGLR(data$Response, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights, verbose, rmExistingFiles, groups)
  } else if (parallelCores <= 1 && inherits(testingSet, 'CrossValidation')) {
    Tab_Pred <- data.frame()
    nCV <- length(testingSet$CrossValidation_list)
    pb <- progress::progress_bar$new(format = 'Fitting Cross-Validation :what  [:bar] Time elapsed: :elapsed', total = nCV, clear = FALSE, show_after = 0)

    for (actual_CV in seq_len(nCV)) {
      if (progressBar) {
        pb$tick(tokens = list(what = paste0(actual_CV, ' out of ', nCV)))
      }

      positionTST <- testingSet$CrossValidation_list[[actual_CV]]
      response_NA <-  data$Response
      response_NA[positionTST] <- NA
      time.init <- proc.time()[3]
      fm <- BGLR(response_NA, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                 verbose = F, rmExistingFiles, groups)

      Tab_Pred <- rbind(Tab_Pred,
                       data.frame(Position = positionTST,
                                  Environment = data$Env[positionTST],
                                  Trait = data$Trait[positionTST],
                                  Partition = actual_CV,
                                  Observed = round(data$Response[positionTST], digits),
                                  Predicted = round(fm$predictions[positionTST], digits)))
    }
    cat('Done!\n')
    out <- list(
      predictions_Summary = Tab_Pred,
      CrossValidation_list = testingSet$CrossValidation_list,
      response = data$Response,
      Design = design,
      nCores = parallelCores,
      response_type =  response_type
    )

    class(out) <- 'BFRCV'
  } else if (parallelCores > 1 && inherits(testingSet, 'CrossValidation')) {
    cl <- snow::makeCluster(parallelCores)
    doSNOW::registerDoSNOW(cl)
    nCV <- length(testingSet$CrossValidation_list)

    pb <- utils::txtProgressBar(max = nCV, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    Tab_Pred <- foreach::foreach(actual_CV = seq_len(nCV), .combine = rbind, .packages = 'GFR', .options.snow = opts) %dopar% {
      positionTST <- testingSet$CrossValidation_list[[actual_CV]]
      response_NA <-  data$Response
      response_NA[positionTST] <- NA
      time.init <- proc.time()[3]
      fm <- BGLR(response_NA, response_type, a, b, ETA, nIter, burnIn, thin, saveAt = paste0('CV', actual_CV,'_'), S0, df0, R2, weights,
                 verbose = F, rmExistingFiles, groups)

      data.frame(Position = positionTST,
                 Environment = data$Env[positionTST],
                 Trait = data$Trait[positionTST],
                 Partition = actual_CV,
                 Observed = round(data$Response[positionTST], digits),
                 Predicted = round(fm$predictions[positionTST], digits))
    }
    cat('Done!\n')
    out <- list(
      predictions_Summary = Tab_Pred,
      CrossValidation_list = testingSet$CrossValidation_list,
      response = data$Response,
      Design = design,
      nCores = parallelCores,
      response_type =  response_type
    )

    class(out) <- 'BFRCV'
  } else {
    positionTST <- testingSet
    response_NA <-  data$Response
    response_NA[positionTST] <- NA
    fm <- BGLR(response_NA, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
               verbose = T, rmExistingFiles, groups)
    Tab_Pred <- rbind(Tab_Pred,
                     data.frame(Position = positionTST,
                                Environment = data$Env[positionTST],
                                Trait = data$Trait[positionTST],
                                Partition = actual_CV,
                                Observed = round(data$Response[positionTST], digits),
                                Predicted = round(fm$predictions[positionTST], digits)))

    out <- list(
      predictions_Summary = Tab_Pred,
      CrossValidation_list = testingSet$CrossValidation_list,
      response = data$Response,
      Design = design,
      nCores = parallelCores,
      response_type =  response_type
    )

    class(out) <- 'BFRCV'
  }
  return(out)
}
