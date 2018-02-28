#' BFR
#'
#' Bayesian Functional Regression
#'
#' More details about the package.
#'
#' @docType package
"_PACKAGE"


#' BFR
#'
#' Bayesian Functional Regression
#'
#' BFR is an modificated version of BGLR which implements a Gibbs sampler for a Bayesian regression model, this new version allows to
#'
#'
#' @param data (\code{data.frame}) the data with the $n$ $Response, also needs $Line and $Env for Cross Validation defined on it (NAs allowed).
#' @param response_type (\code{character}) It can be 'gaussian' or 'ordinal'.
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
#'
#' @seealso \code{\link[BGLR]{BGLR}}
#'
#' examples
#'
#' @export

BFR <- function(data = NULL, response_type = 'gaussian', a=NULL, b=NULL, ETA = NULL, nIter = 1500,
                  burnIn = 500, thin = 5, saveAt = '', S0 = NULL, df0 = 5, R2 = 0.5, weights = NULL,
                  verbose = TRUE, rmExistingFiles = TRUE, groups = NULL, CrossValidation = NULL,
                  set_seed = NULL, dec = 4){
  if (inherits(ETA, 'ETA')) {
    data <- ETA$data
    ETA <- ETA$ETA
  } else {
    data <- validate.dataset(data)
  }

  if (!is.null(CrossValidation)) {

    switch(CrossValidation$Type,
           KFold = {
             PT <- CV.KFold(data, K = CrossValidation$nFolds, set_seed)
             nCV <- CrossValidation$nFolds
           },
           RandomPartition = {
             PT <- CV.RandomPart(data, NPartitions = CrossValidation$nPartitions, PTesting = CrossValidation$pTesting, Traits.testing = CrossValidation$Traits.testing, set_seed)
             nCV <- CrossValidation$nPartitions
           },
           stop(paste0('ERROR: The Cross Validation  ', CrossValidation$Type, " is't implemented"))
    )

    if (verbose) {
      cat("This might be time demanding, let's take sit and a cup of coffe\n")

      pb <- progress::progress_bar$new(format = 'Fitting the :what  [:bar] Time elapsed: :elapsed', total = nCV + 1, clear = FALSE, show_after = 0)
    }

    data$Predictions <- NA
    Tab_Pred <- data.frame()
    ## Init cross validation
    for (i in seq_len(nCV)) {

      if (verbose) {
        pb$tick(tokens = list(what = paste0( i, ' CV of ', nCV)))
      }

      response_NA <-  data$Response
      Pos_NA <- PT$cv[[paste0('partition',i)]]
      response_NA[Pos_NA] <- NA
      time.init <- proc.time()[3]
      fm <- BGLR(response_NA, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                       verbose = F, rmExistingFiles, groups)


      switch(response_type,
        gaussian = {
          predicted <- fm$predictions
          data$Predictions[Pos_NA] <- predicted[Pos_NA]
          Tab <- data.frame(Env = data$Env[Pos_NA], Fold = i, y_p = predicted[Pos_NA],
                            y_o = data$Response[Pos_NA])
          Tab_Pred <- rbind(Tab_Pred, Cor_Env(Tab, Time = proc.time()[3] - time.init ))
        },
        ordinal = {
          predicted <- as.integer(colnames(fm$probs)[apply(fm$probs,1,which.max)])
          data$Predictions[Pos_NA] <- predicted[Pos_NA]

          Tab <- data.frame(Env = data$Env[Pos_NA], Fold = i, y_p = predicted[Pos_NA],
                            y_o = data$Response[Pos_NA] )
          Tab_Pred <- rbind(Tab_Pred, Cor_Env_Ordinal(Tab, Time = proc.time()[3] - time.init))
        },
          stop(paste0('The response_type: ', response_type, " is't implemented"))
      )
    }

    Tab_Pred <- add_mean_amb(Tab_Pred, dec)

    if (verbose) {
      pb$tick(tokens = list(what = paste0(i, ' CV of ', nCV)))
      cat('Done.\n')
    }

    out <- list(
      predictions_Summary = Tab_Pred,
      cv = PT$cv,
      response = data$Response,
      predictions = data$Predictions
    )

    class(out) <- 'BFRCV'
    return(out)
  }else{
    return(BGLR(data$Response, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                      verbose, rmExistingFiles, groups))
  }
}

