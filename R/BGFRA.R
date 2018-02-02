#' @title BGFRA
#'
#' @description Bayesian Genomic Functional Regression Analysis
#'
#'
#' @param data (\code{data.frame}) the data with the $n$ $Response, also needs $Line and $Env for Cross Validation defined on it (NAs allowed).
#' @param response_type (\code{string}) It can be "gaussian" or "ordinal".
#' @param a (\code{numeric}, $n$) Only requiered for censored outcomes. It's a vector specifying lower bounds for censored observation. By default is null.
#' @param b (\code{numeric}, $n$) Only requiered for censored outcomes. It's a vector specifying upper bounds for censored observation. By default is null.
#' @param ETA (\code{list}) Two level list used to specify the regression function.
#' @param nIter (\code{integer}) The number of iterations.
#' @param burnIn (\code{integer}) The number of burn-in.
#' @param thin (\code{integer}) The number of thinning.
#' @param saveAt (\code{string}) This may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.
#' @param S0 (\code{numeric}) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes.
#' @param df0 (\code{numeric}) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes.
#' @param R2 (\code{numeric},$0<R2<1$) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified.
#' @param weights (\code{numeric}, $n$) a vector of weights, may be NULL. The residual variance of each data-point is set to be proportional to the inverse of the squared-weight. Only used with Gaussian outcomes.
#' @param verbose (\code{logical}) By default is \code{TRUE} and shows a fitting model progress bar if Folds <=1 or cross validation progress bar if Folds > 2.
#' @param rmExistingFiles (\code{logical}) By default is \code{TRUE} and removes existing output files from previous runs.
#' @param groups (\code{factor}) A vector of the same lenght of \code{data$response} that associates observations with groups, each group will have an associated variance component for the error term.
#' @param folds (\code{integer}) A $n$ number of the cross validations.
#' @param set_seed (\code{integer}) A seed to replicable research.
#'
#' @details BGFRA is an modificated version of BGLR which implements a Gibbs sampler for a Bayesian regression model, this new version allows to
#'
#'
#' @seealso \code{\link[BGLR]{BGLR}}
#'
#' examples
#'
#' @export

BGFRA <- function(data, response_type = "gaussian", a=NULL, b=NULL, ETA = NULL, nIter = 1500,
                  burnIn = 500, thin = 5, saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                  verbose = TRUE, rmExistingFiles = TRUE, groups=NULL, folds=1, set_seed=NULL){
  if (folds>1) {
    if (is.null(dim(data)[2])||dim(data)[2]<2) {
      stop("To realice Fold Cross-validation BGFRA requieres data param like a data.frame with $Response, $Line and optional the $Env specified on it")
    }

    cat("This might be time demanding, let's take sit and a cup of coffe\n")
    pb <- progress::progress_bar$new(format = ":what [:bar] Time elapsed: :elapsed", total = folds, clear = FALSE)

    ## Partitions
    PT <- crossvalidation(data, folds, set_seed)
    saveFile(PT, paste0(saveAt, "crossValidation_Partitions.RData"),rmExistingFiles)

    data$Predictions <- NA
    Tab_Pred <- data.frame()

    ## Init cross validation
    for (i in 1:folds) {
      pb$tick(tokens = list(what = paste0("Fitting the model - ", i, " CV of ", folds)))

      y_NA <-  data$response
      Pos_NA <- PT$cv[[paste0('partition',i)]]
      y_NA[Pos_NA] <- NA

      fm <- BGLR(y_NA, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                       verbose = F, rmExistingFiles, groups)

      saveFile(fm, paste0(saveAt, "FittedModel_CValidation_Fold-",i,".RData"))

      if (response_type == "gaussian") {
        yHat <- fm$yHat
        data$Predictions <- yHat[Pos_NA]

        if (is.null(data$Env)) {
          Tab_Pred <- rbind(Tab_Pred, data.frame(Fold = i,
                            Cor = cor(yHat[Pos_NA],data$response[Pos_NA]),
                            mean((yHat[Pos_NA]-data$response[Pos_NA])**2)))
        }else{
          Tab <- data.frame(Env = data$Env[Pos_NA], Fold=i, y_p = yHat[Pos_NA], y_o = data$response[Pos_NA] )
          Tab_Pred <- rbind(Tab_Pred, Cor_Env(Tab))
        }

      } else{
        stop(paste0(response_type, " - Not implemented yet"))
      }
    }

    if (!is.null(data$Env)) {
      Tab_Pred <- add_mean_amb(Tab_Pred)
    }

    saveFile(Tab_Pred, paste0(saveAt, "Results.RData"))
    cat("\nDone.")

    out <- list(
      results = Tab_Pred,
      cv = PT$cv,
      response = data$response,
      NA_predictions = data$Predictions
    )

    class(out) <- "BGFRA-CV"
    return(out)
  }else{
    return(BGLR(data$response, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                      verbose, rmExistingFiles, groups))
  }
}

