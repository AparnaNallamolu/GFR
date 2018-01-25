#' @title BGFRA
#'
#' @description Bayesian Genomic Functional Regression Analysis
#'
#' @details This is a generic function
#'
#' @param y \code{data.frame} the data with $Response, $Line and $Env, defined on it (NAs allowed).
#' @param response_type \code{string} can be "gaussian" or "ordinal"
#' @param a,b \code{numeric} .
#' @param ETA \code{list} two level list.
#' @param weights \code{numeric}
#' @param nIter,burnIN,thin \code{integer}
#' @param saveAt \code{string}
#' @param S0,df0 \code{numeric}
#' @param R2 \code{numeric}
#' @param verbose \code{logical}
#' @param rmExistingFiles \code{logical}
#' @param groups \code{factor}
#' @param Folds \code{integer}
#' @param set.seed \code{integer} to reproducible research
#'
#'
#'
#' @export

BGFRA <- function(y, response_type = "gaussian", a=NULL, b=NULL, ETA = NULL, nIter = 1500,
                 burnIn = 500, thin = 5, saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                 verbose = FALSE, rmExistingFiles = TRUE, groups=NULL, Folds=0, set_seed=NULL){

  if (Folds>0) {
    if (is.null(dim(y)[2])||dim(y)[2]<2) {
      stop("To realice Fold Cross-validation BGFRA requieres y param like a data.frame with Line and Env specified on it")
    }

    PT <- crossvalidation(y, Folds, set_seed)

    saveFile(PT, paste0(saveAt, "crossValidation_Partitions.RData"),rmExistingFiles)
    Tab_Pred = data.frame()

    for (i in 1:Folds) {
      y_NA <-  y$response
      Pos_NA = PT$g_Pos_ls[[paste('g',i,sep='')]]
      y_NA[Pos_NA] = NA
      fm <- BGLR::BGLR(y_NA, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                 verbose = F, rmExistingFiles, groups)

      saveFile(fm, paste0(saveAt, "FittedModel_CValidation_Fold-",i,".RData"))

      if (response_type == "gaussian") {
        yHat = fm$yHat
        Tab = data.frame(Env = y$Env[Pos_NA], Fold=i, y_p = yHat[Pos_NA], y_o = y$response[Pos_NA] )
        Tab_Pred = rbind(Tab_Pred, Cor_MSEP_Env_f(Tab))
      } else{
        stop(paste0(response_type, " - Not implemented yet"))
      }
    }
    Tab_Pred <- add_mean_amb(Tab_Pred)
    saveFile(Tab_Pred, paste0(saveAt, "Results.RData"))

    return(Tab_Pred)
  }else{
    return(BGLR::BGLR(y$response, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                    verbose, rmExistingFiles, groups))
  }
}

