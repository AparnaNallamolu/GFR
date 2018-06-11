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
#'
#' @export

# DATA = data.frame()
# multivariateAdjust or multivariateEstimation
# REMOVE response_type (Implementar en data)

BFR <- function(data = NULL, datasetID = 'Line',  Multivariate = "Traditional", a=NULL, b=NULL, ETA = NULL, nIter = 1500, burnIn = 500, thin = 5,
                saveAt = '', S0 = NULL, df0 = 5, R2 = 0.5, weights = NULL, verbose = TRUE, rmExistingFiles = TRUE, groups = NULL,
                CrossValidation = NULL, set_seed = NULL, dec = 4) {
  if (inherits(ETA, 'ETA')) {
    data <- ETA$data
    design <- ETA$Design
    datasetID <- ETA$ID
    ETA <- ETA$ETA
  } else {
    data <- validate.dataset(data, datasetID, orderData = F, Multivariate = Multivariate)
    design <- 'Handmade'
  }

  response_type <- ifelse(is.factor(data$Response), 'ordinal', 'gaussian')
  data$Response <- as.numeric(data$Response)

  Multivariate <- validate.MTME(Multivariate)
  if (!is.null(CrossValidation) && Multivariate == 'Traditional') {

    switch(CrossValidation$Type,
           KFold = {
             if (is.null(CrossValidation$nFolds)) {message('Crossvalidation is used but nFolds is null, by default nFolds is set to 5.')}
             PT <- CV.KFold(data, DataSetID = datasetID, K = CrossValidation$nFolds, set_seed)
             nCV <- CrossValidation$nFolds
           },
           RandomPartition = {
             if (is.null(CrossValidation$NPartitions)) {message('Crossvalidation is used but NPartitions is null, by default NPartitions is set to 10.')}
             PT <- CV.RandomPart(data, NPartitions = CrossValidation$NPartitions, PTesting = CrossValidation$PTesting, Traits.testing = CrossValidation$Traits.testing, set_seed)
             nCV <- CrossValidation$NPartitions
           },
           Error(paste0('The Cross Validation  ', CrossValidation$Type, " is't implemented"))
    )

    if (verbose) {
      Message("This might be time demanding...")

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
      Pos_NA <- PT$CrossValidation_list[[paste0('partition',i)]]
      response_NA[Pos_NA] <- NA
      time.init <- proc.time()[3]
      fm <- BGLR(response_NA, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                       verbose = F, rmExistingFiles, groups)

      #fm <- BGLR(response_NA, ETA, verbose=F)
      switch(response_type,
        gaussian = {
          predicted <- fm$predictions
          data$Predictions[Pos_NA] <- predicted[Pos_NA]
          Tab <- data.frame(Env = data$Env[Pos_NA], Trait = data$Trait[Pos_NA], Fold = i,
                            y_p = predicted[Pos_NA], y_o = data$Response[Pos_NA])
          Tab_Pred <- rbind(Tab_Pred, Cor_Env(Tab, Time = proc.time()[3] - time.init ))
        },
        ordinal = {
          predicted <- as.integer(colnames(fm$probs)[apply(fm$probs,1,which.max)])
          data$Predictions[Pos_NA] <- predicted[Pos_NA]

          Tab <- data.frame(Env = data$Env[Pos_NA], Trait = data$Trait[Pos_NA], Fold = i,
                            y_p = predicted[Pos_NA], y_o = data$Response[Pos_NA] )
          Tab_Pred <- rbind(Tab_Pred, Cor_Env_Ordinal(Tab, Time = proc.time()[3] - time.init, nFolds = nCV))
        },
          Error(paste0('The response_type: ', response_type, " is't implemented"))
      )
    }
    if (response_type == 'gaussian') {
      Tab_Pred <- add_mean_amb(Tab_Pred, dec)
    } else {
      Tab_Pred <- add_mean_amb_Ordinal(Tab_Pred, dec)
    }

    if (verbose) {
      pb$tick(tokens = list(what = paste0(i, ' CV of ', nCV)))
      cat('Done.\n')
    }

    out <- list(
      predictions_Summary = Tab_Pred,
      CrossValidation_list = PT$CrossValidation_list,
      response = data$Response,
      predictions = data$Predictions,
      Design = design
    )

    class(out) <- 'BFRCV'
    return(out)
  } else if (!is.null(CrossValidation) && Multivariate == 'SVD') {

    switch(CrossValidation$Type,
           KFold = {
             if (is.null(CrossValidation$nFolds)) {message('Crossvalidation is used but nFolds is null, by default nFolds is set to 5.')}
             PT <- CV.KFold(data[,1:3], K = CrossValidation$nFolds, set_seed)
             nCV <- CrossValidation$nFolds
           },
           RandomPartition = {
             if (is.null(CrossValidation$NPartitions)) {message('Crossvalidation is used but NPartitions is null, by default NPartitions is set to 10.')}
             PT <- CV.RandomPart(data[,1:3], NPartitions = CrossValidation$NPartitions, PTesting = CrossValidation$PTesting, Traits.testing = CrossValidation$Traits.testing, set_seed)
             nCV <- CrossValidation$NPartitions
           },
           Error(paste0('The Cross Validation  (', CrossValidation$Type, ") is't implemented"))
    )

    if (verbose) {
      Message("This might be time demanding...")
      pb <- progress::progress_bar$new(format = 'Fitting the :what  [:bar] Time elapsed: :elapsed', total = nCV + 1, clear = FALSE, show_after = 0)
    }

    Criteria <- data.frame(matrix(0, ncol = 3, nrow = 9 * nCV))  #Cada columna
    Env.ALL <- data.frame(matrix(NA, nrow = 9, ncol = 5))
    ## Init cross validation
    for (i in seq_len(nCV)) {
      time.init <- proc.time()[3]

      if (verbose) {
        pb$tick(tokens = list(what = paste0( i, ' CV of ', nCV)))
      }

      response <-  as.matrix(data[,-c(1,2)]) #Ignore Line and Env
      nt <- ncol(response)
      nI <- length(unique(data$Env))
      Pos_NA <- PT$CrossValidation_list[[paste0('partition',i)]]
      training_data <- response[-Pos_NA, ]
      SVD_Y <- svd(training_data)
      U <- SVD_Y$u
      V <- SVD_Y$v
      tV <- t(V)

      Ytilde <- response %*% V
      Ytilde[Pos_NA,] <- rep(NA, ncol(response))

      Y_pred <- matrix(NA, nrow = nrow(response), ncol = nt)
      Beta_PC <- matrix(NA, nrow = nI+1, ncol = nt)
      SDBeta_PC <- matrix(NA, nrow = nI+1, ncol = nt)
      Sigma1_PC <- matrix(0, nrow = nt+1, ncol = nt)
      Sigma2_PC <- matrix(0, nrow = nt+1, ncol = nt)
      Sigma3_PC <- matrix(0, nrow = nt+1, ncol = nt)
      SigmaError_PC <- matrix(0, nrow = nt+1, ncol = nt)
      SD1_PC <- matrix(0, nrow = nt+1, ncol = nt)
      SD2_PC <- matrix(0, nrow = nt+1, ncol = nt)
      SD3_PC <- matrix(0, nrow = nt+1, ncol = nt)
      SDError_PC <- matrix(0, nrow = nt+1, ncol = nt)

      for (i in seq_len(nt)) {
        y_i <- Ytilde[, i]
        fm <- BGLR(y_i, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                   verbose = F, rmExistingFiles, groups)

        yhat2 <- fm$predictions
        betas_est <- fm$mu + fm$ETA[[1]]$b
        SDbetas_est <- fm$SD.mu + fm$ETA[[1]]$SD.b
        Y_pred[, i] <- yhat2
        Beta_PC[, i] <- c(betas_est, fm$mu)
        SDBeta_PC[, i] <- c(SDbetas_est, fm$SD.mu)
        Sigma1_PC[i, i] <- fm$ETA[[2]]$varU
        SD1_PC[i, i] <- fm$ETA[[2]]$SD.varU
        Sigma2_PC[i, i] <- fm$ETA[[3]]$varU
        SD2_PC[i, i] <- fm$ETA[[3]]$SD.varU
        # Sigma3_PC[i, i] <- fm$ETA[[4]]$varU
        # SD3_PC[i, i] <- fm$ETA[[4]]$SD.varU
        SigmaError_PC[i, i] <- fm$varE
        SDError_PC[i, i] <- fm$SD.varE

        # switch(response_type,
        #        gaussian = {
        #          predicted <- fm$predictions
        #          # data$Predictions[Pos_NA] <- predicted[Pos_NA]
        #          Tab <- data.frame(Env = data$Env[Pos_NA], Trait = data$Trait[Pos_NA], Fold = i,
        #                            y_p = predicted[Pos_NA], y_o = data$Response[Pos_NA])
        #          Tab_Pred <- rbind(Tab_Pred, Cor_Env(Tab, Time = proc.time()[3] - time.init ))
        #        },
        #        stop(paste0('The response_type: ', response_type, " is't implemented for SVD multivariate method"))
        # )
      }

      Y_pred_Final <- Y_pred %*% tV


      Y_all = cbind(as.matrix(data[,-c(1,2)]), Y_pred_Final, data$Env)
      Y_all_tst = Y_all[Pos_NA, ]
      Y_all_tst
      Data_pred = data.frame(matrix(NA, ncol = 3, nrow = nt * nI))
      Traits = noquote(colnames(data[,-c(1,2)]))
      Env = unique(data$Env)
      Etiquetas = expand.grid(Trait = Traits, Env = Env)

      Data_pred[, 1] = paste(Etiquetas$Trait, Etiquetas$Env, sep = "_")
      for (j in 1:nI) {
        #for (r in 1:nI){
        Env_i = Y_all_tst[Y_all_tst[, 2 * nt + 1] == j, ]
        Cor_all_i = cor(Env_i[, -(2 * nt + 1)])
        cor_U = diag(Cor_all_i[1:(nt), (nt + 1):(2 * nt)])
        PL = numeric()
        for (k in 1:nt) {
          PL = c(PL, mean((Env_i[, k] - Env_i[, nt + k]) ^ 2))
        }
        Pos_L = PL
        for (s in 1:nt) {
          index = nt * (j - 1) + nt
          Data_pred[(((index - nt) + 1):index), 2] = cor_U
          Data_pred[(((index - nt) + 1):index), 3] = c(Pos_L)
        }
      }
      Data_pred
      colnames(Data_pred) = c("Group", "Cor", "MSEP")

      indexN = 9 * (i - 1) + 9
      Criteria[(((indexN - 9) + 1):indexN),] = Data_pred
    }

    Criteria
    GG = noquote(Data_pred[, 1])
    for (i in 1:length(GG)) {
      Criteria.Env_o = Criteria[Criteria[, 1] == GG[i], ]
      Env.ALL[i, ] = c(
        GG[i],
        mean(Criteria.Env_o[, 2], na.rm = TRUE),
        (sd(Criteria.Env_o[, 2], na.rm = TRUE) / sqrt(20)),
        mean(Criteria.Env_o[, 3], na.rm = TRUE),
        (sd(Criteria.Env_o[, 3], na.rm = TRUE) / sqrt(20))
      )

    }
    # Env.ALL <- tidyr::separate(Env.ALL, 'Trait_Env', c("Trait", "Env"), sep = "_")

    if (verbose) {
      pb$tick(tokens = list(what = paste0(i, ' CV of ', nCV)))
      cat('Done.\n')
    }

    out <- list(
      predictions_Summary = Env.ALL,
      CrossValidation_list = PT$CrossValidation_list,
      response = data[, -c(1,2)],
      predictions = Y_pred_Final,
      Design = design,
      Beta_est = betas_est,
      SDBeta_est = SDbetas_est,
      Beta_PC = Beta_PC,
      SDBeta_PC = SDBeta_PC,
      Sigma1_PC = Sigma1_PC,
      SD1_PC = SD1_PC,
      Sigma2_PC = Sigma2_PC,
      SD2_PC = SD2_PC,
      SigmaError_PC = SigmaError_PC,
      SDError_PC = SDError_PC
    )

    class(out) <- 'BFRCV'
    return(out)
  } else {
    return(BGLR(data$Response, response_type, a, b, ETA, nIter, burnIn, thin, saveAt, S0, df0, R2, weights,
                      verbose, rmExistingFiles, groups))
  }
}

validate.MTME <- function(Multivariate){
  if (Multivariate %in% c('SVD', 'Traditional')) {
    return(Multivariate)
  } else {
    Error('The Multivariate parameter not found')
  }
}
