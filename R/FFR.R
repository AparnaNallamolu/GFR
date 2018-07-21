#' FFR
#'
#' Frequentist Functional Regression
#'
#' More details about this function
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach %dopar%
#' @export
FFR <- function(dataset, datasetID = 'Line', X = NULL, Z = NULL, R = NULL, method = "NR", testingSet = NULL, parallelCores = 1, init = NULL, iters = 20, tolpar = 1e-3,
                tolparinv = 1e-6, verbose = FALSE, constraint = TRUE, EIGEND = FALSE, forced = NULL, IMP = FALSE, complete = TRUE, check.model = TRUE, restrained = NULL,
                REML = TRUE, init.equal = TRUE, set_seed = NULL){
  if (inherits(Z, 'ETA')) {
    X <- Z$ETA[['X']]
    Z <- Z$ETA
    Z['X'] <- NULL
    design <- Z$Design
    datasetID <- Z$ID
    REML <- Z$REML
  } else {
    data <- validate.dataset(dataset, datasetID, orderData = F)
    design <- 'Handmade'
  }

  nCores <- detectCores()
  if(nCores < parallelCores){
    Message(paste0('The number of cores available is less than the specified in the parallelCores parameter we will use only ', nCores,'.'), verbose)
    parallelCores <- nCores
  }

  ########        INIT VALUES       #########
  # get_var_fm <- ifelse(is.factor(data$Response), 'fm$probs', 'fm$predictions')
  # data$Response <- as.numeric(data$Response)
  time.end <- c()

  if (is.null(testingSet)) {
    out <- mmer(Y = dataset$Response, X = X, Z = Z, R = R, method = method, init = init, iters = iters, tolpar = tolpar,
                tolparinv = tolparinv, verbose = verbose, constraint = constraint, EIGEND = EIGEND,
                forced = forced, IMP = IMP, complete = complete, check.model = check.model, restrained = restrained,
                REML = REML, init.equal = init.equal)
  } else if (parallelCores <= 1 && inherits(testingSet, 'CrossValidation')) {
    Tab_Pred <- data.frame()
    nCV <- length(testingSet$CrossValidation_list)
    pb <- progress::progress_bar$new(format = 'Fitting Cross-Validation :what  [:bar] :percent;  Time elapsed: :elapsed', total = nCV, clear = FALSE, show_after = 0)

    for (actual_CV in seq_len(nCV)) {
      if (progressBar) {
        pb$tick(tokens = list(what = paste0(actual_CV, ' out of ', nCV)))
      }

      positionTST <- testingSet$CrossValidation_list[[actual_CV]]
      response_NA <-  data$Response
      response_NA[positionTST] <- NA
      time.init <- proc.time()[3]
      fm <- mmer(Y = dataset$Response, X = X, Z = Z, R = R, method = method, init = init, iters = iters, tolpar = tolpar,
                 tolparinv = tolparinv, verbose = verbose, constraint = constraint, EIGEND = EIGEND,
                 forced = forced, IMP = IMP, complete = complete, check.model = check.model, restrained = restrained,
                 REML = REML, init.equal = init.equal)
      if(is.factor(data$Response)) {
        predictions <- as.integer(colnames(fm$probs)[apply(fm$probs,1,which.max)])
      } else{
        predictions <- fm$predictions
      }

      Tab_Pred <- rbind(Tab_Pred,
                        data.frame(Position = positionTST,
                                   Environment = data$Env[positionTST],
                                   Trait = data$Trait[positionTST],
                                   Partition = actual_CV,
                                   Observed = round(data$Response[positionTST], digits),
                                   Predicted = round(predictions[positionTST], digits)))
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

    class(out) <- 'FFRCV'
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
      fm <- mmer(Y = dataset$Response, X = X, Z = Z, R = R, method = method, init = init, iters = iters, tolpar = tolpar,
                 tolparinv = tolparinv, verbose = verbose, constraint = constraint, EIGEND = EIGEND,
                 forced = forced, IMP = IMP, complete = complete, check.model = check.model, restrained = restrained,
                 REML = REML, init.equal = init.equal)

      if(is.factor(data$Response)) {
        predictions <- as.integer(colnames(fm$probs)[apply(fm$probs,1,which.max)])
      } else{
        predictions <- fm$predictions
      }

      data.frame(
        Position = positionTST,
        Environment = data$Env[positionTST],
        Trait = data$Trait[positionTST],
        Partition = actual_CV,
        Observed = round(data$Response[positionTST], digits),
        Predicted = round(predictions[positionTST], digits)
      )
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

    class(out) <- 'FFRCV'
  } else {
    positionTST <- testingSet$CrossValidation_list[[actual_CV]]
    response_NA <-  data$Response
    response_NA[positionTST] <- NA
    time.init <- proc.time()[3]
    fm <- mmer(Y = dataset$Response, X = X, Z = Z, R = R, method = method, init = init, iters = iters, tolpar = tolpar,
               tolparinv = tolparinv, verbose = verbose, constraint = constraint, EIGEND = EIGEND,
               forced = forced, IMP = IMP, complete = complete, check.model = check.model, restrained = restrained,
               REML = REML, init.equal = init.equal)

    if(is.factor(data$Response)) {
      predictions <- as.integer(colnames(fm$probs)[apply(fm$probs,1,which.max)])
    } else{
      predictions <- fm$predictions
    }
    Tab_Pred <- rbind(Tab_Pred,
                      data.frame(Position = positionTST,
                                 Environment = data$Env[positionTST],
                                 Trait = data$Trait[positionTST],
                                 Partition = actual_CV,
                                 Observed = round(data$Response[positionTST], digits),
                                 Predicted = round(predictions[positionTST], digits)))

    out <- list(
      predictions_Summary = Tab_Pred,
      CrossValidation_list = testingSet$CrossValidation_list,
      response = data$Response,
      Design = design,
      nCores = parallelCores,
      response_type =  response_type
    )

    class(out) <- 'FFRCV'
  }
  return(out)

}
