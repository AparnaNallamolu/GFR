#' FFR
#'
#' Frequentist Functional Regression
#'
#' More details about this function
#' @export
FFR <- function(dataset, datasetID = 'Line', X = NULL, Z = NULL, R = NULL, method = "NR", CrossValidation = NULL, init = NULL, iters = 20, tolpar = 1e-3,
                tolparinv = 1e-6, verbose = FALSE, constraint = TRUE, EIGEND = FALSE,
                forced = NULL, IMP = FALSE, complete = TRUE, check.model = TRUE, restrained = NULL,
                REML = TRUE, init.equal = TRUE, set_seed=NULL){
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

  if (!is.null(CrossValidation)) {

    switch(CrossValidation$Type,
           KFold = {
             if (is.null(CrossValidation$nFolds)) {message('Crossvalidation is used but nFolds is null, by default nFolds is set to 5.')}
             nCV <- CrossValidation$nFolds
           },

           RandomPartition = {
             if (is.null(CrossValidation$NPartitions)) {message('Crossvalidation is used but nFolds is null, by default nFolds is set to 10.')}
             PT <- CV.RandomPart(data, NPartitions = CrossValidation$NPartitions, PTesting = CrossValidation$PTesting, Traits.testing = CrossValidation$Traits.testing, set_seed)
             nCV <- CrossValidation$NPartitions
           },
           Error(paste0('ERROR: The Cross Validation  ', CrossValidation$Type, " is't implemented"))
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
      fm <- mmer(Y = dataset$Response, X = X, Z = Z, R = R, method = method, init = init, iters = iters, tolpar = tolpar,
                 tolparinv = tolparinv, verbose = verbose, constraint = constraint, EIGEND = EIGEND,
                 forced = forced, IMP = IMP, complete = complete, check.model = check.model, restrained = restrained,
                 REML = REML, init.equal = init.equal)

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
               Tab_Pred <- rbind(Tab_Pred, Cor_Env_Ordinal(Tab, Time = proc.time()[3] - time.init))
             },
             Error(paste0('The response_type: ', response_type, " is't implemented"))
      )
    }

    Tab_Pred <- add_mean_amb(Tab_Pred, dec)

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

    class(out) <- 'FFRCV'
    return(out)
  }else{
    return(mmer(Y = dataset$Response, X = X, Z = Z, R = R, method = method, init = init, iters = iters, tolpar = tolpar,
                tolparinv = tolparinv, verbose = verbose, constraint = constraint, EIGEND = EIGEND,
                forced = forced, IMP = IMP, complete = complete, check.model = check.model, restrained = restrained,
                REML = REML, init.equal = init.equal))
  }
}
