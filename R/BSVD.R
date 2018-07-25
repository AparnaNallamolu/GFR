BSVD <- function(data = NULL, datasetID = 'Line', GenomicMatrix = NULL, priorType = 'BRR', nIter = 1500, burnIn = 500, thin = 5,
                 saveAt = '', verbose = TRUE,progressBar = TRUE, CrossValidation = NULL, set_seed = NULL, digits = 4) {
  if (progressBar) {
    Message("This might be time demanding...", verbose)
    pb <- progress::progress_bar$new(format = 'Fitting Cross-Validation :what [:bar] :percent; Time elapsed: :elapsed', total = nCV, clear = FALSE, show_after = 0)
  }

  Criteria <- data.frame(matrix(0, ncol = 3, nrow = 9 * nCV))  #Cada columna
  Env.ALL <- data.frame(matrix(NA, nrow = 9, ncol = 5))


  priorType <- validate.prior(priorType)

  Env_prior <- ifelse(length(unique(data$Env)) < 10, "FIXED", priorType)

  if (!is.null(GenomicMatrix)) {
    L <- t(chol(GenomicMatrix))
    XL <- model.matrix(~0+as.factor(data[, datasetID])) %*% L
    XL <- kronecker(X = GenomicMatrix, Y = diag(length(unique(data$Env))))
    rownames(XL) <- data[, datasetID]
    XL <- validate.GenomicMatrix(XL, data[, datasetID])
    Line_prior <- 'RKHS'
    LineK <- TRUE
  } else {
    XL <- model.matrix(~0+as.factor(data[, datasetID]))
    Line_prior <- priorType
    LineK <- FALSE
  }

  data <- tidyr::spread(data, 'Trait', 'Response')

  ETA <- list(Env = ETAList(data$Env, Env_prior),
              Line = ETAList(XL, Line_prior, TRUE, LineK),
              LinexEnv = ETAList(XL, priorType, interaction1 = data$Env, likeKernel = TRUE, withK = T))

  ## Init cross validation
  for (i in seq_len(nCV)) {
    time.init <- proc.time()[3]

    if (progressBar) {
      pb$tick(tokens = list(what = paste0(actual_CV, ' out of ', nCV)))
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
}
