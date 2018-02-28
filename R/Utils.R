#' @title Post-Processing
#'
#' @description Postprocesamiento
#' @param Tab \code{data.frame}
#'
#' @export
#'

Cor_Env <- function(Tab, Time){
  Envs <- unique(Tab$Env)
  Res <- data.frame(Fold = NA, Env = Envs, Pearson = NA, SE_Pearson = NA, MSEP = NA,  SE_MSEP = NA, Time = NA)
  Res$Time[1] <- Time
  for (i in 1:length(Envs)) {
    Tab_i <- Tab[Tab$Env == Envs[i],]
    Cor <- cor(Tab_i$y_p, Tab_i$y_o , use = "pairwise.complete.obs")
    MSEP <- mean((Tab_i$y_p - Tab_i$y_o)**2, na.rm = T)
    Res$Pearson[i] <- Cor
    Res$MSEP[i] <- MSEP
    Res$Fold[i] <- Tab$Fold[i]
  }

  # Res <- rbind(Res, data.frame(Fold = Tab$Fold[i], Env = "All", Cor = cor(Tab$y_p, Tab$y_o), MSEP =  mean((Tab$y_p-Tab$y_o)**2)))
  return(Res)
}

Cor_Env_Ordinal <- function(Tab, Folds = 1, Time){
  Envs <- unique(Tab$Env)
  Res <- data.frame(Fold = NA, Env = Envs, Pearson = NA, SE_Pearson = NA, MSEP = NA, SE_MSEP = NA, Time = NA)
  Res$Time[1] <- Time
  for (i in 1:length(Envs)) {
    Tab_i <- Tab[Tab$Env == Envs[i],]
    tabl <- table(Tab_i$y_p, Tab_i$y_o)
    prop.tabl <- prop.table(tabl)
    Cor <- sum(diag(prop.tabl))
    Res$Pearson[i] <- Cor
    Res$MSEP[i] <- NA
    Res$Fold[i] <- Tab$Fold[i]
  }
  return(Res)
}

saveFile <- function(Data, fname, rmExistingFiles=TRUE) {
  if (rmExistingFiles) {
    unlink(fname)
  }else{
      fname <- strsplit(fname,".", fixed = T)[[1]]
      fname <- append(fname,date(),after = length(fname) - 1)
  }
  save(Data, file = paste(fname, collapse = '.'))
}

add_mean_amb <- function(Tab_Pred, dec = 4){
  Envs <- unique(Tab_Pred$Env)
  Res <- data.frame(Fold = NA, Env = Envs, Pearson = NA, SE_Pearson = NA, MSEP = NA, SE_MSEP = NA, Time = NA)
  Res$Time[1] <- mean(Tab_Pred$Time, na.rm = T)
  for (i in seq_len(length(Envs))) {
    Res$Pearson[i] <- mean(Tab_Pred$Pearson[which(Tab_Pred$Env == Envs[i])], na.rm = T)
    Res$MSEP[i] <- mean(Tab_Pred$MSEP[which(Tab_Pred$Env == Envs[i])], na.rm = T)
    sd_Pearson <- sd(Tab_Pred$Pearson[which(Tab_Pred$Env == Envs[i])], na.rm = T)
    sd_MSEP <- sd(Tab_Pred$MSEP[which(Tab_Pred$Env == Envs[i])], na.rm = T)

    Res$SE_Pearson[i] <- sd_Pearson / sqrt(length(unique(Tab_Pred$Fold)))
    Res$SE_MSEP[i] <- sd_MSEP / sqrt(length(unique(Tab_Pred$Fold)))
    Res$Fold[i] <- "Average_all"
  }

  Tab_Pred$Pearson <- round(Tab_Pred$Pearson, dec)
  Tab_Pred$MSEP <- round(Tab_Pred$MSEP, dec)
  Res[, -(1:2)] <- round(Res[, -(1:2)], dec)
  return(rbind(Tab_Pred, Res))
}
