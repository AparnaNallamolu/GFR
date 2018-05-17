#' cleanDat function
#'
#' Easy way to delete dat files
#'
#' @param forceClean By default is FALSE, else, the function not requiere the confirmation by the user to delete de .dat files.
#'
#' @return
#' logical vector, of the files deleted.
#'
#'
#' @export
cleanDat <- function(noConfirm = F) {
  path <- dir()
  files <- path[which(grepl('.dat', path) == TRUE)]
  if (length(files) != 0) {
    if (noConfirm) {
      Status <- file.remove(files)
      names(Status) <- files
      return(Status)
    }
    Message("Are you sure that you want to delete this files?: y/n ")
    cat(files, sep = '\t')
    response <- scan(what = character(length = 1),n = 1,quiet = TRUE)
    if (response == 'y') {
      Status <- file.remove(files)
      names(Status) <- files
      return(Status)
    } else {
      Message("The deletion was cancelled by the user.")
    }
  } else{
    Message("No .dat files found in this directory")
  }
}


Cor_Env <- function(Tab, Time){
  Envs <- unique(Tab$Env)
  Traits <- unique(Tab$Trait)
  com <- expand.grid(Envs, Traits)
  Res <- data.frame(Fold = NA, Env = com[,1], Trait = com[, 2], Pearson = NA, SE_Pearson = NA, MSEP = NA,  SE_MSEP = NA, Time = NA)
  Res$Time[1] <- Time
  for (i in seq_len(length(com[, 1]))) {
    pos_env_trait <- intersect(which(Tab$Env == com[i,1]), which(Tab$Trait == com[i, 2]))

    Tab_i <- Tab[pos_env_trait, ]
    Res$Pearson[i] <- cor(Tab_i$y_p, Tab_i$y_o , use = "pairwise.complete.obs")
    Res$MSEP[i] <- mean((Tab_i$y_p - Tab_i$y_o)**2, na.rm = T)
    Res$Fold[i] <- Tab$Fold[i]
  }

  # Res <- rbind(Res, data.frame(Fold = Tab$Fold[i], Env = "All", Cor = cor(Tab$y_p, Tab$y_o), MSEP =  mean((Tab$y_p-Tab$y_o)**2)))
  return(Res)
}


add_mean_amb <- function(Tab_Pred, dec = 4){
  Envs <- unique(Tab_Pred$Env)
  Traits <- unique(Tab_Pred$Trait)

  com <- expand.grid(Envs, Traits)
  Res <- data.frame(Fold = NA, Env = com[,1], Trait = com[, 2], Pearson = NA, SE_Pearson = NA, MSEP = NA,  SE_MSEP = NA, Time = NA)
  Res$Time[1] <- mean(Tab_Pred$Time, na.rm = T)
  for (i in seq_len(length(com[, 1]))) {
    pos_env_trait <- intersect(which(Tab_Pred$Env == com[i,1]), which(Tab_Pred$Trait == com[i, 2]))

    Res$Pearson[i] <- mean(Tab_Pred$Pearson[pos_env_trait], na.rm = T)
    Res$MSEP[i] <- mean(Tab_Pred$MSEP[pos_env_trait], na.rm = T)
    sd_Pearson <- sd(Tab_Pred$Pearson[pos_env_trait], na.rm = T)
    sd_MSEP <- sd(Tab_Pred$MSEP[pos_env_trait], na.rm = T)

    Res$SE_Pearson[i] <- sd_Pearson / sqrt(length(unique(Tab_Pred$Fold)))
    Res$SE_MSEP[i] <- sd_MSEP / sqrt(length(unique(Tab_Pred$Fold)))
    Res$Fold[i] <- "Average_all"
  }

  Tab_Pred$Pearson <- round(Tab_Pred$Pearson, dec)
  Tab_Pred$MSEP <- round(Tab_Pred$MSEP, dec)
  Res[, -(1:3)] <- round(Res[, -(1:3)], dec)
  return(rbind(Tab_Pred, Res))
}

Cor_Env_Ordinal <- function(Tab, Time, nFolds){
  Envs <- unique(Tab$Env)
  Traits <- unique(Tab$Trait)
  com <- expand.grid(Envs, Traits)
  Res <- data.frame(Fold = NA, Env = com[,1], Trait = com[, 2], PCC = NA, SE_PCC = NA, MSEP = NA,  SE_MSEP = NA, Time = NA)
  Res$Time[1] <- Time
  for (i in seq_len(length(com[, 1]))) {
    pos_env_trait <- intersect(which(Tab$Env == com[i,1]), which(Tab$Trait == com[i, 2]))

    Tab_i <- Tab[pos_env_trait, ]
    tabl <- table(factor(Tab_i$y_p, levels = sort(unique(Tab_i$y_o))), as.factor(Tab_i$y_o))
    prop.tabl <- prop.table(tabl)
    Res$PCC[i] <- sum(diag(prop.tabl))
    Res$MSEP[i] <- mean((Tab_i$y_p - Tab_i$y_o)**2, na.rm = T)
    Res$Fold[i] <- Tab$Fold[i]
  }
  return(Res)
}


add_mean_amb_Ordinal <- function(Tab_Pred, dec = 4){
  Envs <- unique(Tab_Pred$Env)
  Traits <- unique(Tab_Pred$Trait)

  com <- expand.grid(Envs, Traits)
  Res <- data.frame(Fold = NA, Env = com[,1], Trait = com[, 2], PCC = NA, SE_PCC = NA, MSEP = NA,  SE_MSEP = NA, Time = NA)
  Res$Time[1] <- mean(Tab_Pred$Time, na.rm = T)
  for (i in seq_len(length(com[, 1]))) {
    pos_env_trait <- intersect(which(Tab_Pred$Env == com[i,1]), which(Tab_Pred$Trait == com[i, 2]))
    PCC <- mean(Tab_Pred$PCC[pos_env_trait], na.rm = T)
    Res$PCC[i] <- PCC
    Res$MSEP[i] <- mean(Tab_Pred$MSEP[pos_env_trait], na.rm = T)
    sd_MSEP <- sd(Tab_Pred$MSEP[pos_env_trait], na.rm = T)

    Res$SE_PCC[i] <- sqrt((PCC*(1 - PCC))/length(unique(Tab_Pred$Fold)))
    Res$SE_MSEP[i] <- sd_MSEP / sqrt(length(unique(Tab_Pred$Fold)))
    Res$Fold[i] <- "Average_all"
  }

  Tab_Pred$PCC <- round(Tab_Pred$PCC, dec)
  Tab_Pred$MSEP <- round(Tab_Pred$MSEP, dec)
  Res[, -(1:3)] <- round(Res[, -(1:3)], dec)
  return(rbind(Tab_Pred, Res))
}
