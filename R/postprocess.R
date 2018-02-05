#' @title Post-Processing
#'
#' @description Postprocesamiento
#' @param Tab \code{data.frame}
#'
#' @export
#'

Cor_Env <- function(Tab){
  Envs = unique(Tab$Env)
  Res = data.frame(Fold = NA, Env = Envs, Cor = NA, MSEP = NA)
  for (i in 1:length(Envs)) {
    Tab_i = Tab[Tab$Env == Envs[i],]
    Res$Cor[i] = cor(Tab_i$y_p,Tab_i$y_o)
    Res$MSEP[i] = mean((Tab_i$y_p - Tab_i$y_o)**2)
    Res$Fold[i] = Tab$Fold[i]
  }

  # Res <- rbind(Res, data.frame(Fold = Tab$Fold[i], Env = "All", Cor = cor(Tab$y_p, Tab$y_o), MSEP =  mean((Tab$y_p-Tab$y_o)**2)))
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

add_mean_amb <- function(Tab_Pred){
  Envs = unique(Tab_Pred$Env)
  Res = data.frame(Fold = NA, Env = Envs, Cor = NA, MSEP = NA)
  if (is.na(Envs[1])) {
    Res$Cor[1] <- mean(Tab_Pred$Cor)
    Res$MSEP[1] <- mean(Tab_Pred$MSEP)
    Res$Fold[1] <- "Average_all"
  }else{
    for (i in 1:length(Envs)) {
      Res$Cor[i] <- mean(Tab_Pred$Cor[which(Tab_Pred$Env == Envs[i])])
      Res$MSEP[i] <- mean(Tab_Pred$MSEP[which(Tab_Pred$Env == Envs[i])])
      Res$Fold[i] <- "Average_all"
    }
  }
  return(rbind(Tab_Pred,Res))
}

