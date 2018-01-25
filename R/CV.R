#' @title Cross-Validation
#'
#' @description Example
#'
#' @param dat_F \code{data.frame}
#' @param K \code{integer}
#'
#' #rm(list=ls(all=TRUE))
#' #K Numero de grupos
#' #dat_F$Line and dat_F$Env
#' @export
#'
crossvalidation <- function(dat_F, K = 5, set_seed=NULL) {
  if(!is.null(set_seed)){
    set.seed(set_seed)
  }
  UL <- unique(dat_F$Line)
  #Number of sites where each line appear
  n_UL <- length(UL)
  nSLA <- rep(NA, n_UL)

  nEAL <- table(dat_F[, c('Line')])#Number of Sites that appear  each line
  L_nE <- data.frame(Line = names(nEAL), nE = c(nEAL))

  #A list of Positions in data set dat_F that will conform the groups
  g_list <- vector('list', K)
  names(g_list) <- paste('g', 1:K, sep = '')

  #Lines that will appear in all groups because
  # only appear in only one Site
  Pos1 <- which(L_nE$nE == 1)
  Pos_1_dat_F <- match(L_nE$Line[Pos1], dat_F$Line)
  #dat_F[Pos_1_dat_F,]

  #Tama?o de cada partici?n sin considerar las lineas
  # que se incluir?n por defaul (las que aparecen en un solo ambiente)
  n <- dim(dat_F)[1]
  nR <- n - length(Pos1)
  ifelse(nR %% K == 0, ng <- rep(nR / K, K), ng <-
           rep(trunc(nR / K), K) +
           c(rep(1, nR - trunc(nR / K) * K), rep(0, K - (nR - trunc(
           nR / K ) * K))))
  #ng
  Pos_all <- 1:n
  #---------------------------------------------------------------
  #First group
  #---------------------------------------------------------------
  if (length(Pos1) == 0){
    dat_F_k <- dat_F
  }
  else{
    dat_F_k <- dat_F[-Pos_1_dat_F,]
  }
  #Lineas ?nicas restantes
  UL_k <- unique(dat_F_k$Line)
  Pos_R_k <- rep(NA, length (UL_k))
  for (j in 1:length(UL_k)) {
    Pos_j_k <-  which(dat_F$Line == UL_k[j])
    Pos_R_k[j] <- sample(Pos_j_k, 1)
  }
  Pos_R_k <- Pos_R_k

  Pos_k_2_dat_F <- sample(Pos_all[-c(Pos_1_dat_F, Pos_R_k)], ng[1])#-length(Pos_1_dat_F))
  g_list[[1]] <- c(Pos_1_dat_F, Pos_k_2_dat_F)
  #---------------------------------------------------------------
  #Group 2,3, .., K
  #---------------------------------------------------------------
  for (k in 2:(K - 1))  {
    #Assigned positions
    Pos_k_a_R = unique(unlist(g_list[1:(k - 1)]))
    dat_F_k = dat_F[-Pos_k_a_R,]
    UL_k = unique(dat_F_k$Line)
    #A las lineas que no aparecen en el grupo k-1, se remueve un
    # site donde aparencen para garantizar que ?stas aparezcan
    # en al menos un site
    UL_k = UL_k[(UL_k %in% dat_F[Pos_k_a_R,]$Line) == FALSE]
    if (length(UL_k) > 0)
    {
      #Posiciones de lineas a mantener fuera del grupo k
      Pos_R_k = rep(NA, length (UL_k))

      for (j in 1:length(UL_k))
      {
        Pos_j_k =  which((dat_F$Line == UL_k[j]))#[-Pos_k_a_R])
        if (length(Pos_j_k) > 1)
        {
          Pos_R_k[j] = sample(Pos_j_k, 1)
        }
        #cat('L',length(Pos_j_k),'\n')
      }
      Pos_R_k = na.omit(Pos_R_k)
      Pos_k_2_dat_F = sample(Pos_all[-c(Pos_k_a_R, Pos_R_k)], ng[k])#-length(Pos_1_dat_F))
      g_list[[k]] = c(Pos_1_dat_F, Pos_k_2_dat_F)
    }
    else{
      Pos_k_2_dat_F = sample(Pos_all[-c(Pos_k_a_R)], ng[k])#-length(Pos_1_dat_F))
      g_list[[k]] = c(Pos_1_dat_F, Pos_k_2_dat_F)
    }
  }
  k = K
  Pos_k_a_R = unique(unlist(g_list[1:(k - 1)]))
  Pos_k_2_dat_F = sample(Pos_all[-c(Pos_k_a_R)], ng[k])#-length(Pos_1_dat_F))
  g_list[[k]] = c(Pos_1_dat_F, Pos_k_2_dat_F)

  n_CL = length(Pos_1_dat_F)
  list(g_Pos_ls = g_list,
       ng = ng + n_CL,
       n_CL =  n_CL)#nc_CL # Number of common lines

}
