#' @title Fourier.Basis
#'
#' @description
#'
#' @param Bands \code{matrix}
#' @param Wavelenghts \code{vector}
#' @param n.basis \code{integer}
#' @param interaction \code{matrix}
#' @export
#'
Fourier.Basis <- function(Bands, Wavelengths = NULL, n.basis = 1, interaction = NULL){
  if (is.null(Wavelengths)) {
    stop("Wavelengths names was not provided")
  }

  bspF <- fda::create.fourier.basis(range(c(Wavelengths)), nbasis = n.basis, period = diff(range(c(Wavelengths))))
  n.ind <- dim(Bands)[1]
  W.bF <- matrix(NA, nrow = n.ind,ncol = n.basis)
  for (h in 1:n.ind) {
    smf <- fda::smooth.basisPar(argvals = c(Wavelengths), y = as.numeric(Bands[h,]), lambda = 0, fdobj = bspF)
    cv_sp_pn <- smf$fd$coefs# Coeficientes cj directamente
    I_KL <- fda::inprod(bspF, bspF)
    xt_h <- t(I_KL%*%cv_sp_pn)
    W.bF[h,] <-  xt_h
  }

  W.bF <- ifelse(is.null(interaction), W.bF, model.matrix(~0+interaction:W.bF)) #EnvxBand
  return(W.bF)
}
