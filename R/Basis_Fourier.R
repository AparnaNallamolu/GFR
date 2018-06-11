#' @title Fourier.Basis
#'
#' @description Lorem ip sum
#'
#' @param Bands (\code{matrix}) Lorem ip sum
#' @param Wavelengths (\code{integer}) Lorem ip sum
#' @param nBasis (\code{integer}) Lorem ip sum
#' @param interaction (\code{matrix}) Lorem ip sum
#' @export
#'
Fourier.Basis <- function(Bands, Wavelengths = NULL, nBasis = 1, interaction = NULL){
  if (is.null(Wavelengths)) {
    Error('Wavelengths names was not provided, use Wavelengths parameter to provided it')
  }

  bspF <- fda::create.fourier.basis(range(c(Wavelengths)), nbasis = nBasis, period = diff(range(c(Wavelengths))))
  n.ind <- dim(Bands)[1]
  W.bF <- matrix(NA, nrow = n.ind,ncol = nBasis)
  for (h in seq_len(n.ind)) {
    smf <- fda::smooth.basisPar(argvals = c(Wavelengths), y = as.numeric(Bands[h,]), lambda = 0, fdobj = bspF)
    cv_sp_pn <- smf$fd$coefs# Coeficientes cj directamente
    I_KL <- fda::inprod(bspF, bspF)
    xt_h <- t(I_KL %*% cv_sp_pn)
    W.bF[h,] <-  xt_h
  }

  if (is.null(interaction)) {
    return(W.bF)
  } else {
    return(model.matrix(~0+interaction:W.bF))
  }
}
