#' @title Bspline.Basis
#'
#' @description Lorem ip sum
#'
#' @param Bands (\code{matrix}) Lorem ip sum
#' @param Wavelenghts (\code{integer}) Lorem ip sum
#' @param nBasis (\code{integer}) Lorem ip sum
#' @param interaction (\code{matrix}) Lorem ip sum
#' @export
#'
Bspline.Basis <- function(Bands, Wavelengths = NULL, nBasis = 1, interaction = NULL){
  if (is.null(Wavelengths)) {
    stop("Wavelengths names was not provided")
  }

  bspl <- fda::create.bspline.basis(range(c(Wavelengths)), nbasis = nBasis, breaks = NULL, norder = 4)
  n.ind <- dim(Bands)[1]
  W.bs <- matrix(NA, nrow = n.ind, ncol = nBasis)
  for (h in 1:n.ind) {
    smf <- fda::smooth.basisPar(argvals = c(Wavelengths), y = as.numeric(Bands[h,]), lambda = 0, fdobj = bspl)
    cv_sp_pn <- smf$fd$coefs
    I_KL <- fda::inprod(bspl, bspl)
    xt_h <- t(I_KL %*% cv_sp_pn)
    W.bs[h,] <- xt_h
  }


  if (is.null(interaction)) {
    return(W.bs)
  } else {
    return(model.matrix(~0+interaction:W.bs))#interactionxBand
  }
}
