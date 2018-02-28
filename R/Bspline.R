#' @title Bspline.Basis
#'
#' @description Solo es una prueba
#'
#' @param Bands \code{matrix}
#' @param Wavelenghts \code{vector}
#' @param n.basis \code{integer}
#' @param interaction \code{matrix}
#' @export
#'
Bspline.Basis <- function(Bands, Wavelengths = NULL, n.basis = 1, interaction = NULL){
  if (is.null(Wavelengths)) {
    stop("Wavelengths names was not provided")
  }

  bspl <- fda::create.bspline.basis(range(c(Wavelengths)), nbasis = n.basis, breaks = NULL, norder = 4)
  n.ind <- dim(Bands)[1]
  W.bs <- matrix(NA, nrow = n.ind, ncol = n.basis)
  for (h in 1:n.ind) {
    smf <- fda::smooth.basisPar(argvals = c(Wavelengths), y = as.numeric(Bands[h,]), lambda = 0, fdobj = bspl)
    cv_sp_pn <- smf$fd$coefs
    I_KL <- fda::inprod(bspl, bspl)
    xt_h <- t(I_KL %*% cv_sp_pn)
    W.bs[h,] <- xt_h
  }


  W.bs <- ifelse(is.null(interaction), W.bs, model.matrix(~0+interaction:W.bs)) #EnvxBand
  return(W.bs)
}
