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
  X.Fu <- matrix(NA, nrow = n.ind,ncol = n.basis)
  for (h in 1:n.ind) {
    smf <- fda::smooth.basisPar(argvals = c(Wavelengths), y = as.numeric(Bands[h,]), lambda = 0.1, fdobj = bspF, Lfdobj = 2)
    cv_sp_pn <- smf$fd$coefs
    I_KL <- fda::inprod(bspF, bspF)
    xt_h <- t(I_KL %*% cv_sp_pn)
    X.Fu[h,] <- xt_h
  }

  X.Fu <- ifelse(is.null(interaction), X.Fu, model.matrix(~0+X.Fu:interaction))
  return(X.Fu)
}
