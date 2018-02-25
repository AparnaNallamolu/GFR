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
  X.FDA <- matrix(NA, nrow = n.ind, ncol = n.basis)
  for (h in 1:n.ind) {
    smf <- fda::smooth.basisPar(argvals = c(Wavelengths), y = as.numeric(Bands[h, ]),  lambda = 0.1, fdobj = bspl, Lfdobj = 2)
    I_KL <- fda::inprod(bspl, bspl)
    X.FDA[h,] <- t(I_KL %*% smf$fd$coefs)
  }

  X.FDA <- ifelse(is.null(interaction), X.FDA, model.matrix(~0+X.FDA:interaction))
  return(X.FDA)
}
