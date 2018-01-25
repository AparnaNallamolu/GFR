#' @title Simple Interaction Matrix
#'
#' @description example
#' @param X \code{matrix}
#'
#' @export
#'
simpleinteractionMatrix <- function(X){
  X <- as.matrix(X)
  return(model.matrix(~0+X))
}

#' @title double Interaction Matrix
#'
#' @description example
#' @param X1 \code{matrix}
#' @param X2 \code{matrix}
#' @export
#'
doubleinteractionMatrix <- function(X1,X2){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  return(model.matrix(~0+X1:X2))
}

#' @title triple Interaction Matrix
#'
#' @description example
#' @param X1 \code{matrix}
#' @param X2 \code{matrix}
#' @param X3 \code{matrix}
#' @export
#'
tripleInteractionMatrix <- function(X1,X2,X3){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  X3 <- as.matrix(X3)
  return(model.matrix(~0+X1:X2:X3))
}
