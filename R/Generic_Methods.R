#' @title Summary.BFR
#'
#' @description Solo es una prueba
#'
#' @param object \code{BFR object} Objeto BFR, resultado de ejecutar BFR
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom stats cor
#'
#' @export
summary.BFR <- function(object,...){

    if (!inherits(object, "BFR")) stop("This function only works for objects of class 'BFR'")

    tmp <- paste('--------------------> Summary of data & model <--------------------')
    cat(tmp,'\n\n')

    tmp <- paste(' Number of phenotypes=', sum(!is.na(object$response)))
    cat(tmp,'\n')

    cat(' Min (TRN)= ', min(object$response,na.rm = TRUE),'\n')
    cat(' Max (TRN)= ', max(object$response,na.rm = TRUE),'\n')
    cat(' Variance of phenotypes (TRN)=', round(var(object$response, na.rm = TRUE),4), '\n')
    cat(' Residual variance=', round(object$varE, 4), '\n')

    n <- length(object$response)

    if (any(is.na(object$response))) {
     		tst <- which(is.na(object$response))
     		cat(' N-TRN=', n - length(tst), ' N-TST=', length(tst), '\n')
     		cat(' Correlation TRN=', round(cor(object$response[-tst], object$predictions[-tst]), 4),'\n')
   }else{
       cat(' N-TRN=',n,'  N-TST=0', '\n\n')
   }

   cat('\n')
   cat(' -- Linear Predictor -- \n')
   cat('\n')
   cat(' Intercept included by default\n')

  for (k in 1:length(object$ETA)) {
    if (object$ETA[[k]]$model == "FIXED") {
      if (!is.null(names(object$ETA)[k])) {
        cat(" Coefficientes in ETA[",k,"] (",names(object$ETA)[k],") were asigned a flat prior\n")
      }else{
        cat(" Coefficientes in ETA[",k,"] (no-name) are asigned a flat prior\n")
      }
    }else{
      if (object$ETA[[k]]$model == "RKHS") {
        if (!is.null(names(object$ETA)[k])) {
          cat(" Coefficientes in ETA[",k,"] (",names(object$ETA)[k],") were assumed to be normally distributed with zero mean and \n covariance (or its eigendecoposition) provided by user \n")
        }else{
          cat(" Coefficientes in ETA[",k,"] (no-name) were assumed to be normally distributed with zero mean and \n covariance (or its eigendecoposition) provided by user \n")
        }
      }else{
        if (!is.null(names(object$ETA)[k])) {
				  cat(" Coefficientes in ETA[",k,"] (",names(object$ETA)[k],") modeled as in ", object$ETA[[k]]$model,"\n")
			  }else{
				cat(" Coefficientes in ETA[",k,"] (no-name) modeled as in ", object$ETA[[k]]$model,"\n")
			  }
		  }
	  }
  }

   cat('\n------------------------------------------------------------------\n');
}

#' @title Summary
#'
#' @description Summary of BFRCV object
#'
#' @param object \code{BFRCV object} BFRCV object, result of use the BFR() function
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
summary.BFRCV <- function(object,...){
  if (!inherits(object, "BFRCV")) stop("This function only works for objects of class 'BFRCV'")
  return(object$predictions_Summary)
}



#' @title residuals.BFR
#'
#' @description Solo es una prueba
#'
#' @param object \code{BFR object} Objeto BFR, resultado de ejecutar BFR
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
residuals.BFR <- function(object,...) {
    if (!inherits(object, "BFR")) stop("This function only works for objects of class 'BFR'")
	object$response - object$predictions
}

#' @title predict.BFR
#'
#' @description Solo es una prueba
#'
#' @param object \code{BFR object} Objeto BFR, resultado de ejecutar BFR
#' @param newdata \code{data.frame} Predictions
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
predict.BFR <- function(object,newdata,...){
    if (!inherits(object, "BFR")) stop("This function only works for objects of class 'BFR'")
	object$predictions
}


#' @title plot.BFR
#'
#' @description Solo es una prueba
#'
#' @param x \code{BFR object} Objeto BFR, resultado de ejecutar BFR
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom graphics plot abline
#' @export
plot.BFR <- function(x, ...){
  ### Check that object is compatible
  if (!inherits(x, "BFR")) stop("This function only works for objects of class 'BFR'")

  limits <- range(c(x$response, x$predictions), na.rm = TRUE)
  plot(x$response, x$predictions, main = "Training", xlim = limits, ylim = limits, xlab = 'Response', ylab = 'Prediction', ...);
  abline(a = 0, b = 1, lty = 3)
}

#' @title boxplot.BFRCV
#'
#' @description Solo es una prueba
#'
#' @param x \code{BFRCV object} Objeto BFRCV, resultado de ejecutar BFR() con el parametro folds > 2
#' @param select \code{string} Pearson or MSEP
#' @param ordered \code{logic} TRUE or FALSE
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom graphics boxplot
#' @export
boxplot.BFRCV <- function(x, select = 'Pearson', ordered = TRUE, ...){
  ### Check that object is compatible
  if (!inherits(x, "BFRCV")) stop("This function only works for objects of class 'BFRCV'")

  results <- x$predictions_Summary

  if (select == "Pearson") {
    plot.y <- results$Pearson
    ylab <- "Pearson's Correlation"
  } else if (select == "MSEP") {
    plot.y <- results$MSEP
    ylab <- "MSEP Average"
  }

  if (length(unique(results$Env)) > 1) {
    results$TxE <- paste0(results$Trait, '_', results$Env)

    if (ordered) {
      results$TxE  <- with(results, reorder(TxE , Pearson, median, na.rm = T))
    }
    boxplot(plot.y ~ results$TxE, col = "grey", ylab = ylab, ...)
  }else{
    boxplot(plot.y, col = "grey", xlab = 'Environment', ylab = ylab, ...)
  }


}

#' @title Plot BFRCV graph
#'
#' @description Plot from BFRCV object
#'
#' @param x \code{BFRCV object} BFRCV object, result of use the BFR() function
#' @param select \code{character} By default ('Pearson'), plot the Pearson Correlations of the BFR Object, else ('MSEP'), plot the MSEP of the BFRCV Object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom graphics arrows axis plot
#' @export
plot.BFRCV <- function(x, select = 'Pearson', ...){
  ### Check that object is compatible
  if (!inherits(x, "BFRCV")) stop("This function only works for objects of class 'BFRCV'")

  results <- x$predictions_Summary[which(x$predictions_Summary$Fold == 'Average_all'), ]

  results <- results[order(results[, select]),]

  if (select == "Pearson") {
    results$SE <- 1.96 * results$SE_Pearson
    ylab <- "Pearson's Correlation"
  } else if (select == "MSEP") {
    results$SE <- 1.96 * results$SE_MSEP[which(results$Fold == 'Average_all')]
    ylab <- select
  }
  x.labels <- paste0(results$Trait, '_', results$Env)
  plot.x <- 1:length(x.labels)
  plot(plot.x, results[, select], ylim = range(c(results[, select] - results$SE, results[, select] + results$SE)),
       type = 'p', ylab = ylab, xlab = '', xaxt = "n", las = 2)
  axis(1, at = plot.x, labels = x.labels, las = 2)
  arrows(plot.x, results[, select] - results$SE, plot.x, results[, select] + results$SE, code = 3, length = 0.02, angle = 90)
}

