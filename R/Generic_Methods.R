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

    if (!inherits(object, "BFR")) Error("This function only works for objects of class 'BFR'")

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


#' @title Summary.MTMECV
#'
#' @description Solo es una prueba
#'
#' @param object \code{MTMECV object} Objeto MTMECV, resultado de ejecutar MTME()
#' @param ... Further arguments passed to or from other methods.
#' @param information compact, extended, complete
#'
#' @importFrom stats cor
#' @importFrom dplyr summarise group_by select '%>%' mutate_if funs
#'
#' @export
summary.MTMECV <- function(object, information = 'compact', digits = 4, ...){
  if (!inherits(object, "MTMECV")) Stop("This function only works for objects of class 'MTMECV'")

  object$results %>%
    group_by(Environment, Trait, Partition) %>%
    summarise(Pearson = cor(Predicted, Observed, use = 'pairwise.complete.obs'),
              MSEP = mean((Predicted - Observed)^2, na.rm = T)) %>%
    select(Environment, Trait, Partition, Pearson, MSEP) %>%
    mutate_if(is.numeric, funs(round(., digits))) %>%
    as.data.frame() -> presum

  presum %>%  group_by(Environment, Trait) %>%
    summarise(SE_MSEP = sd(MSEP, na.rm = T)/sqrt(n()), MSEP = mean(MSEP, na.rm = T),
              SE_Pearson = sd(Pearson, na.rm = T)/sqrt(n()), Pearson = mean(Pearson, na.rm = T))  %>%
    select(Environment, Trait, Pearson, SE_Pearson, MSEP, SE_MSEP) %>%
    mutate_if(is.numeric, funs(round(., digits))) %>%
    as.data.frame() -> finalSum

  out <- switch(information,
                compact = finalSum,
                complete = presum,
                extended = {
                  finalSum$Partition <- 'All'
                  presum$Partition <- as.character(presum$Partition)
                  presum$SE_Pearson <- NA
                  presum$SE_MSEP <- NA
                  rbind(presum, finalSum)
                }
  )
  return(out)
}


#' @title Summary.BFRCV
#'
#' @description Solo es una prueba
#'
#' @param object \code{BFRCV object} Objeto BFRCV, resultado de ejecutar BFRCV()
#' @param ... Further arguments passed to or from other methods.
#' @param information compact, extended, complete
#'
#' @importFrom stats cor
#' @importFrom dplyr summarise group_by select '%>%' mutate_if funs
#'
#' @export
summary.BFRCV <- function(object, information = 'compact', digits = 4, ...){
  if (!inherits(object, "BFRCV")) Stop("This function only works for objects of class 'BFRCV'")

  switch (object$response_type,
          gaussian = {
            object$predictions_Summary %>%
              group_by(Environment, Trait, Partition) %>%
              summarise(Pearson = cor(Predicted, Observed, use = 'pairwise.complete.obs'),
                        MSEP = mean((Predicted - Observed)^2, na.rm = T)) %>%
              select(Environment, Trait, Partition, Pearson, MSEP) %>%
              mutate_if(is.numeric, funs(round(., digits))) %>%
              as.data.frame() -> presum

            presum %>%  group_by(Environment, Trait) %>%
              summarise(SE_MSEP = sd(MSEP, na.rm = T)/sqrt(n()), MSEP = mean(MSEP, na.rm = T),
                        SE_Pearson = sd(Pearson, na.rm = T)/sqrt(n()), Pearson = mean(Pearson, na.rm = T))  %>%
              select(Environment, Trait, Pearson, SE_Pearson, MSEP, SE_MSEP) %>%
              mutate_if(is.numeric, funs(round(., digits))) %>%
              as.data.frame() -> finalSum

            if(information == 'extended') {
              finalSum$Partition <- 'All'
              presum$Partition <- as.character(presum$Partition)
              presum$SE_Pearson <- NA
              presum$SE_MSEP <- NA
              extended <- rbind(presum, finalSum)
            }
          },
          ordinal = {
            object$predictions_Summary %>%
              group_by(Environment, Trait, Partition) %>%
              summarise(CC = sum(diag(prop.table(table(factor(Predicted, levels = sort(unique(Observed))), as.factor(Observed))))),
                        MSEP = mean((Predicted - Observed)**2, na.rm = T)) %>%
              select(Environment, Trait, Partition, CC, MSEP) %>%
              as.data.frame() -> presum

            presum %>%
              group_by(Environment, Trait) %>%
              summarise(SE_MSEP = sd(MSEP, na.rm = T), MSEP = mean(MSEP, na.rm = T),
                        SE_CC_2 = sd(CC, na.rm = T)/sqrt(n()),
                        CC = mean(CC, na.rm = T),
                        SE_CC = sqrt((CC*(1 - CC))/n())) %>%
              select(Environment, Trait, CC, SE_CC, SE_CC_2, MSEP, SE_MSEP) %>%
              as.data.frame() -> finalSum

            if(information == 'extended') {
              finalSum$Partition <- 'All'
              presum$Partition <- as.character(presum$Partition)
              presum$SE_CC <- NA
              presum$SE_CC_2 <- NA
              presum$SE_MSEP <- NA
              extended <- rbind(presum, finalSum)
            }
          },
          Error('The $Response type of the dataset was not detected correctly, please check it.')
  )

  out <- switch(information,
                compact = finalSum,
                complete = presum,
                extended = extended
  )
  return(out)
}


#' Print BFRCV information object
#'
#' @param x object a
#' @param ...  more objects
#'
#' @return test
#' @export
#'
print.BFRCV <- function(x, ...){
  cat('Fitted Model with: \n',
      'eta = \n', 'usign ',
      x$nIter, ' iterations, burning the first ', x$burnIn, ' and thining every ', x$thin, '\n',
      'Average runtime by CV', mean(x$executionTime) ,'\n\n',
      'Summary of the predictions: \n')

  print.data.frame(summary(x, 'compact', digits = 3), print.gap = 2L, quote = FALSE)

  cat('\n Use str() function to found more datailed information.')
  invisible(x)
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
    if (!inherits(object, "BFR")) Error("This function only works for objects of class 'BFR'")
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
    if (!inherits(object, "BFR")) Error("This function only works for objects of class 'BFR'")
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
  if (!inherits(x, "BFR")) Error("This function only works for objects of class 'BFR'")

  limits <- range(c(x$response, x$predictions), na.rm = TRUE)
  plot(x$response, x$predictions, main = "Training", xlim = limits, ylim = limits, xlab = 'Response', ylab = 'Prediction', ...);
  abline(a = 0, b = 1, lty = 3)
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
  if (!inherits(x, "BFRCV")) Error("This function only works for objects of class 'BFRCV'")

  results <- summary(x, 'complete')

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
  if (!inherits(x, "BFRCV")) Error("This function only works for objects of class 'BFRCV'")

  results <- summary(x, 'complete')

  switch (select,
          Pearson = {
            plot.y <- results$Pearson
            ylab <- "Pearson's Correlation"
          }, MSEP = {
            plot.y <- results$MSEP
            ylab <- "MSEP Average"
          }, CC = {
            plot.y <- results$CC
            ylab <- "Classification correct average"
          },
          Error('Error in select parameter.')
  )

  if (length(unique(results$Env)) > 1) {
    results$TxE <- paste0(results$Trait, '_', results$Env)

    if (ordered && select != 'MSEP') {
      results$TxE  <- with(results, reorder(TxE , Pearson, median, na.rm = T))
    } else if (ordered && select == 'MSEP') {
      results$TxE  <- with(results, reorder(TxE , MSEP, median, na.rm = T))
    }

    boxplot(plot.y ~ results$TxE, col = "grey", ylab = ylab, ...)
  }else{
    boxplot(plot.y, col = "grey", ylab = ylab, ...)
  }
}

#' @title boxplot.MTMECV
#'
#' @description Solo es una prueba
#'
#' @param x \code{MTMECV object} Objeto MTMECV, resultado de ejecutar MTME()
#' @param select \code{string} Pearson or MSEP
#' @param ordered \code{logic} TRUE or FALSE
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom graphics boxplot
#' @export
boxplot.MTMECV <- function(x, select = 'Pearson', ordered = TRUE, ...){
  ### Check that object is compatible
  if (!inherits(x, "MTMECV")) Error("This function only works for objects of class 'MTMECV'")

  results <- summary(x, 'complete')

  switch (select,
          Pearson = {
            plot.y <- results$Pearson
            ylab <- "Pearson's Correlation"
          }, MSEP = {
            plot.y <- results$MSEP
            ylab <- "MSEP Average"
          }, CC = {
            plot.y <- results$CC
            ylab <- "Classification correct average"
          },
          Error('Error in select parameter.')
  )

  if (length(unique(results$Env)) > 1) {
    results$TxE <- paste0(results$Trait, '_', results$Env)

    if (ordered && select != 'MSEP') {
      results$TxE  <- with(results, reorder(TxE , Pearson, median, na.rm = T))
    } else if (ordered && select == 'MSEP') {
      results$TxE  <- with(results, reorder(TxE , MSEP, median, na.rm = T))
    }

    boxplot(plot.y ~ results$TxE, col = "grey", ylab = ylab, ...)
  }else{
    boxplot(plot.y, col = "grey", ylab = ylab, ...)
  }
}


#' @title Plot MTMECV graph
#'
#' @description Plot from MTMECV object
#'
#' @param x \code{MTMECV object} MTMECV object, result of use the MTME() function
#' @param select \code{character} By default ('Pearson'), plot the Pearson Correlations of the MTME Object, else ('MSEP'), plot the MSEP of the MTMECV Object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom graphics arrows axis plot
#' @export
plot.MTMECV <- function(x, select = 'Pearson', ...){
  ### Check that object is compatible
  if (!inherits(x, "MTMECV")) Error("This function only works for objects of class 'MTMECV'")

  results <- summary(x)
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
       type = 'p', ylab = ylab, xlab = '', xaxt = "n", ...)
  axis(1, at = plot.x, labels = x.labels, las = 2)
  arrows(plot.x, results[, select] - results$SE, plot.x, results[, select] + results$SE, code = 3, length = 0.02, angle = 90)
}

