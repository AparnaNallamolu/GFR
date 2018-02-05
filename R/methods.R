#' @title Summary.BGFRA
#'
#' @description Solo es una prueba
#'
#' @param object \code{BGFRA object} Objeto BGFRA, resultado de ejecutar BGFRA
#'
#' @export
summary.BGFRA=function(object,...){

    if(!inherits(object, "BGFRA")) stop("This function only works for objects of class `BGFRA'")

    tmp<-paste('--------------------> Summary of data & model <--------------------')
    cat(tmp,'\n\n')

    tmp<-paste(' Number of phenotypes=', sum(!is.na(object$response)))
    cat(tmp,'\n')

    cat(' Min (TRN)= ', min(object$response,na.rm=TRUE),'\n')
    cat(' Max (TRN)= ', max(object$response,na.rm=TRUE),'\n')
    cat(' Variance of phenotypes (TRN)=', round(var(object$response,na.rm=TRUE),4),'\n')
    cat(' Residual variance=',round(object$varE,4),'\n')

    n<-length(object$response)

    if(any(is.na(object$response)))
    {
     		tst<-which(is.na(object$response))

     		cat(' N-TRN=',n-length(tst), ' N-TST=',length(tst),'\n')

     		cat(' Correlation TRN=',round(cor(object$response[-tst],object$predicted[-tst]),4),'\n')

   }else{
       cat(' N-TRN=',n,'  N-TST=0', '\n\n')
   }

   cat('\n')
   cat(' -- Linear Predictor -- \n')
   cat('\n')
   cat(' Intercept included by default\n')

  for(k in 1:length(object$ETA)){
    if(object$ETA[[k]]$model=="FIXED"){


      if(!is.null(names(object$ETA)[k])){
        cat(" Coefficientes in ETA[",k,"] (",names(object$ETA)[k],") were asigned a flat prior\n")
      }else{
        cat(" Coefficientes in ETA[",k,"] (no-name) are asigned a flat prior\n")
      }


    }else{
      if(object$ETA[[k]]$model=="RKHS"){


        if(!is.null(names(object$ETA)[k])){
          cat(" Coefficientes in ETA[",k,"] (",names(object$ETA)[k],") were assumed to be normally distributed with zero mean and \n covariance (or its eigendecoposition) provided by user \n")
        }else{
          cat(" Coefficientes in ETA[",k,"] (no-name) were assumed to be normally distributed with zero mean and \n covariance (or its eigendecoposition) provided by user \n")
        }


      }else{


        if(!is.null(names(object$ETA)[k])){
				  cat(" Coefficientes in ETA[",k,"] (",names(object$ETA)[k],") modeled as in ", object$ETA[[k]]$model,"\n")
			  }else{
				cat(" Coefficientes in ETA[",k,"] (no-name) modeled as in ", object$ETA[[k]]$model,"\n")
			  }


		  }
	  }
  }

   cat('\n------------------------------------------------------------------\n');
}

#' @title residuals.BGFRA
#'
#' @description Solo es una prueba
#'
#' @param object \code{BGFRA object} Objeto BGFRA, resultado de ejecutar BGFRA
#'
#'
#' @export
residuals.BGFRA <- function(object,...) {
    if(!inherits(object, "BGFRA")) stop("This function only works for objects of class `BGFRA'")
	object$response-object$predicted
}

#' @title predict.BGFRA
#'
#' @description Solo es una prueba
#'
#' @param object \code{BGFRA object} Objeto BGFRA, resultado de ejecutar BGFRA
#'
#' @export
predict.BGFRA <- function(object,newdata,...){
    if (!inherits(object, "BGFRA")) stop("This function only works for objects of class `BGFRA'")
	object$predicted
}


#' @title plot.BGFRA
#'
#' @description Solo es una prueba
#'
#' @param object \code{BGFRA object} Objeto BGFRA, resultado de ejecutar BGFRA
#'
#' @export
plot.BGFRA <- function(x, ...){
  ### Check that object is compatible
  if(!inherits(x, "BGFRA")) stop("This function only works for objects of class `BGFRA'")

  limits<-range(c(x$response,x$predicted),na.rm=TRUE)
  plot(x$response,x$predicted,main="Training",xlim=limits,ylim=limits,xlab='Response',ylab='Prediction');
  abline(a=0,b=1,lty=3)
}

#' @title boxplot.BGFRA-CV
#'
#' @description Solo es una prueba
#'
#' @param object \code{BGFRA-CV object} Objeto BGFRA-CV, resultado de ejecutar BGFRA() con el parametro folds > 2
#'
#' @export
boxplot.BGFRACV <- function(x, ordered=F,...){
  ### Check that object is compatible
  if(!inherits(x, "BGFRACV")) stop("This function only works for objects of class `BGFRA'")

  results <- x$results

  if(length(unique(results$Env))>1){
    if(ordered){
      results$Env <- with(results, reorder(Env , Cor, median, na.rm=T))
    }
    boxplot(results$Cor ~ results$Env, col="grey", xlab='Environment', ylab='Correlation')
  }else{
    boxplot(results$Cor, col="grey", xlab='Environment', ylab='Correlation')
  }
}


