#' RegressionInterpolation
#'
#' @param Y Response variable with missing value (NA) (numeric)
#' @param X Explanatory variable (numeric)
#'
#' @param CorrectionType "quadratic" and/or "linear" (character)
#'
#' @param Local If FALSE, compute the regression on all the Y and X values, if
#'   TRUE compute the regression between the 2 surrounding values of the missing
#'   value (NA) (default: FALSE) (logical)
#'
#' @details The variables X and Y must be of the same length.
#' Because local regression (*Local* = TRUE) is done with the 2 framing values,
#' this regression is always linear and not quadratic.
#'
#' @return Y (numeric) with interpolated missing values depending on the form of
#'   model chosen (*CorrectionType*)
#'
#'@importFrom stats lm poly
#' @export
#'
#' @examples
#'
#' DBH = c(34.5, NA, 34.0, 34.6, 35.0, 34.9, NA)
#' Time = c(1998, 2008, 2016, 2017, 2018, 2019, 2021)
#' plot(Time, DBH)
#'
#'
#' # Compute the corrected DBH
#' DBHCor <- RegressionInterpolation(Y = DBH, X = Time, Local = FALSE)
#'
#' DBHCorLocal <- RegressionInterpolation(Y = DBH, X = Time, Local = TRUE)
#'
#' plot(Time, DBHCor)
#' plot(Time, DBHCorLocal)
#'
RegressionInterpolation <- function(
    Y,
    X,
    CorrectionType = "linear",
    Local = FALSE
){
  
  #### Arguments check ####
  
  # X/Y (numeric)
  if(!all(unlist(lapply(list(X, Y),
                        inherits, c("numeric", "integer"))))){
    stop("The 'X' and 'Y' variables of the 'RegressionInterpolation' function must be numeric")
  }
  # X & Y of the same lenght
  if(length(X) > length(Y)) Y[(length(Y)+1):length(X)] <- NA
  # stop("The variables X and Y must be of the same length ('RegressionInterpolation()' function)")
  
  if(Local == FALSE){
    # CorrectionType (character)
    if(!any(any(CorrectionType %in% "quadratic") || any(CorrectionType %in% "linear")))
      stop("The 'CorrectionType' argument value must be 'quadratic' and/or 'linear'")
  }
  
  #### Function ####
  
  miss <- which(is.na(Y)) # DBH = NA -> value to replace
  
  # miss: nb of the missing value(s) (which values) (vector)
  
  Y[miss] <- sapply(miss, function(i) { # i = each value de miss
    
    if(Local == TRUE){
      
      yxval <- Select2FramingValues(i = i, Y = Y, X = X)
      yval <- yxval$yval # yval: the 2 Y values around i
      xval <- yxval$xval # xval: the 2 X values corresponding to the two yval values
      
    }else{
      yval = Y # all Y values
      xval = X # all X values
    }
    
    if(length(which(!is.na(Y))) < 2){
      
      yi <- Y[!is.na(Y)]
      
    }else if("quadratic" %in% CorrectionType & length(which(!is.na(yval))) > 2){
      
      # Degree 2 polynomial regression (= quadratic)
      reg <- lm(yval ~ poly(xval, degree = 2, raw = TRUE))$coef # 'degree' must be less than number of unique points
      yi <- reg[1] + reg[2] * X[i] + reg[3] * X[i]^2 # (y = b + ax + cx^2),  DBHi = b + a YEARi + c YEARi^2
      
      # plot(X, Y)
      # curve(reg[1] + reg[2] * x + reg[3] *x^2, add=T)
      
    }else{ # "linear"
      
      # Linear regression: DBH ~ years
      reg <- lm(yval ~ xval)$coef # extract the coefs
      yi <- reg[1] + reg[2] * X[i] # (y = b + ax),  DBHi = b + a YEARi
      
      # plot(X, Y)
      # curve(reg[1] + reg[2] * x, add=T)
      
    }
    
    return(yi) # DBH of i, yi -> Y[miss]
    
  }) # sapply end (for each i)
  
  return(as.numeric(unlist(Y))) # corrected DBHs
  
}