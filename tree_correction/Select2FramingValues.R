#' Select the 2 framing values around a given value
#'
#' @param i Position of the value from which to find the framing values
#'
#' @param Y Vector of values of the variable from which to select the 2 values
#'   around i (numeric)
#' @param X Vector of values of the explanatory variable of the variable Y (thus
#'   also to filter) (numeric)
#'
#' @return 2 vectors of 2 values each (list):
#'  - yval: the 2 Y values around i
#'  - xval: the 2 X values corresponding to the two yval values
#'
#' @export
#'
#' @examples
#' DBHCor = c(17.7, 26.2, NA, 34.6, 34.6, 34.6, 38.0)
#' Time = c(1998, 2008, 2016, 2017, 2018, 2019, 2021)
#'
#' Select2FramingValues(i = 3, Y = DBHCor, X = Time)
#'
Select2FramingValues <- function(
    i,
    Y,
    X
){
  
  # miss <- which(is.na(Y)) # DBH = NA -> value to replace
  
  # miss: nb of the missing value(s) (which values) (vector)
  
  # for(i in miss) { # i = each NA value (in miss)
  
  # i=3
  
  # Case 1: i is the first value of the serie
  if (i < min(which(!is.na(Y)))) {
    ## Keep the 2 next values
    yval <- Y[which(!is.na(Y))[1:min(2, sum(!is.na(Y)))]] # Keep 2 Y values from 1st non-NA to the 2nd if exists
    xval <- X[which(!is.na(Y))[1:min(2, sum(!is.na(Y)))]] # Keep X values linked to theses 2 Y values
  }
  
  # Case 2: i is the last value of the serie
  else if (i > max(which(!is.na(Y)))) {
    ## Keep the 2 previous values
    yval <- Y[which(!is.na(Y))[(sum(!is.na(Y)) - 1):sum(!is.na(Y))]] # Keep the of the second last non-NA value to the last
    xval <- X[which(!is.na(Y))[(sum(!is.na(Y)) - 1):sum(!is.na(Y))]] # X values linked
    
    # yval <- yval[!is.na(yval)] # s'assurer que ce sont des non-NA
    # xval <- xval[!is.na(yval)]
  }
  
  # Case 3: i is in the middle of the serie
  else {
    ## Keep the 2 values around
    yval <-
      Y[c(max( which( !is.na( Y[1:(i - 1)] ) ) ), # take the last value non-NA before NA
          i + min( which( !is.na( Y[(i + 1):length(Y)] ) ) ))] # and the next non-NA value after the NA
    # prendre la valeur la plus récente non-NA avant le NA, et la prochaine valeur non-NA après le NA
    # = prendre les valeurs encadrantes non-NA les plus proches du NA
    
    xval <-
      X[c(max( which( !is.na( Y[1:(i - 1)] ) ) ),
          i + min( which( !is.na( Y[(i + 1):length(Y)] ) ) ))] # X values linked
  }
  
  # } # loop end (for each i)
  
  return(list(yval = yval, xval = xval))
  
}