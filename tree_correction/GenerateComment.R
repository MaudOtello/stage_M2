#' GenerateComment
#'
#' @description For all rows concerned by the condition, the function add a
#'   string in the chosen column of a data.table. If there is already a
#'   value, the string is pasted after the "/" separator.
#'
#' @param Data Dataset (data.table)
#' @param condition Vector of logicals
#' @param comment The string to add in the column (character)
#' @param column Column name to fill (character)
#'
#' @return The input data.table with the filled column.
#'
#' @export
#' @examples
#' \dontrun{
#' library(data.table)
#' dt <- data.table(A = c(2, 1, 3), B = c(6, 10, 6))
#'
#' dt[A == 2,
#'    "information" := paste0("A = 2")] # 1st comment
#'
#' GenerateComment(dt, condition = dt[,B] == 6, comment = "B = 6", column = "information")
#'}
GenerateComment <- function(Data, condition, comment, column = "Comment"){
  
  if(!column %in% names(Data)) Data[, column] <- ""
  
  # Apply the function 'CommentByRow' by row
  for (r in 1:nrow(Data[condition,])) {
    Data[condition,][r,] <- CommentByRow(Data[condition,][r,], comment, column)
  }
  
  return(Data) # in data.table
}


#' CommentByRow
#'
#' @description Add a string in a chosen column of a dataset.
#' If there is already a value, the string is pasted after the "/" separator.
#'
#' @param row A 1 row dataset (data.frame)
#' @param comment The string to add in the chosen column (character)
#' @param column Column name to fill (character)
#'
#' @return The input data.frame with the filled chosen column.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dt <- data.frame(A = c(2), B = c(6))
#'
#' dt <- CommentByRow(dt, comment = "A = 2", column = "information") # 1st comment
#'
#' CommentByRow(dt, comment = "B = 6", column = "information")
#'}
CommentByRow <- function(row, comment, column = "Comment"){
  
  setDF(row) # as data.frame, because it is easier to code in R base than in data.table because of the lazy evaluation
  
  if(!column %in% names(row)) row[, column] <- ""
  
  row[, column] <- ifelse(row[, column] == "",
                          comment,
                          paste(row[, column], comment, sep ="/"))
  return(row)
}

# A <- paste0("A")
# B <- c(A, paste0("B"))
#
# # 1st example
# dt <- data.table(A = c(2, 1), B = c(6, 10))
#
# dt[A == 2,
#    "Comment" := paste0("A = 2")] #  1st comment
# dt[B == 6,
#    "Comment" := paste0("B = 6")] # new comment
#
# dt[B == 6,
#    ("Comment") := c(paste0("B = 6"))] # comments
#
# dt[B == 6]$Comment <- paste0(dt[B == 6]$Comment, ",", "B = 6")
#
# # 2nd example
# dt <- data.table(A = c(2, 1, 3), B = c(6, 10, 6))
#
# dt[A == 2,
#    "Comment" := paste0("A = 2")] # 1st comment
#
# dt[B == 6,
#    "Comment" := paste0(dt[B == 6]$Comment, ",", "B = 6")] # add a comment without errase the first