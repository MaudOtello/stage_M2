#' UniqueMeasurement
#'
#' @description Of the repeated measurements of an individual in the same census
#'   year, keep only the measurements (rows)
#'   taken at the **highest POM** (*KeepMeas = "MaxHOM"*),
#'   and/or the **most recent measurement** (same year but more recent date)
#'   (*KeepMeas = "MaxDate"*)
#'
#' @param Data Dataset (data.table)
#'
#' @param KeepMeas In case of **multiple measurements** in the same census year:
#' Possible values: "MaxHOM", "MaxDate" (character).
#'   - "MaxHOM": keep the measurement taken at the **highest HOM/POM**
#'   - "MaxDate": keep the **most recent measurement** (same
#'                year but more recent date)
#'
#' @param ID Column name indicating the identifier of the individual (character)
#'
#' @return Dataset (data.table) with 1 measurement (1 row) per Year.
#'
#' @export
#'
#' @examples
#'
#' library(data.table)
#'
#' Data <- data.table( # 7 rows
#' IdStem = c(rep("a", 4), rep("b", 2), "c"),
#' Year = c(rep(2000, 3), 2001, rep(2000, 3)),
#' Date = as.Date(c("2000-01-10", "2000-01-20", "2000-01-30", "2001-01-10",
#'                  rep("2000-01-10", 3))),
#' HOM = c(rep(1, 3), 2, 1, 2, 1)
#' )
#' Data
#'
#' UniqueMeasurement(TestData)
#'
UniqueMeasurement <- function(
    Data,
    KeepMeas = c("MaxHOM", "MaxDate"), # 1st HOM, then Date, or only one
    ID = "IdStem"
){
  
  Data <- unique(Data)   # if there are duplicate rows, delete them
  
  Rslt <- data.table()
  
  # Separate duplicated and non duplicated --------------------------------------------------------------------------------
  
  DuplicatedID <- Data[duplicated(Data[, list(get(ID), Year)]), list(get(ID), Year)] # duplicated ID-Years
  
  if(nrow(DuplicatedID) > 0){
    
    DuplicatedID[, IDYear := paste(V1, Year, sep = "/")] # code to detect
    
    Data[, IDYear := paste(get(ID), Year, sep = "/")] # code to detect
    
    DuplicatedRows <- Data[IDYear %in% DuplicatedID[, IDYear]] # rows with duplicated ID-Years
    
    UniqueData <- Data[!DuplicatedRows, on = .NATURAL] # rows that do not have duplicate measurements
    
    # nrow(Data) == nrow(DuplicatedRows) + nrow(UniqueData) # to check
    
    
    # Keep the row with the upper POM -------------------------------------------------------------------------------------
    
    if("MaxHOM" %in% KeepMeas){
      
      
      # POM or HOM? ---------------------------------------------------------------------------------------------------
      # If no POM take HOM
      if((!"POM" %in% names(Data) | all(is.na(Data$POM))) &
         ("HOM" %in% names(Data) & any(!is.na(Data$HOM))) ){ POMv <- "HOM"
         
      }else{ POMv <- "POM"}
      
      if(!any(c("POM", "HOM") %in% names(Data)) | (all(is.na(Data$POM)) &  all(is.na(Data$HOM))) )
        message("You have chosen to make corrections to the measurements taken at the highest POM,
        but 'POM' and HOM' columns are empty for ", IdStem,".
             The 'MaxHOM' criterion can therefore not be taken into account.")
      
      
      DuplicatedRows[, MaxHOM := max(get(POMv)), by = list(get(ID), Year)] # compute the max HOM/POM
      NowUniqueData <- DuplicatedRows[get(POMv) == MaxHOM] # keep only this measurement
      NowUniqueData[, MaxHOM := NULL] # remove obsolete col
      
      if(nrow(rbind(UniqueData, NowUniqueData)) != nrow(unique(Data, by = c(ID, "Year"))) ){ # if there are still duplicates at the same POM
        
        if("MaxDate" %in% KeepMeas) DuplicatedRows <- NowUniqueData # to continue with MaxDate
        
        if(!"MaxDate" %in% KeepMeas)
          stop("There are still several measurements per census year despite the choice of the maximum POM.
                You can also select the most recent measurement (same Year but more recent Date) with KeepMeas = 'MaxDate'")
        
      }else{ Rslt <- rbind(UniqueData, NowUniqueData) }
      
    } # end "MaxHOM"
    
    
    # Keep the row with the upper date ------------------------------------------------------------------------------------
    
    if("MaxDate" %in% KeepMeas & nrow(Rslt) == 0){
      
      DuplicatedRows[, MaxDate := max(Date), by = list(get(ID), Year)]
      NowUniqueData <- DuplicatedRows[Date == MaxDate]
      NowUniqueData[, MaxDate := NULL]
      
      if(nrow(rbind(UniqueData, NowUniqueData)) == nrow(unique(Data, by = c(ID, "Year"))) ){ # if there are no more duplicates
        
        Rslt <- rbind(UniqueData, NowUniqueData)
        
      }else{
        
        DuplicatedID <- Data[duplicated(Data[, list(get(ID), Year)]), list(get(ID), Year)]
        
        setnames(DuplicatedID, "V1", paste("Duplicated", ID, sep = ""))
        
        if(nrow(DuplicatedID) > 0){
          
          b <- capture.output(DuplicatedID)
          c <- paste(b, "\n", sep = "")
          # cat("Your data set is:\n", c, "\n")
          
          stop("Duplicated ", ID, "(s) in a census:\n", c, "\n")
          
        }
        
        if(!"MaxHOM" %in% KeepMeas){
          stop("There are still several measurements per census year despite the choice of the maximum POM.
                You can also select the measurement taken at the highest POM with KeepMeas = 'MaxHOM'")
          
        }else{ stop("There are still several measurements per census year
             despite the choice of the maximum POM and the most recent measurement")}
      }
      
    } # end "MaxDate"
    
    Data <- Rslt
    
  } # end if duplicated
  
  return(Data)
  
}