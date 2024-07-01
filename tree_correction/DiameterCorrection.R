#' Diameter correction
#'
#' @param Data Dataset (data.frame or data.table)
#'   The dataset must contain the columns:
#'   - `IdTree` or `IdStem` (character)
#'   - `ScientificName_TreeDataCor` (character)
#'   - `Diameter` (numeric)
#'   - `Year` (numeric)
#'   - **`POM` (Point Of Measurement) (factor)** or
#'     **`HOM` (Height Of Measurement) (numeric)** if you want to correct from
#'      the **"POM change"**
#'   If you want to apply the **"phylogenetic hierarchical"** correction, the
#'   dataset must also contain the columns:
#'   - `Genus_TreeDataCor` (character)
#'   - `Family_TreeDataCor` (character)
#'
#' @param KeepMeas In case of **multiple diameter measurements** in the same
#'   census year:
#' Possible values: "MaxHOM", "MaxDate" (character).
#'   - "MaxHOM": apply the correction to the measurement taken at the
#'               **highest POM**
#'   - "MaxDate": apply the correction to the **most recent measurement** (same
#'                year but more recent date)
#'
#' @param DefaultHOM Default Height Of Measurement in meter (Default: 1.3 m)
#'   (numeric, 1 value)
#'
#' @param MaxDBH Maximum possible DBH (Diameter at the default HOM) of your
#'   stand in cm (numeric, 1 value)
#'
#' @param PositiveGrowthThreshold in cm/year: a tree
#'   widening by more than this value is considered abnormal (numeric, 1 value)
#'
#' @param NegativeGrowthThreshold in cm/census: the possible
#'   positive measurement error (+n) cannot be corrected until the growth
#'   appears abnormal, but a negative measurement error can be allowed until -n
#'   (a tree does not decrease). Thus the positive measurement error (+n) is
#'   "compensated". (numeric, 1 value)
#'
#' @param Pioneers Scientific names of the pioneer species of the site, as in
#'   the `ScientificName_TreeDataCor` column (characters vector)
#'
#' @param PioneersGrowthThreshold in cm/year: a tree of a pioneer species that
#'   widens by more than this value is considered abnormal (numeric, 1 value)
#'
#' @param WhatToCorrect Possible values: "POM change", "punctual", "shift"
#'   (character)
#'   - "POM change": detect POM change in the column `POM` and correct the
#'                   Diameter values from it.
#'   - "punctual": detect if the error is punctual and correct it by
#'                 interpolation.
#'   - "shift": detect if there is a shift of several Diameter values and
#'              links them to the 1st measurements set.
#'
#' @param CorrectionType Possible values: "individual", "phylogenetic
#'   hierarchical" (character).
#'   - "individual": replace abnormal growth by interpolation from the
#'                   individual values.
#'   - "phylogenetic hierarchical": replace abnormal growth with the average
#'          growth of other trees in the dataset, at the specific, genus, family
#'          or stand level, within a DBH range of x cm (*DBHRange* argument).
#'          If the number of these trees < n (*MinIndividualNbr* argument)
#'          at the specific level, we switch to the genus level etc.
#'
#' @param DBHRange DBH range in cm to take into account to select other trees in
#'   the dataset to apply "phylogenetic hierarchical" correction (Default: 10
#'   cm) (numeric, 1 value)
#' @param MinIndividualNbr Minimum number of individuals to take into account in
#'   "phylogenetic hierarchical" correction (Default: 5) (numeric, 1 value)
#'
#' @param OtherCrit Other criteria to select the individuals used for the
#'   calculation of the mean growth in the "phylogenetic hierarchical"
#'   correction. Give the name of the column(s) for which the individuals must
#'   have the same value as the tree to correct (e.g. c("Plot", "Subplot"))
#'   (character)
#'
#' @param Digits Number of decimal places to be used in the `DBHCor` column
#'   (Default: 1L) (integer)
#'
#' @param DBHCorForDeadTrees (logical) TRUE: return DBHCor also for dead trees.
#'   FALSE: do not return DBHCor for dead trees. In this case it is advisable to
#'   have corrected the tree life status with the *StatusCorrection()* function.
#'
#' @param coef description... (numeric)
#'
#' @param DetectOnly TRUE: Only detect errors, FALSE: detect and correct errors
#'   (Default: FALSE) (logical)
#'
#' @return Fill the *Comment* column with error type informations. If
#'   *DetectOnly* = FALSE, add columns:
#'   - *Diameter_TreeDataCor*: corrected trees diameter at default HOM
#'   - *DiameterCorrectionMeth* = "local linear regression","weighted
#'       mean"/phylogenetic hierarchical("species"/"genus"/"family"/"stand")/
#'       "shift realignment"/"Same value".
#'   - *POM_TreeDataCor* (factor): POM value at which the corrected diameters are proposed.
#'       Corresponds to the 1st POM value at which the stem was measured.
#'   - *HOM_TreeDataCor* (numeric): HOM value at which the corrected diameters
#'       are proposed. Corresponds to the 1st HOM value at which the stem was
#'       measured.
#'
#' @details When there is only 1 `Diameter` value for a tree/stem,
#'   `Diameter_TreeDataCor` takes the original `Diameter` value. If this value
#'   is 0 or > MaxDBH, `Diameter_TreeDataCor` takes NA. Diameters not linked to
#'   an IdTree/IdStem or to a Census Year are not processed.
#'   Punctual error correction only with linear regression and not quadratic,
#'   because punctual errors are corrected from a local regression with the 2
#'   framing values.
#'
#' @importFrom utils capture.output
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' library(data.table)
#' data(TestData)
#'
#' Rslt <- DiameterCorrection(
#'  TestData,
#'   WhatToCorrect = c("POM change", "punctual", "shift"),
#'     CorrectionType = c("linear", "individual"),
#'     MinIndividualNbr = 1, Digits = 2L)
#'
#' DiameterCorrectionPlot(Rslt, OnlyCorrected = TRUE)
#'
DiameterCorrection <- function(
    Data,
    
    KeepMeas = c("MaxHOM", "MaxDate"),
    
    DefaultHOM = 1.3,
    MaxDBH = 500,
    PositiveGrowthThreshold = 5,
    NegativeGrowthThreshold = -2,
    
    Pioneers = NULL,
    PioneersGrowthThreshold = 7.5,
    
    WhatToCorrect = c("POM change", "punctual", "shift"),
    CorrectionType = c("individual", "phylogenetic hierarchical"),
    
    DBHRange = 10,
    MinIndividualNbr = 5,
    OtherCrit = NULL,
    Digits = 1L,
    
    DBHCorForDeadTrees = TRUE,
    
    coef = 0.9,
    
    DetectOnly = FALSE
){
  
  #### Arguments check ####
  # Data
  if (!inherits(Data, c("data.table", "data.frame")))
    stop("Data must be a data.frame or data.table")
  
  # IdStem or IdTree? ---------------------------------------------------------------------------------------
  # If no IdStem take IdTree
  if((!"IdStem" %in% names(Data) | all(is.na(Data$IdStem))) &
     ("IdTree" %in% names(Data) & any(!is.na(Data$IdTree))) ){ ID <- "IdTree"
     
  }else{ ID <- "IdStem"}
  
  if(!any(c("IdStem", "IdTree") %in% names(Data)) | (all(is.na(Data$IdStem)) &  all(is.na(Data$IdTree))) )
    stop("The 'IdStem' or 'IdTree' column is missing in your dataset")
  # ---------------------------------------------------------------------------------------------------------
  
  # Diameter column exists
  if(!"Diameter" %in% names(Data))
    stop("The 'Diameter' column does't exist in the dataset")
  
  # DefaultHOM/Min-MaxDBH/Positive-Negative-PioneersGrowthThreshold/DBHRange/MinIndividualNbr (numeric, 1 value)
  if(!all(unlist(lapply(list(DefaultHOM, MaxDBH,
                             PositiveGrowthThreshold, NegativeGrowthThreshold, PioneersGrowthThreshold,
                             DBHRange, MinIndividualNbr),
                        length)) %in% 1) |
     !all(unlist(lapply(list(PositiveGrowthThreshold, NegativeGrowthThreshold, DefaultHOM, PioneersGrowthThreshold),
                        inherits, c("numeric", "integer")))))
    stop("The 'PositiveGrowthThreshold', 'NegativeGrowthThreshold', 'PioneersGrowthThreshold' and 'DefaultHOM' arguments
         of the 'DiameterCorrection' function must be 1 numeric value each")
  
  # Pioneers (characters vector)
  if(!inherits(Pioneers, "character") & !is.null(Pioneers))
    stop("'Pioneers' argument must be a characters vector, or NULL")
  
  # WhatToCorrect
  if(!any(c("POM change","punctual", "shift") %in% WhatToCorrect))
    stop("The 'WhatToCorrect' argument value must be among 'POM change', 'punctual' and 'shift'")
  
  # CorrectionType
  if(!any(c("linear", "quadratic", "individual", 'phylogenetic hierarchical') %in% CorrectionType))
    stop("The 'CorrectionType' argument value must be among
         'quadratic', 'linear', 'individual' and 'phylogenetic hierarchical'")
  
  # Digits
  if(!inherits(Digits, "integer") & Digits != as.integer(Digits))  {
    warning(paste0("The 'Digits' argument must be an integer. Value entered (", Digits, ")  coerced to ", as.integer(Digits), "."))
    Digits <- as.integer(Digits)
  }
  
  # DetectOnly (logical)
  if(!inherits(DetectOnly, "logical"))
    stop("The 'DetectOnly' argument must be a logical")
  
  # Taper before if 'HOM' in the dataset and not 'TaperDBH_TreeDataCor'
  if(any(!is.na(Data$HOM)) & !"TaperDBH_TreeDataCor" %in% names(Data)) # HOM exists?
    message("You have the 'HOM' information in your dataset.
            We advise you to correct your diameters also with the 'taper' correction (TaperCorrection() function)")
  
  # If 'POM' 'POM change' correction is advised
  if((all(is.na(Data$HOM)) | !"HOM" %in% names(Data)) &
     any(!is.na(Data$POM)) & !any(WhatToCorrect %in% "POM change")) # POM exists?
    message("You have the 'POM' information in your dataset.
            We advise you to correct your diameters also from the 'POM change' ('WhatToCorrect' argument)")
  
  # 'POM change' correction needs 'POM' or 'HOM' values
  if(!any(c("POM", "HOM") %in% names(Data)) | (all(is.na(Data$POM)) &  all(is.na(Data$HOM))) )
    stop("You have chosen to make a 'POM change' correction,
       but you do not have the necessary 'POM' or HOM' column in your dataset or they are empty")
  
  # In data.table
  setDT(Data)
  
  Data <- unique(Data)   # if there are duplicate rows, delete them
  
  # Dataset with the dead trees if no correction wanted for them --------------------------------------------
  if("LifeStatus_TreeDataCor" %in% names(Data)){ Status <- "LifeStatus_TreeDataCor"
  }else if ("LifeStatusCor" %in% names(Data)){ Status <- "LifeStatusCor"
  }else if ("LifeStatus" %in% names(Data)){ Status <- "LifeStatus"
  }else{stop("You have chosen DBHCorForDeadTrees = FALSE.
             To apply this choice the dataset must contain the column
             'LifeStatus_TreeDataCor', 'LifeStatusCor' or 'LifeStatus'")}
  
  if(DBHCorForDeadTrees == FALSE){
    DeadTrees <- Data[get(Status) == FALSE]
    Data <- Data[get(Status) == TRUE | is.na(get(Status))] # AliveTrees
    # nrow(Data) == nrow(AliveTrees) + nrow(DeadTrees) # to check
  }
  
  # Remove duplicated measurements per Year because different POM or Date -----------------------------------
  CompleteData <- copy(Data)
  Data <- UniqueMeasurement(Data, KeepMeas = KeepMeas, ID = ID)
  
  DuplicatedRows <- CompleteData[!Data, on = .NATURAL] # rows removed
  
  # Remove duplicated measurements (randomly)
  # Data <- Data[!duplicated(Data[, list((get(ID), Year)], fromLast = TRUE)] # keep the last measurement
  
  
  if(!"Comment" %in% names(Data)) Data[, Comment := ""]
  if(DetectOnly %in% FALSE){
    if(!"DiameterCorrectionMeth" %in% names(Data)) Data[, DiameterCorrectionMeth := ""]
  }
  
  # If no diameter value, write a comment
  Data <- GenerateComment(Data,
                          condition = is.na(Data[, Diameter]),
                          comment = "Missing value in 'Diameter'")
  
  #### Function ####
  
  # Order IDs and times in ascending order ----------------------------------------------------------------------------
  Data <- Data[, Year := as.numeric(Year)]
  Data <- Data[order(get(ID), Year)]
  
  # IDs vector --------------------------------------------------------------------------------------------------------
  Ids <- as.vector(na.omit(unique(Data[!is.na(Year), get(ID)]))) # Tree Ids
  
  # Dataset with the rows without IDS ----------------------------------------------------------------------------------
  DataIDNa <- Data[is.na(get(ID))]
  
  # Dataset with the rows without Year ----------------------------------------------------------------------------------
  DataYearNa <- Data[is.na(Year)]
  
  # Taper correction ------------------------------------------------------------------------------------------------------
  # if("taper" %in% CorrectionType) {
  #   Data <- TaperCorrection(Data,
  #                           DefaultHOM = DefaultHOM,
  #                           TaperParameter = TaperParameter, TaperFormula = TaperFormula,
  #                           DetectOnly = DetectOnly)
  # }
  
  
  
  # Apply for all the trees -----------------------------------------------------------------------------------------------
  # i = "100635"
  Data <- do.call(rbind, lapply(Ids, function(i) DiameterCorrectionByTree(
    DataTree = Data[get(ID) %in% i & !is.na(Year)], # per ID, all censuses
    Data = Data,
    
    DefaultHOM = DefaultHOM,
    MaxDBH = MaxDBH,
    PositiveGrowthThreshold = PositiveGrowthThreshold,
    NegativeGrowthThreshold = NegativeGrowthThreshold,
    
    Pioneers = Pioneers,
    PioneersGrowthThreshold = PioneersGrowthThreshold,
    
    WhatToCorrect = WhatToCorrect,
    CorrectionType = CorrectionType,
    
    DBHRange = DBHRange,
    MinIndividualNbr = MinIndividualNbr,
    OtherCrit =  OtherCrit,
    
    Digits = Digits,
    
    coef = coef,
    
    DetectOnly = DetectOnly
  )
  )) # do.call apply the 'rbind' to the lapply result
  
  # Re-put the rows duplicated, or without ID or Year -----------------------------------------------------------------
  Data <- rbindlist(list(Data, DuplicatedRows, DataIDNa, DataYearNa), use.names = TRUE, fill = TRUE)
  
  # Re-put the dead trees in the dataset (if there are not corrected by choice)
  if(DBHCorForDeadTrees == FALSE){
    Data <- rbindlist(list(Data, DeadTrees), use.names = TRUE, fill = TRUE)
  }
  
  # Order IDs and times in ascending order ----------------------------------------------------------------------------
  Data <- Data[order(get(ID), Year)]
  
  if(DetectOnly %in% FALSE){
    # Rename correction columns
    setnames(Data, c("DBHCor", "POMCor", "HOMCor"),
             c("Diameter_TreeDataCor", "POM_TreeDataCor", "HOM_TreeDataCor"),
             skip_absent = TRUE)
  }
  
  
  return(Data)
  
}

#' DiameterCorrectionByTree
#'
#' @param DataTree A dataset corresponding to a single tree/stem's (1 IdTree/IdStem)
#'   measurements (data.frame or data.table).
#'
#' @param Data Complete dataset (data.table) used if the "phylogenetic
#'   hierarchical" correction (*CorrectionType* argument) is chosen.
#'   The dataset must contain the columns:
#'   - `IdStem` or `IdTree` (character)
#'   - `ScientificNameCor` (character)
#'   - `GenusCor` (character)
#'   - `FamilyCor` (character)
#'   - `Diameter` (numeric)
#'   - `Year` (numeric)
#'
#' @param DefaultHOM Default Height Of Measurement in meter (Default: 1.3 m)
#'   (numeric, 1 value)
#'
#' @param MaxDBH Maximum possible DBH of your stand in cm (numeric, 1 value)
#'
#' @param PositiveGrowthThreshold in cm/year: a tree
#'   widening by more than x cm/year is considered abnormal (numeric, 1 value)
#'
#' @param NegativeGrowthThreshold in cm/census: the possible
#'   positive measurement error (+n) cannot be corrected until the growth
#'   appears abnormal, but a negative measurement error can be allowed until -n
#'   (a tree does not decrease). Thus the positive measurement error (+n) is
#'   "compensated". (numeric, 1 value)
#'
#' @param Pioneers Scientific names of the pioneer species of the site, as in
#'   the 'ScientificNameCor' column (characters vector)
#'
#' @param PioneersGrowthThreshold in cm/year: a tree of a pioneer species that
#'   widens by more than x cm/year is considered abnormal (numeric, 1 value)
#'
#' @param WhatToCorrect  c("POM change", "punctual", "shift") (character)
#'   - "POM change": detect POM change in the column 'POM' and correct
#'                   the Diameter values from it.
#'   - "punctual": detect if the error is punctual and correct it by
#'                 interpolation.
#'   - "shift": detect if there is a shift of several 'Diameter' values and
#'              links them to the trust measurements set.
#'
#' @param CorrectionType Possible values: "individual", "phylogenetic
#'   hierarchical" (character).
#'   - "individual": replace abnormal growth by interpolation from the
#'                   individual values.
#'   - "phylogenetic hierarchical": replace abnormal growth with the average
#'          growth of other trees in the dataset, at the specific, genus, family
#'          or stand level, within a DBH range of x cm (*DBHRange* argument).
#'          If the number of these trees < n (*MinIndividualNbr* argument)
#'          at the specific level, we switch to the genus level etc.
#'
#' @param DBHRange DBH range in cm to take into account to select other trees in
#'   the dataset to apply "phylogenetic hierarchical" correction (Default: 10
#'   cm) (numeric, 1 value)
#' @param MinIndividualNbr Minimum number of individuals to take into account in
#'   "phylogenetic hierarchical" correction (Default: 5) (numeric, 1 value)
#'
#' @param OtherCrit Other criteria to select the individuals used for the
#'   calculation of the mean growth in the "phylogenetic hierarchical"
#'   correction. Give the name of the column(s) for which the individuals must
#'   have the same value as the tree to correct (e.g. c("Plot", "Subplot"))
#'   (character)
#'
#' @param Digits Number of decimal places to be used in the 'DBHCor' column
#'   (Default: 1L) (integer)
#'
#' @param coef description... (numeric)
#'
#' @param DetectOnly TRUE: Only detect errors, FALSE: detect and correct errors
#'   (Default: FALSE) (logical)
#'
#' @return Fill the *Comment* column with error type informations. If
#'   *DetectOnly* = FALSE, add columns:
#'   - *DBHCor*: corrected trees diameter at default HOM
#'   - *DiameterCorrectionMeth* = "local linear regression"/"weighted
#'   mean"/phylogenetic hierarchical("species"/"genus"/"family"/"stand")/"shift
#'   realignment"/"Same value".
#'
#' @details When there is only 1 `Diameter` value for a tree/stem, `DBHCor`
#'   takes the original `Diameter` value. If this value is 0 or > MaxDBH,
#'   `DBHCor` takes NA.
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' library(data.table)
#' data(TestData)
#'
#'  DataTree <- data.table(IdStem = "c",
#'       ScientificName = "A",
#'       Year = c(seq(2000,2008, by = 2), 2012, 2014,2016, 2020), # 9 Diameter values
#'       Diameter = c(13:16, 16-4, (16-4)+2, (16-4)+3, 15-4, (15-4)+2), # 0.5 cm/year
#'       POM = as.factor(c(0, 0, 0, 0, 1, 1, 1, 2, 2)))
#'
#' Rslt <- DiameterCorrectionByTree(
#'   DataTree, TestData,
#'   WhatToCorrect = c("POM change", "punctual", "shift"),
#'   CorrectionType = c("linear", "individual")
#'   )
#' setnames(Rslt, "POMCor", "POM_TreeDataCor")
#' DiameterCorrectionPlot(Rslt, CorCol = "DBHCor")
#'
DiameterCorrectionByTree <- function(
    DataTree,
    Data,
    
    DefaultHOM = 1.3,
    MaxDBH = 500,
    PositiveGrowthThreshold = 5,
    NegativeGrowthThreshold = -2,
    
    Pioneers = NULL,
    PioneersGrowthThreshold = 7.5,
    
    WhatToCorrect = c("POM change", "punctual", "shift"),
    CorrectionType = "individual",
    
    DBHRange = 10,
    MinIndividualNbr = 5,
    OtherCrit = NULL,
    Digits = 1L,
    
    coef = 0.9,
    
    DetectOnly = FALSE
){
  
  # print(unique(DataTree[, get(ID)])) # to debug
  
  
  #### Arguments check ####
  # DataTree
  if (!inherits(DataTree, c("data.table", "data.frame")))
    stop("DataTree must be a data.frame or data.table")
  
  # IdStem or IdTree? ---------------------------------------------------------------------------------------
  # If no IdStem take IdTree
  if((!"IdStem" %in% names(DataTree) | all(is.na(DataTree$IdStem))) &
     ("IdTree" %in% names(DataTree) & any(!is.na(DataTree$IdTree))) ){ ID <- "IdTree"
     
  }else{ ID <- "IdStem"}
  
  if(!any(c("IdStem", "IdTree") %in% names(DataTree)) | (all(is.na(DataTree$IdStem)) &  all(is.na(DataTree$IdTree))) )
    stop("The 'IdStem' or 'IdTree' column is missing in your dataset")
  # ---------------------------------------------------------------------------------------------------------
  
  
  # if there are several IdTrees
  if(length(unique(DataTree[, get(ID)])) != 1){
    stop("DataTree must correspond to only 1 same tree/stem so 1 same ",ID,"
    (the ",ID,"s: " ,paste0(unique(DataTree[, get(ID)]), collapse = "/"),"")
  }
  
  # In data.table
  setDT(DataTree)
  
  
  # Arrange year in ascending order
  DataTree <- DataTree[order(Year)] # data.table::order
  
  if(!"TaperDBH_TreeDataCor" %in% names(DataTree)){
    if(DetectOnly %in% FALSE){
      if("POM" %in% names(DataTree)) DataTree[, POMCor := POM[!is.na(POM)][1]] # Corrected diameter is at the 1st POM
      if("HOM" %in% names(DataTree)) DataTree[, HOMCor := HOM[!is.na(HOM)][1]] # Corrected diameter is at the 1st  HOM
    }
  }
  
  
  # If not enough Diameter values
  if(sum(!is.na(DataTree$Diameter)) > 1){
    
    #### Function ####
    
    # Pioneers species case
    
    if(!is.null(Pioneers) & PioneersGrowthThreshold != PositiveGrowthThreshold){
      
      ## ScientificName_TreeDataCor or ScientificName?
      if("ScientificName_TreeDataCor" %in% names(DataTree)){
        SfcName <- "ScientificName_TreeDataCor"
        
      }else if(!"ScientificName_TreeDataCor" %in% names(DataTree) & "ScientificName" %in% names(DataTree)){
        SfcName <- "ScientificName"
        
      }else if(!any(c("ScientificName_TreeDataCor", "ScientificName") %in% names(DataTree)))
        
        stop("There are no 'ScientificName_TreeDataCor' or 'ScientificName' column.
           It is not possible to take into account the pioneer character of species in the diameter correction.
             If you do not want to take into account the pioneer character in the diameter correction,
             leave the argument Pioneers = NULL.")
      
      if(any(na.omit(unique(DataTree[, get(SfcName)]) == Pioneers))){ # if it's a pioneer species
        
        PositiveGrowthThreshold <- PioneersGrowthThreshold # take the Pioneers growth threshold
      }
      
    } # end Pioneers criteria
    
    # If the taper correction has been made, start from it
    if("TaperDBH_TreeDataCor" %in% names(DataTree)) DBHCor <- Diameter <- DataTree[, TaperDBH_TreeDataCor]
    if(!"TaperDBH_TreeDataCor" %in% names(DataTree)) DBHCor <- Diameter <- DataTree[, Diameter]
    
    Time <- DataTree[, Year]
    
    # DBH = 0 is impossible
    DBHCor[DBHCor == 0] <- NA
    
    # DBH > MaxDBH -> DBH = NA
    DBHCor[DBHCor > MaxDBH] <- NA
    
    # Correction with POM ---------------------------------------------------------------------------------------------------
    if("POM change" %in% WhatToCorrect){
      
      # POM or HOM?
      # If no POM take HOM
      if((!"POM" %in% names(DataTree) | all(is.na(DataTree$POM))) &
         ("HOM" %in% names(DataTree) & any(!is.na(DataTree$HOM))) ){ POMv <- "HOM"
         
      }else{ POMv <- "POM"}
      
      if(!any(c("POM", "HOM") %in% names(DataTree)) | (all(is.na(DataTree$POM)) &  all(is.na(DataTree$HOM))) )
        message("You have chosen to make a 'POM change' correction,
        but 'POM' and HOM' columns are empty for ", unique(DataTree[, get(ID)]),".
             This correction will therefore not be applied to this stem.")
      
      
      ## POM change detection -----------------------------------------------------------------------------------------------
      if(any(!is.na(DataTree[, get(POMv)]))) { # POM exists?
        
        # Check the POM value over the time
        # If POM decreases comment it
        DataTree <- GenerateComment(DataTree,
                                    condition = as.numeric(rownames(DataTree)) %in%
                                      (which(diff(as.numeric(DataTree[, get(POMv)])) < 0) +1),
                                    comment = "POM decrease")
        
        # POM change detection
        POMChange <- NA  # 1st val = NA because it's the default POM
        for( n in (2:(length(DataTree[, get(POMv)]))) ){
          POMChange <- c(POMChange, DataTree[, get(POMv)][n-1] != DataTree[, get(POMv)][n]) # (TRUE = POM change)
        }
        
        raised = which(POMChange)-1 # which are TRUE (-1 to be in cresc indice and not DBH indice)
        
        
        if(length(raised) != 0){ # if there are POM changes
          
          DataTree <- GenerateComment(DataTree,
                                      condition = as.numeric(rownames(DataTree)) %in% (raised+1), # +1 to be in DBH indice
                                      comment = paste0("POM change"))
          
          if(DetectOnly %in% FALSE){
            
            # Compute diameter incrementation without the inits shift
            cresc <- ComputeIncrementation(Var = DBHCor, Type = "annual", Time = Time)
            cresc_abs <- ComputeIncrementation(Var = DBHCor, Type = "absolute", Time = Time)
            # Remove incr between 2 shifts (take growth only intra seq)
            cresc[raised]  <- NA
            cresc_abs[raised] <- NA
            
            # Put NA if other abnormal incrementation
            AbnormalCrescs <- (cresc >= PositiveGrowthThreshold | cresc_abs < NegativeGrowthThreshold)
            cresc[AbnormalCrescs]  <- NA
            cresc_abs[AbnormalCrescs]  <- NA
            
            if(length(cresc[!is.na(cresc)]) > 0){ # if cresc != NA
              
              if("individual" %in% CorrectionType) {
                
                IndCorRslt <- IndividualDiameterShiftCorrection(DataTree = DataTree,
                                                                DBHCor = DBHCor, Time = Time,
                                                                cresc = cresc, cresc_abs = cresc_abs,
                                                                cresc_abn = raised,
                                                                coef = coef)
                
                DataTree <- IndCorRslt$DataTree
                DBHCor <- IndCorRslt$DBHCor
                
                
              } # end individual correction
              
            } # end if cresc != NA
            
            if(!"individual"%in% CorrectionType & "phylogenetic hierarchical" %in% CorrectionType){
              
              DataTree <- PhylogeneticHierarchicalCorrection(
                DataTree = DataTree,
                Data = Data,
                cresc = cresc, cresc_abs = cresc_abs, cresc_abn = raised,
                DBHCor = DBHCor, Time = Time,
                PositiveGrowthThreshold = PositiveGrowthThreshold,
                NegativeGrowthThreshold = NegativeGrowthThreshold,
                DBHRange = DBHRange, MinIndividualNbr = MinIndividualNbr, OtherCrit = OtherCrit, coef = coef)
              
              DBHCor <- DataTree[,DBHCor]
            }
            
            ## 3. + trunk width reduction factor (if POM change (only?)) ------------------------------------------------------
            
          } # End correction "POM change"
          
        }# if there are POM changes
      }# if there are POMs
    }# Correction with POM
    
    
    
    # Punctual/shift error detection  + replace with NA if punctual ---------------------------------------------------------
    if(any("punctual" %in% WhatToCorrect | "shift" %in% WhatToCorrect)){
      
      if(length(DBHCor[!is.na(DBHCor)]) == 2){ # if only 2 non-NA values
        
        TwoValCorRslt <- TwoValDiameterCor(DataTree = DataTree,
                                           Data = Data,
                                           DBHCor = DBHCor, Time = Time,
                                           CorrectionType = CorrectionType,
                                           PositiveGrowthThreshold = PositiveGrowthThreshold,
                                           NegativeGrowthThreshold = NegativeGrowthThreshold,
                                           DBHRange = DBHRange,
                                           MinIndividualNbr = MinIndividualNbr,
                                           OtherCrit = OtherCrit,
                                           DetectOnly = DetectOnly)
        
        DataTree <- TwoValCorRslt$DataTree
        DBHCor <- TwoValCorRslt$DBHCor
        
      }else{
        
        DBHCor <- PunctualErrorDetection(
          DBHCor = DBHCor, Time = Time,
          PositiveGrowthThreshold = PositiveGrowthThreshold, NegativeGrowthThreshold = NegativeGrowthThreshold,
          DetectOnly = DetectOnly)
        # ça serait bien de renvoyer qqchose si un shift est detecté pour être plus secure (y refléchir)
        
        if("DBHCor" %in% names(DataTree)){
          DataTree[, DBHCor := NULL] # remove the DBHCor col to avoid conflict
        }
        
        DataTree[,DBHCor := DBHCor]
        
        DataTree <- GenerateComment(DataTree,
                                    condition = (is.na(DataTree[,DBHCor]) & !is.na(DataTree[,Diameter])),
                                    comment = paste0("Abnormal diameter value (punctual error)"))
        
        if(DetectOnly %in% TRUE) DataTree[,DBHCor := NULL] # remove the DBHCor col if we detect only
      }
    }
    
    # Shift Correction ------------------------------------------------------------------------------------------------------
    if("shift" %in% WhatToCorrect){
      ## Init shift detection si PunctualErrorDetection() ne s'en est pas chargé --------------------------------------------
      ### Compute diameter incrementation without the inits shift
      cresc <- ComputeIncrementation(Var = DBHCor, Type = "annual", Time = Time)
      cresc_abs <- ComputeIncrementation(Var = DBHCor, Type = "absolute", Time = Time)
      
      ### Detect abnormal growth --------------------------------------------------------------------------------------------
      cresc_abn <- which(cresc >= PositiveGrowthThreshold | cresc_abs < NegativeGrowthThreshold) # abnormal values indices
      # le retour à la normale est considéré comme une erreur (perte excessive)
      
      if(length(cresc_abn) != 0) { # if there are abnormal values
        
        if("DBHCor" %in% names(DataTree)){
          DataTree[, DBHCor := NULL] # remove the DBHCor col to avoid conflict
        }
        
        DataTree[,DBHCor := DBHCor]
        
        DataTree <- GenerateComment(DataTree,
                                    condition = as.numeric(rownames(DataTree)) %in% (cresc_abn+1),
                                    comment = paste0("Abnormal diameter value (shift error)"))
        
        if(DetectOnly %in% TRUE) DataTree[,DBHCor := NULL] # remove the DBHCor col if we detect only
        
        
        if(DetectOnly %in% FALSE){
          
          # Remove abnormal growths
          cresc[cresc_abn] <- NA
          cresc_abs[cresc_abn] <- NA
          
          if(length(cresc[!is.na(cresc)]) > 0){ # if cresc != NA
            
            if("individual" %in% CorrectionType) {
              
              IndCorRslt <- IndividualDiameterShiftCorrection(DataTree = DataTree,
                                                              DBHCor = DBHCor, Time = Time,
                                                              cresc = cresc, cresc_abs = cresc_abs,
                                                              cresc_abn = cresc_abn,
                                                              coef = coef)
              
              DataTree <- IndCorRslt$DataTree
              DBHCor <- IndCorRslt$DBHCor
              
              
            } # end individual correction
            
          } # end if cresc != NA
          
          if(!"individual"%in% CorrectionType & "phylogenetic hierarchical" %in% CorrectionType){
            DataTree <- PhylogeneticHierarchicalCorrection(
              DataTree = DataTree,
              Data = Data,
              cresc = cresc, cresc_abs = cresc_abs, cresc_abn = cresc_abn,
              DBHCor = DBHCor, Time = Time,
              PositiveGrowthThreshold = PositiveGrowthThreshold,
              NegativeGrowthThreshold = NegativeGrowthThreshold,
              DBHRange = DBHRange, MinIndividualNbr = MinIndividualNbr, OtherCrit = OtherCrit, coef = coef)
            
            DBHCor <- DataTree[,DBHCor]
          }
          
          ## 3. + trunk width reduction factor (if POM change (only?)) ----------------------------------------------------------
        } # End shift correction
      }
    }
    
    
    if(DetectOnly %in% FALSE & "punctual" %in% WhatToCorrect & any(is.na(DBHCor))){ # Na to be replaced
      
      # Compute diameter incrementation without the abnormal values
      cresc <- ComputeIncrementation(Var = DBHCor, Type = "annual", Time = Time)
      cresc_abs <- ComputeIncrementation(Var = DBHCor, Type = "absolute", Time = Time)
      
      # Put NA if other abnormal incrementation
      AbnormalCrescs <- (cresc >= PositiveGrowthThreshold | cresc_abs < NegativeGrowthThreshold)
      AbnormalCrescs <- which(AbnormalCrescs)
      cresc[AbnormalCrescs]  <- NA
      cresc_abs[AbnormalCrescs]  <- NA
      DBHCor[AbnormalCrescs +1] <- NA
      
      i <- which(is.na(DBHCor)) # id of all the NA to interpolate
      
      # Check that only non-abnormal growths are kept
      if(length(which(cresc[!is.na(cresc)] >= PositiveGrowthThreshold | cresc_abs[!is.na(cresc_abs)] < NegativeGrowthThreshold))==0){
        
        # Replace NA by the correction ------------------------------------------------------------------------------------------
        # Regression only with 2 values around the NA (local)
        DBHCor <- RegressionInterpolation(Y = DBHCor, X = Time, CorrectionType = CorrectionType, Local = TRUE) # Compute the corrected cresc
        
        # Add the column with the correction method  ------------------------------------------------------------------------
        # Punctual error correction only with linear regression and not quadratic,
        # because punctual errors are corrected from a local regression with the 2 framing values.
        
        DataTree <- GenerateComment(DataTree,
                                    condition = as.numeric(rownames(DataTree)) %in% (i),
                                    comment = "local linear regression",
                                    column = "DiameterCorrectionMeth")
        
        
      }else{warning("There are still abnormal growths. Either the selected methods are insufficient
                    or the method needs to be improved")}
      
    }
    
    if(DetectOnly %in% FALSE){
      # Check that there are no more abnormal growths -----------------------------------------------------------------------------
      cresc <- ComputeIncrementation(Var = DBHCor, Type = "annual", Time = Time)
      cresc_abs <- ComputeIncrementation(Var = DBHCor, Type = "absolute", Time = Time)
      
      if(any(na.omit(cresc >= PositiveGrowthThreshold | cresc_abs < NegativeGrowthThreshold))){
        
        warning("There are still abnormal growths for the tree/stem ", unique(DataTree[, get(ID)]),". Either the selected methods are insufficient
                    or the method needs to be improved")
      }
    }
    
    
    # 'DBHCor' vector in DataTree -------------------------------------------------------------------------------------------
    if(DetectOnly %in% FALSE){
      
      if("DBHCor" %in% names(DataTree)){
        DataTree[, DBHCor := NULL] # remove the DBHCor col to avoid conflict
      }
      
      DataTree[, DBHCor := round(DBHCor, digits = Digits)] }
    
  }else if (sum(!is.na(DataTree$Diameter)) < 2 & DetectOnly %in% FALSE){ # if only 1 DBH value
    
    if("TaperDBH_TreeDataCor" %in% names(DataTree)) DataTree[, DBHCor := TaperDBH_TreeDataCor] # keep taper Diameter
    if(!"TaperDBH_TreeDataCor" %in% names(DataTree))  DataTree[, DBHCor := Diameter] # keep original Diameter
    
    
    # DBH = 0 or > MaxDBH is impossible
    DataTree[DBHCor == 0 | DBHCor > MaxDBH, DBHCor := NA]
    
  }
  
  return(DataTree)
}
