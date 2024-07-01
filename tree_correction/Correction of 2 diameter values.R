#' Correction of 2 diameter values
#'
#' @param DataTree A dataset corresponding to a single tree/stem's (1
#'   IdTree/IdStem) measurements (data.table)
#'   If 'phylogenetic hierarchical' correction is asked, the dataset must
#'   contain the columns:
#'   - `IdStem` or `IdTree`(character)
#'   - `ScientificNameCor` (character)
#'   - `GenusCor` (character)
#'   - `FamilyCor` (character)
#'   - `Diameter` (numeric)
#'   - `Year` (numeric)
#'
#' @param Data Complete dataset (data.table)
#'   The dataset must contain the columns:
#'   - `IdStem` or `IdTree` (character)
#'   - `ScientificNameCor` (character)
#'   - `GenusCor` (character)
#'   - `FamilyCor` (character)
#'   - `Diameter` (numeric)
#'   - `Year` (numeric)
#'
#' @param DBHCor Diameter vector in cm (numeric)
#' @param Time Time vector in years (numeric)
#'
#' @param CorrectionType Possible values:
#' - NULL: The 1st non-NA diameter value is assigned to the 2nd non-NA diameter
#'         value
#' - "phylogenetic hierarchical": replace abnormal growth with the average
#'       growth of other trees in the dataset, at the specific, genus, family or
#'       stand level, within a DBH range of x cm (*DBHRange* argument). If the
#'       number of these trees < n (*MinIndividualNbr* argument) at the specific
#'       level, we switch to the genus level etc.
#'
#' @param PositiveGrowthThreshold in cm/year : a tree
#'   widening by more than x cm/year is considered abnormal (numeric, 1 value)
#'
#' @param NegativeGrowthThreshold in cm/census : The possible
#'   positive measurement error (+n) cannot be corrected until the growth
#'   appears abnormal, but a negative measurement error can be allowed until -n
#'   (a tree does not decrease). Thus the positive measurement error (+n) is
#'   "compensated". (numeric, 1 value)
#'
#' @param DBHRange DBH range in cm to take into account to select other trees in
#'   the dataset to apply "phylogenetic hierarchical" correction (Default: 10
#'   cm) (numeric, 1 value)
#'
#' @param MinIndividualNbr Minimum number of individuals to take into account in
#'   "phylogenetic hierarchical" correction (Default: 5) (numeric, 1 value)
#'
#' @param OtherCrit Other criteria to select the individuals used for the
#'   calculation of the mean growth. Give the name of the column(s) for which the
#'   individuals must have the same value as the tree to correct (e.g. c("Plot",
#'   "Subplot")) (character)
#'
#' @param DetectOnly TRUE: Only detect errors, FALSE: detect and correct errors
#'   (Default: FALSE) (logical)
#'
#' @return List of 2 objects:
#' - DBHCor (numeric vector): corrected by the method
#' - DataTree (data.table): with the *DiameterCorrectionMeth* column filled with
#'    correction method.
#'
#' @export
#'
#' @examples
#' library(data.table)
#' data(TestData)
#'
#' DataTree <- TestData[IdTree %in% "100771"]
#' DataTree <- DataTree[order(Year)] # order de dt
#'
#' DataTree$Diameter <- c(NA, 13, NA, 8, NA)
#'
#' DBHCor <- DataTree$Diameter
#' Time <- DataTree$Year
#'
#' TwoValDiameterCor(DataTree,
#'                   Data = TestData,
#'                   DBHCor, Time,
#'                   CorrectionType = "phylogenetic hierarchical",
#'                   DBHRange = 10,
#'                   MinIndividualNbr = 1, OtherCrit = NULL)
#'
TwoValDiameterCor <- function(
    DataTree,
    Data,
    DBHCor,
    Time,
    
    CorrectionType = NULL,
    
    PositiveGrowthThreshold = 5,
    NegativeGrowthThreshold = -2,
    
    DBHRange,
    MinIndividualNbr,
    OtherCrit,
    
    DetectOnly = FALSE
){
  
  # Compute diameter incrementation -------------------------------------------------------------------------
  cresc <- ComputeIncrementation(Var = DBHCor, Type = "annual", Time = Time)
  cresc_abs <- ComputeIncrementation(Var = DBHCor, Type = "absolute", Time = Time)
  
  # Detect abnormal growth ----------------------------------------------------------------------------------
  cresc_abn <- which(cresc >= PositiveGrowthThreshold | cresc_abs < NegativeGrowthThreshold) # abnormal values indices
  
  # Correction
  if(DetectOnly %in% FALSE){
    if(length(cresc_abn) != 0) {
      
      
      DataTree <- GenerateComment(DataTree,
                                  condition = as.numeric(rownames(DataTree)) %in% (which(!is.na(DBHCor))[2]),
                                  comment = "Abnormal diameter value")
      
      
      #### Punctual correction ####
      
      if(!"phylogenetic hierarchical" %in% CorrectionType){
        
        
        
        DBHCor[!is.na(DBHCor)][2] <- DBHCor[!is.na(DBHCor)][1] # trust the 1st value
        
        DataTree <- GenerateComment(DataTree,
                                    condition = as.numeric(rownames(DataTree)) %in% (which(!is.na(DBHCor))[2]),
                                    comment = "Same value",
                                    column = "DiameterCorrectionMeth")
        
      } # end punctual correction
      
      
      #### Phylogenetic hierarchical correction ####
      
      if("phylogenetic hierarchical" %in% CorrectionType){
        DataTree <- PhylogeneticHierarchicalCorrection(
          DataTree = DataTree,
          Data = Data,
          cresc = cresc, cresc_abs = cresc_abs, cresc_abn = cresc_abn,
          DBHCor = DBHCor, Time = Time,
          PositiveGrowthThreshold = PositiveGrowthThreshold,
          NegativeGrowthThreshold = NegativeGrowthThreshold,
          DBHRange = DBHRange, MinIndividualNbr = MinIndividualNbr, OtherCrit = OtherCrit)
        
        DBHCor <- DataTree[,DBHCor]
        
      } # end phylogenetic hierarchical correction
      
    } # if cresc abnormal
    
  } # end correction
  
  return(list(DBHCor = DBHCor,
              DataTree = DataTree))
  
}