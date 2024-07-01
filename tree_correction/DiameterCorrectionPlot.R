#' Plot diameter correction result
#'
#' @param Data Dataset (data.frame or data.table)
#'   The dataset must contain the columns:
#'   - `IdStem` (character)
#'   - `Year` (numeric)
#'   - `Diameter` (numeric)
#'   - `DBHCor` (numeric)
#'   - `HOM` (Height Of Measurement) (numeric)
#'   - `HOMCor` (Corrected Height Of Measurement) (numeric)
#'
#' @param OnlyCorrected TRUE: plot only corrected stems, FALSE: plot all stems
#'   (logical)
#'
#' @param SeveralWindows TRUE: return each page in a new window (better
#'   visualisation in Rstudio), FALSE: return each page in the same window
#'   (needed to save all the pages) (logical)
#'
#' @param CorCol Corrected column name (character)
#'
#' @return The plots of the initial diameter values and proposed corrections, by
#'   IdStem.
#'
#' @importFrom ggplot2 ggplot geom_point geom_line aes theme_minimal
#'   position_nudge scale_colour_manual labs vars
#' @importFrom ggforce facet_wrap_paginate
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices dev.new
#'
#' @export
#'
#' @examples
#'
#'\dontrun{
#' pdf("DiameterCorrectionPlots_TestData.pdf", width = 25, height = 10)
#' DiameterCorrectionPlot(Rslt, OnlyCorrected = TRUE, SeveralWindows = FALSE)
#' dev.off()
#'}
#'
DiameterCorrectionPlot <- function(
    Data,
    OnlyCorrected = FALSE,
    CorCol = "Diameter_TreeDataCor",
    # InitialCol = "Diameter"
    SeveralWindows = TRUE
){
  
  #### Arguments check ####
  
  # Data ---------------------------------------------------------------------------------------------------------------
  if(!inherits(Data, c("data.table", "data.frame")))
    stop("Data must be a data.frame or data.table")
  
  # IdStem or IdTree? ---------------------------------------------------------------------------------------
  # If no IdStem take IdTree
  if((!"IdStem" %in% names(Data) | all(is.na(Data$IdStem))) &
     ("IdTree" %in% names(Data) & any(!is.na(Data$IdTree))) ){ ID <- "IdTree"
     
  }else{ ID <- "IdStem"}
  
  if(!any(c("IdStem", "IdTree") %in% names(Data)) | (all(is.na(Data$IdStem)) &  all(is.na(Data$IdTree))) )
    stop("The 'IdStem' or 'IdTree' column is missing in your dataset")
  # ---------------------------------------------------------------------------------------------------------
  
  
  # POM or HOM? ----------------------------------------------------------------------------------------------------------
  # If no HOM take POM
  if((!"HOM" %in% names(Data) | all(is.na(Data$HOM))) &
     ("POM" %in% names(Data) & any(!is.na(Data$POM))) ){ POMv <- "POM"
     
  }else{ POMv <- "HOM"}
  
  if((!"HOM_TreeDataCor" %in% names(Data) | all(is.na(Data$HOM_TreeDataCor))) &
     ("POM_TreeDataCor" %in% names(Data) & any(!is.na(Data$POM_TreeDataCor))) ){ POMcorv <- "POM_TreeDataCor"
     
  }else{ POMcorv <- "HOM_TreeDataCor"}
  
  # Columns --------------------------------------------------------------------------------------------------------------
  # IdStem, Year, Diameter, DBHCor, HOM, HOMCor
  if(!all(c("Year", "Diameter", CorCol, POMv, POMcorv) %in% names(Data)))
    stop(paste0("'Year', 'Diameter', '",CorCol,"', '",POMv,"', ",POMcorv,"' should be columns of Data"))
  
  
  #### Function ####
  
  # Order IDs and times in ascending order ----------------------------------------------------------------------------
  Data <- Data[order(get(ID), Year)]
  
  if(OnlyCorrected == TRUE){
    # Only corrected stems ----------------------------------------------------------------------------------------------
    IDCor <- Data[Diameter != get(CorCol), get(ID)] #  corrected stems
    
    DataCor <- Data[get(ID) %in% IDCor] #  corrected stems
    
  }else{
    DataCor <- Data
    IDCor <- Data[, get(ID)]
  }
  
  
  # Define nrow and ncol for the facet
  n <- length(unique(IDCor))
  if(n<3) { i = 1
  }else{ i = 3}
  
  # Plot --------------------------------------------------------------------------------------------------------------
  
  if(SeveralWindows == TRUE)
    dev.new()
  
  for(p in seq_len(ceiling(length(unique(IDCor))/9))){
    print(ggplot(DataCor) +
            aes(x = Year) +
            
            # Duplicated measurement
            {if(nrow(subset(DataCor, !is.na(Diameter) & is.na(get(CorCol)))) > 0)
              geom_point(data = subset(DataCor, !is.na(Diameter) & is.na(get(CorCol))),
                         aes(y = Diameter,
                             color = 'Duplicated measurement'),
                         shape = "circle", size = 3.9) } +
            # Initial
            geom_point(data = subset(DataCor, !is.na(Diameter)),
                       aes(y = Diameter,
                           color = ifelse(Diameter != get(CorCol), 'Initial', 'Conserved')),
                       shape = "circle", size = 3.9) +
            geom_line(data = subset(DataCor, !is.na(Diameter)),
                      aes(y = Diameter, color = ifelse(Diameter != get(CorCol), 'Initial', 'Conserved'))) +
            
            ggrepel::geom_text_repel(data = subset(DataCor, (!is.na(Diameter) & !is.na(get(POMv)))),
                                     aes(y = Diameter, label = get(POMv), colour = "HOM"),
                                     point.size = 3.9, size = 3, direction = "y") +
            
            # Corrected
            geom_line(data = subset(DataCor, !is.na(get(CorCol))),
                      aes(y = get(CorCol),
                          color = ifelse(Diameter != get(CorCol), 'Corrected', 'Conserved'))) +
            geom_point(data = subset(DataCor, !is.na(get(CorCol))),
                       aes(y = get(CorCol),
                           color = ifelse(Diameter != get(CorCol) | is.na(Diameter), 'Corrected', 'Conserved')),
                       shape = "circle", size = 3.9) +
            
            ggrepel::geom_text_repel(data = subset(DataCor,
                                                   (!is.na(get(CorCol)) & !is.na(get(POMcorv)) & (Diameter != get(CorCol)) | is.na(Diameter))),
                                     aes(y = get(CorCol), label = get(POMcorv), colour = "HOM"),
                                     point.size = 3.9, size = 3, direction = "y") +
            ggrepel::geom_text_repel(data = subset(DataCor, (!is.na(get(CorCol)) & DiameterCorrectionMeth != "")),
                                     aes(y = get(CorCol), label = DiameterCorrectionMeth, colour = "Methode"),
                                     point.size = 10, size = 3) +
            
            # Colours
            scale_colour_manual(name = "Status", values = c("Initial" = "red",
                                                            "Corrected" = "forestgreen",
                                                            "Conserved" = "black",
                                                            {if(nrow(subset(DataCor, !is.na(Diameter) & is.na(get(CorCol)))) > 0)
                                                              "Duplicated measurement" = "grey" },
                                                            "Methode" = "purple",
                                                            "HOM" = "blue")) +
            theme_minimal() +
            
            # Titles
            labs(
              # title =  paste("ID: ",unique(DataCor[, get(ID)]),""),
              x = "Year", y = "Diameter (cm)") +
            
            
            ggforce::facet_wrap_paginate(vars(get(ID), ScientificName),
                                         scales = "free",
                                         ncol = min(n,3), nrow = i, page = p)
    )
    
    if(SeveralWindows == TRUE & p < ceiling(length(unique(IDCor))/9))
      dev.new()
  }
  
  # return(Pl)
  
}