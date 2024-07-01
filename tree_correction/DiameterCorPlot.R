library(tidyverse)
library(data.table)

Data <- setDT(Juveniles_esp_6_corr)
OnlyCorrected = T
CorCol = "Diameter_TreeDataCor"
SeveralWindows = F


# Order IDs and times in ascending order ----------------------------------------------------------------------------
ID <- "IdTree"
POMv <- "POM"
POMcorv <- "POM_TreeDataCor"
Data <- Data[order(get(ID), Year)]

if(OnlyCorrected == TRUE){
  # Only corrected stems ----------------------------------------------------------------------------------------------
  Data[,Diameter := round(Diameter, 1L)] # Digits
  Data[,Diameter_TreeDataCor := round(Diameter_TreeDataCor, 1L)] # Digits
  
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

pdf("DiameterCorrectionPlots_6sp.pdf", width = 25, height = 10)

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


dev.off()
