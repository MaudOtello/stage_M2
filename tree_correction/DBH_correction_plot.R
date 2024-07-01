##Title : DBH à corriger
## Autor : Maud OTELLO 
## Date : 11/06/2024

# Objectif du script faire apparaître graphiquement les corrections apportées.

# ---- Packages utilisées ----
library(tidyverse) # Pour le langage et les fonction associées
library(sf)        # Pour manipuler les shp

# ---- (A sauter) Préparation des données ----
# # C'est pour simplement que tu vois les différentes manips et filtres réalisées sur le jeu de données
# # Importation des données
# Guyafor_Juveniles <- read.csv("2024ParacouJuveniles(v2).csv", sep=";", comment.char="#")
# 
# # Mise en forme du df
# Juveniles <- Guyafor_Juveniles %>%
#   filter(Protocole == "Regeneration")%>%
#   mutate(scientificName = paste(Genus, Species, sep = " "),
#          DBH = Circ/pi)
# 
# # Application des filtres
# Juveniles <- Juveniles %>%
#   filter(!is.na(Xutm) & 
#            Plot %in% c(1, 6, 11) & 
#            CodeAlive == 1) %>%
#   st_as_sf(coords = c("Xutm", "Yutm"), 
#            crs = st_crs(2972))
# 
# #importation des fichiers shapefiles
# Placettes <- st_read('RegenerationPlots.shp')
# # Sélection des centres 'réels' 
# PlacettesP1_6_11 <- Placettes %>% 
#   filter(Plot %in% c(1,6,11) &
#            CoordType == 'RealCoord')  
# 
# # Réajustement du crs 
# st_crs(PlacettesP1_6_11) <- st_crs(2972) 
# 
# # buffer 5m autour du centre des placettes 
# PlacettesP1_6_11_buffer5m <- st_buffer(PlacettesP1_6_11,5)  
# 
# # Fusion des placettes aux Juveniles
# Juveniles_placette <- st_intersection(PlacettesP1_6_11_buffer5m, Juveniles)
# 
# # Pas besoin de la spatialisation anymore + Guyafor pas 2002 + pas 2003, P1 & 6
# Juveniles_placette <- Juveniles_placette %>%
#   filter(CensusYear != 2002)%>%
#   filter(CensusYear != 2003)%>%
#   st_drop_geometry() # Geometry inutile à présent, retour à un tibble

# ----- Correction des DBH ----- 
# Importation fichier
Juveniles_placette <- read.csv("Juveniles_placette.csv", sep=";", comment.char="#")

# Sélection des espèces d'intérêts
esp_6 <- c("Dicorynia guianensis", "Pradosia cochlearia", "Iryanthera sagotiana", "Qualea rosea", "Vouacapoua americana", "Jacaranda copaia")

Juveniles_placette <- Juveniles_placette %>%
  filter(scientificName %in% esp_6) 

# Préparation des données pour correction
Juveniles_acorr <- Juveniles_placette %>%
  select(Plot, idTree, idRegePlot,scientificName, DBH,CensusYear,POM,CodeAlive)

colnames(Juveniles_acorr) <- c("Plot","IdTree","idRegePlot", "ScientificName_TreeDataCor", "Diameter", "Year", "POM","LifeStatus")

Juveniles_acorr <- Juveniles_acorr %>%
  mutate(
    IdTree = as.character(IdTree),
    ScientificName_TreeDataCor = as.character(ScientificName_TreeDataCor),
    Diameter = as.numeric(Diameter),
    Year = as.factor(Year),
    POM = as.numeric(POM)
  )
# application fonction  
Juveniles_corr <- DiameterCorrection(Juveniles_acorr, NegativeGrowthThreshold = 0,
                                     WhatToCorrect = c("punctual", "shift"),
                                     CorrectionType = c("individual"), 
                                     DBHCorForDeadTrees = F)

# Cette fonction arrondie au 1 chiffre après la virgule, voir comment changer ça

# ----- Plot de la correction des diamètres ----
# Obligation d'utiliser les POM les mettre à 130 directement
