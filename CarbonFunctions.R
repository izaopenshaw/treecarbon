#=============== Calculate Carbon using Biomass ==========================
#'
#' @title Estimate Tree Carbon using Biomass package functions
#' @description Using the Biomass package to calculate carbon
#' @param df dataframe containing columns; Genus, Species, DBH (in cm),
#' Height (in m). Height is optional, but df must include this column.
#' @param coords either a vector of coordinates of the site or a matrix of
#' coordinates for each tree of longitude and latitude
#' @param region of the World. See ?getWoodDensity for the list of regions
#' @param output.all if TRUE outputs all data from processing, else just outputs carbon figures
#' @returns  your dataframe back with added columns of carbon estimates. If
#' output.all = FALSE, then returns columns 'Genus_corrected','Species_corrected',
#' 'Family','Latitude','Longitude','DBH','AGB_Biomass_kg'. If output.all = TRUE
#' then additionally returns columns 'Wood_Density', 'Wood_Density_sd',
#' 'Height_est', 'RSE' (Residual Standard Error of the model), 'Height_1' (which
#' is inputed height filled in with Height estimate where missing)
#'
# WD (in g/cm3) to add optionally or lookup from global tables?

biomass <- function(df, coords, region = "World", output.all = TRUE){

  # Call BIOMASS package
  if(nchar(system.file(package='BIOMASS')) == 0 ){
    install.packages("BIOMASS",dependencies = TRUE)}
  library("BIOMASS")

  # Correct name
  correct <- correctTaxo(genus = df$Genus, species = df$Species)
  df$Genus_corrected   <- correct$genusCorrected
  df$Species_corrected <- correct$speciesCorrected
  df$Modified          <- correct$nameModified

  # Get Wood Density
  # data("wdData") # Global wd database
  wd <- getWoodDensity(df$Genus_corrected, df$Species_corrected, region=region)
  df$Wood_Density    <- wd$meanWD
  df$Wood_Density_sd <- wd$sdWD
  df$Family          <- wd$family

  # Get height estimates using Chave E
  h <- retrieveH(D = as.numeric(df$DBH), coord = coords)
  df$RSE        <- h$RSE
  df$Height_est <- h[["H"]]

  # Combine height data with height estimates
  df$Height_1 <- NA
  h.data <- df[!is.na(df$Height),]
  h.est  <- df[is.na(df$Height),]
  h.data$Height_1  <- h.data$Height
  h.est$Height_1   <- h.est$Height_est

  df <- rbind(h.data, h.est)

  # Calculate Carbon using Biomass package
  df$AGB_Biomass_kg <- computeAGB(D=as.numeric(df$DBH),
                                  WD=as.numeric(df$Wood_Density),
                                  H=df$Height_1)*1000

  if(output.all == FALSE){
    df <- df[c('Genus_corrected','Species_corrected','Family','Latitude',
               'Longitude','DBH','AGB_Biomass_kg')]}


  return(df)
}

#=============== Calculate Carbon using allodb ==========================
# ==== Inputs:
# df: data frame containing columns; DBH (in cm), Genus_corrected, Species_corrected (from output of biomass function)
# coords: either the coordinates of the site, a vector of longitude and latitude
#         or a matrix of coordinates for each tree
# ==== Optional Inputs:
# output.all: if TRUE outputs the coefficients of the model, a*DBH^b+e {e ~ N(0,sigma^2}
# new.eqtable: a subset or extension of the allometric equation table. Create with new_equations
allodb <- function(df, coords, output.all = TRUE, new.eqtable = NULL){

  # Call allodb package
  if(nchar(system.file(package='allodb')) == 0 ){
    remotes::install_github("ropensci/allodb")}
  library("allodb")

  # Biomass for all data. By default all equations will be used
  df$AGB_allodb_kg <- get_biomass(dbh = as.numeric(df$DBH),
                               genus = df$Genus_corrected,
                               species = df$Species_corrected,
                               coords = coords,
                               new_eqtable = new.eqtable)

  if(output.all == TRUE){
    df$allodb_a <- df$allodb_b <- df$allodb_sigma <- NA

    # Get parameters and sigma: AGB= a*DBH^b+e {e ~ N(0,sigma^2}
    params <- est_params(genus = df$Genus_corrected,
                         species = df$Species_corrected,
                         coords = coords)
    # Create a name column in params and df
    params$name <- paste(params$genus, params$species)
    df$Name <- paste(df$Genus_corrected, df$Species_corrected)
    df$Name <- as.factor(df$Name)
    # Create an empty df with same columns
    df.p <- df[-c(1:nrow(df)),]
    for(i in 1:length(levels(df$Name))){
      name <- levels(df$Name)[i]
      trees <- df[df$Name == name,]
      trees$allodb_a <- params$a[params$name == name]
      trees$allodb_b <- params$b[params$name == name]
      trees$allodb_sigma <-  params$sigma[params$name == name]
      df.p <- rbind(df.p, trees)
    }
    df <- df.p
  }

  return(df)
}



