#=============== Calculate Carbon using Biomass ==========================
# ==== Inputs:
# df: data frame containing columns; Genus, Species, DBH (in cm), Height (even if NA for all entries)
# coords: either the coordinates of the site, a vector of longitude and latitude
#         or a matrix of coordinates for each tree
# ==== Optional Inputs:
# region: see ?getWoodDensity for the list of regions
# output.all: if TRUE outputs all data from processing

biomass <- function(df, coords, region = "World", output.all = TRUE){
  
  # Call BIOMASS package
  if(nchar(system.file(package='BIOMASS')) == 0 ){ 
    install.packages("BIOMASS",dependencies = TRUE)}
  library("BIOMASS")

  # Correct name
  correct <- correctTaxo(genus = df$Genus, species = df$Species)
  df$Genus.corrected   <- correct$genusCorrected
  df$Species.corrected <- correct$speciesCorrected
  df$Modified          <- correct$nameModified
  
  # Get Wood Density  # data("wdData") # Global wd database
  wd <- getWoodDensity(df$Genus.corrected, df$Species.corrected, region=region) 
  df$Wood.Density    <- wd$meanWD
  df$Wood.Density.sd <- wd$sdWD
  df$Family          <- wd$family
  
  # Get height estimates using Chave E
  h <- retrieveH(D = as.numeric(df$DBH), coord = coords)
  df$Height.est <- h[["H"]]
  df$RSE        <- h$RSE
  
  # Combine height data with height estimates
  df$Height.1 <- NA
  h.data <- df[!is.na(df$Height),]
  h.est  <- df[is.na(df$Height),]
  h.data$Height.1  <- h.data$Height
  h.est$Height.1   <- h.est$Height.est

  df <- rbind(h.data, h.est)
  
  # Calculate Carbon using Biomass package.
  df$AGB.Biomass.kg <- computeAGB(D=as.numeric(df$DBH),
                                  WD=as.numeric(df$Wood.Density),
                                  H=df$Height.1)*1000
  
  if(output.all == FALSE){
    df <- df[c('Genus.corrected','Species.corrected','Family','Latitude',
               'Longitude','DBH','AGB.Biomass.kg')]}
  
  
  return(df)
}

#=============== Calculate Carbon using allodb ==========================
# ==== Inputs:
# df: data frame containing columns; Genus.corrected, Species.corrected, DBH (in cm)
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
  df$AGB.allodb <- get_biomass(dbh = as.numeric(df$DBH), 
                               genus = df$Genus.corrected,
                               species = df$Species.corrected,
                               coords = coords, 
                               new_eqtable = new.eqtable)
  
  if(output.all == TRUE){
    df$allodb_a <- df$allodb_b <- df$allodb_sigma <- NA
    
    # Get parameters and sigma: AGB= a*DBH^b+e {e ~ N(0,sigma^2}
    params <- est_params(genus = df$Genus.corrected,
                         species = df$Species.corrected,
                         coords = coords)
    # Create a name column in params and df
    params$name <- paste(params$genus, params$species)
    df$Name <- paste(df$Genus.corrected, df$Species.corrected)
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



