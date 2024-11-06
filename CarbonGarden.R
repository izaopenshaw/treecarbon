#================ Carbon Garden================
#================ Packages ================
library(stringr)
library(ggplot2)

# For updates, remove and reinstall
remove.packages("WoodlandCarbonCode")
devtools::install_github("izaopenshaw/WoodlandCarbonCode", force = TRUE) # force = TRUE to install updates
library(WoodlandCarbonCode)

#================ Import ================
setwd("C:/Users/ipr10kg/OneDrive - The Royal Botanic Gardens, Kew/Documents/Methodology")
df <- read.csv("OG_Data/CarbonGarden.csv")
df$DBH <- as.numeric(str_trim(df$DBH..cm.., side = c("both", "left", "right")))
df$Height <- as.numeric(str_trim(df$Height..m.., side = c("both", "left", "right")))
df$Itrees_kg <- as.numeric(str_trim(df$Itrees.carbon..kg.., side = c("both", "left", "right")))

#================ Calculate Carbon ================
coords <- c(-0.088837,51.071610)
df <- biomass(df, coords, region="Europe", output.all = TRUE)
df <- allodb(df, coords, output.all = TRUE)

df$CVF[df$Classification=="broadleaf"] <- 0.488
df$CVF[df$Classification=="conifer"] <- 0.508 # Error = 0.6
# Thomas, Sean C., and Adam R. Martin. "Carbon content of tree tissues: a synthesis." Forests 3.2 (2012): 332-352. https://www.mdpi.com/1999-4907/3/2/332.

df$AGC_Biomass_kg <- df$AGB_Biomass_kg*df$CVF
df$AGC_allodb_kg  <- df$AGB_allodb_kg*df$CVF
df$allodb_sd <- df$allodb_sigma*df$CVF

sp_codes <- lookspcode(as.character(df$Name), df$Classification, returnv = 'single')
carbon <- fc_agc_error(sp_codes$spcode, df$DBH, df$Height)

df$AGC_WCC <- carbon$AGC*1000
df$WCC_sigma <- carbon$sig_AGC*1000

# Replace negative values in df$AGC_WCC with NA
df$AGC_WCC[df$AGC_WCC < 0] <- NA

#================ Plot ================
# Plot all
carbon_data <- data.frame(
  Carbon_Type = c("Itrees", "Biomass", "allodb", "WCC"),
  Carbon_Value = c(sum(df$Itrees_kg , na.rm = TRUE),
                   sum(df$AGC_Biomass_kg, na.rm = TRUE),
                   sum(df$AGC_allodb_kg, na.rm = TRUE),
                   sum(df$AGC_WCC, na.rm = TRUE)),
  Sigma = c(mean(df$allodb_sigma, na.rm = TRUE),
            mean(df$allodb_sigma, na.rm = TRUE),
            mean(df$allodb_sigma, na.rm = TRUE),
            mean(sqrt(df$WCC_sigma), na.rm = TRUE))
)

ggplot(carbon_data, aes(x = Carbon_Type, y = Carbon_Value, fill = Carbon_Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
#  geom_errorbar(aes(ymin = Carbon_Value - Sigma, ymax = Carbon_Value + Sigma), width = 0.2) +
  labs(title = "Carbon Garden Total Carbon",
       x = "Carbon Type",
       y = "Total Carbon (kg)") +
  theme_minimal()

# Plot one tree
i <- 1
tree_data <- df[i,]

carbon_data <- data.frame(
  Carbon_Type = c("Itrees_kg ", "AGC_Biomass", "AGC_allodb", "AGC_WCC"),
  Carbon_Value = c(sum(df$Itrees_kg , na.rm = TRUE),
                   sum(df$AGC_Biomass_kg, na.rm = TRUE),
                   sum(df$AGC_allodb_kg, na.rm = TRUE),
                   sum(df$AGC_WCC, na.rm = TRUE)),
  Sigma = c(mean(df$allodb_sigma, na.rm = TRUE),
            mean(df$allodb_sigma, na.rm = TRUE),
            mean(df$allodb_sigma, na.rm = TRUE),
            mean(df$allodb_sigma, na.rm = TRUE))
)

# Plot the data for the selected tree
ggplot(carbon_data, aes(x = Carbon_Type, y = Carbon_Value, fill = Carbon_Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = Carbon_Value - Sigma, ymax = Carbon_Value + Sigma), width = 0.2) +
  labs(
    title = paste("Total Carbon for Tree ID:", i),
    x = "Carbon Type",
    y = "Carbon (kg)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

