# SuessR Package


# Function to calculate the Laws correction. Used within the SuessR() and SuessR.custom() functions.

laws.fun <- function(e1 = 1, e2 = 26.5,e.1 = -1, laws.CO2, P, sst, C, b) {
  e2 + e1 - e.1 - (1/(1+((laws.CO2*P)/ (0.5*C*(1+b))))) * ((e2 - e.1)/ (b+1))
}


# Function to calculate and apply Suess and Laws corrections.

#' @title Calculate Suess and Laws corrections for data from a built-in region
#' @description The SuessR() function generates region-specific Suess, Laws, and net (Suess + Laws) corrections for
#' d13C data input by the user. The net correction is then used to calculate the corrected d13C data, which are supplied
#' in the output. This function is specifically for data from regions currently built into SuessR
#' ("Bering", "Aleutians", and "Gulf of Alaska" as of February 2020).
#' @param data A dataset including sample ID ('id'), year of sample collection ('year'), uncorrected d13C data ('d13C'),
#' and region ('region').
#' @param correct.to The year to which the d13C data will be corrected. Defaults to 1850, pre-Suess effect.
#' @details Figure out what should go in the 'DETAILS' section for these functions.
#' @return The output of this function is a data frame that includes the sample ID ('id'), year ('year'),
#' uncorrected d13C ('d13c.uncor'), Laws correction ('Laws.cor'), Suess Correction ('Suess.cor'), net correction
#' ('net.cor'), and corrected d13C ('d13c.cor') for each sample.
#' @examples
#' example.data <- data.frame(id = c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6"),
#'                            year = c(2017, 2017, 2017, 1977, 1977, 1977),
#'                            d13c = c(-12, -12, -12, -12, -12, -12),
#'                            region = c("Bering", "Aleutians", "Gulf of Alaska", "Bering", "Aleutians", "Gulf of Alaska"))
#'
#' SuessR(data = example.data)



SuessR <- function(data, correct.to = 1850) {

  # Make sure the data is a dataframe
  if(is.matrix(data)) data <- as.data.frame(data)

  # Screen input data to check for (and remove) missing values, with warning
  if(anyNA(data)) warning("Rows with missing values were deleted")
  data <- na.omit(data)

  # Screen input data to make sure all years provided are integers between 1850 and the current year -1
  current.year <- as.numeric(substr(date(), 21,24))-1
  if(any(!(data$year %in% 1850:current.year)))  stop(paste("Year must be an integer between 1850 and ", current.year))

  # Screen input data to make sure the 'correct.to' year is an integer between 1850 and the current year -1
  current.year <- as.numeric(substr(date(), 21,24))-1
  if(any(!(correct.to %in% 1850:current.year)))  stop(paste("correct.to must be an integer between 1850 and ", current.year))

  # Screen input data to check for non-negative d13c values
  if(any(data$d13c >= 0)) warning(paste("Some d13c value(s) are non-negative"))

  # Screen input data to check for (and remove) missing values, with warning
  if(any(!(data$region %in% SuessR.reference.data$region))) stop(paste("Unrecognized region:", unique(data$region[!(data$region %in% SuessR.reference.data$region)]), "  "))



  # Create output data frame
  SuessR.out <- data.frame(id = data$id, year = data$year)

  ref <- SuessR.reference.data

  #ln(K0)
  ref$lnK0 <- with(ref,
                   -58.0931 + 90.5069 * (100/(sst+273.15))
                   + (22.294 * log((sst+273.15)/100))
                   + (sss * (0.027766 + -0.025888 * ((sst+273.15)/100)
                           + 0.005058 * ((sst+273.15)/100)^2)))

  #Ocean Increase
  n <- nrow(ref)
  ref$oi <- c(NA, with(ref, (CO2atm[-1] - CO2atm[-n])*0.4))

  #~f(CO2)ocean
  ref$fCO2 <- c(284.25, cumsum(ref$oi[-1]) + 284.25)

  #CO2aq
  ref$CO2aq <- exp(ref$lnK0)*ref$fCO2

  # Laws expression for a given year
  ref$laws.current <- with(ref, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))

  data <- merge(data, ref, c("region", "year"), sort = F)

  dat1850 <- ref[ref$year==1850,]
  dat1850$laws1850 <- with(dat1850, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat1850[,c("region", "laws1850")], "region", sort = F)
  dat.correct.to <- ref[ref$year==correct.to,]
  dat.correct.to$laws.correct.to <- with(dat.correct.to, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat.correct.to[,c("region", "laws.correct.to")], "region", sort = F)
  data$Laws.cor <- round(with(data, (laws.current - laws1850) - (laws.correct.to - laws1850)),2)

  SuessR.out <- data[,c("id","year", "d13c", "Laws.cor")]
  names(SuessR.out)[3] <- "d13c.uncor"
  SuessR.out$Suess.cor  <- round(with(data, 0.014 * exp((year-1850)*0.027)
                                - 0.014*exp((correct.to-1850)*0.027)),2)
  SuessR.out$net.cor    <- SuessR.out$Suess.cor + SuessR.out$Laws.cor
  SuessR.out$d13c.cor   <- data$d13c + SuessR.out$net.cor
  SuessR.out <- SuessR.out[order(rownames(data)),]
  print(SuessR.out)
}



#' @title Calculate Suess and Laws corrections for data from a custom region
#' @description The SuessR.custom() function generates region-specific Suess, Laws, and net (Suess + Laws) corrections for
#' d13C data input by the user. The net correction is then used to calculate the corrected d13C data, which are supplied
#' in the output. This function is specifically for data from regions not currently built into SuessR
#' (i.e., different from "Bering", "Aleutians", and "Gulf of Alaska" as of February 2020).
#'
#' @param data A dataframe including sample ID, year of sample collection, uncorrected d13C
#' data, and region. This function is specifically for data from regions not
#' currently built into SuessR. Columns must be named 'id', 'year', 'd13c',
#' and 'region'.
#' @param custom.region.data A dataframe containing environmental data for the custom region from which the
#' samples originated. Must contain columns titled 'year', 'region', 'r', 'C', P', "mu', 'beta', 'sst', 'S',
#' and 'CO2atm'. See details for information on how to supply these parameters appropriately.
#' @param correct.to The year to which the d13C data will be corrected. Defaults to 1850, pre-Suess effect.
#' @details Figure out what should go in the 'DETAILS' section for these functions.
#' @return Figure out what should go in the 'VALUES' section for these functions.
#' @examples
#' example.region.data <- data.frame(year = seq(from = 1850, to = 2019, by = 1),
#'                                                          region = rep("Example Region", 170),
#'                                                          esub1 = rep(1, 170),
#'                                                          esubneg1 = rep(1, 170),
#'                                                          esub2 = rep(26.5, 170),
#'                                                          r = rep(5, 170),
#'                                                          C = rep(1.09e-09, 170),
#'                                                          P = rep(4.57e-10, 170),
#'                                                          mu = rep(0.5, 170),
#'                                                          beta = rep(0.2, 170),
#'                                                          sst = seq(5.9, 6.6, 170),
#'                                                          S = seq(32.3, 32.7, 170),
#'                                                          CO2atm = SuessR.reference.data$CO2atm[1:170])
#'
#' example.custom.data <- data.frame(id = c("Sample 1", "Sample 2", "Sample 3", "Sample 4"),
#'                                                           year = c(1850, 1900, 1950, 2000),
#'                                                           d13c = c(-12, -12, -12, -12),
#'                                                           region = rep("Example Region", 4)
#'
#' SuessR.custom(data = example.custom.data, custom.region.data = example.region.data, correct.to = 1850)


# Second function that allows users to add a custom region

SuessR.custom <- function(data, custom.region.data, correct.to = 1850) {


  # Make sure the data is a dataframe
  if(is.matrix(data)) data <- as.data.frame(data)

  # Screen input data to check for (and remove) missing values, with warning
  if(anyNA(data)) warning("Rows with missing values were deleted")
  data <- na.omit(data)

  # Screen input data to make sure all years provided are integers between 1850 and the current year -1
  current.year <- as.numeric(substr(date(), 21,24))-1
  if(any(!(data$year %in% 1850:current.year)))  stop(paste("Year must be an integer between 1850 and ", current.year))

  # Screen input data to make sure the 'correct.to' year is an integer between 1850 and the current year -1
  current.year <- as.numeric(substr(date(), 21,24))-1
  if(any(!(correct.to %in% 1850:current.year)))  stop(paste("correct.to must be an integer between 1850 and ", current.year))

  # Screen input data to check for non-negative d13c values
  if(any(data$d13c >= 0)) warning(paste("Some d13c value(s) are non-negative"))




  # Create output data frame
  SuessR.out <- data.frame(id = data$id, year = data$year)

  ref <- rbind(SuessR.reference.data, custom.region.data)


  # Screen input data to check for (and remove) missing values, with warning
  if(any(!(data$region %in% ref$region))) stop(paste("Unrecognized region: "), data$region[!(data$region %in% ref$region)])



  #ln(K0)
  ref$lnK0 <- with(ref,
                   -58.0931 + 90.5069 * (100/(sst + 273.15))
                   + (22.294 * log((sst + 273.15)/100))
                   + (sss * (0.027766 + -0.025888 * ((sst+273.15)/100)
                           + 0.005058 * ((sst + 273.15)/100)^2)))

  #Ocean Increase
  n <- nrow(ref)
  ref$oi <- c(NA, with(ref, (CO2atm[-1] - CO2atm[-n])*0.4))

  #~f(CO2)ocean
  ref$fCO2 <- c(284.25, cumsum(ref$oi[-1]) + 284.25)

  #CO2aq
  ref$CO2aq <- exp(ref$lnK0)*ref$fCO2

  # Laws expression for a given year
  ref$laws.current <- with(ref, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))

  data <- merge(data, ref, c("region", "year"), sort = F)

  dat1850 <- ref[ref$year==1850,]
  dat1850$laws1850 <- with(dat1850, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat1850[,c("region", "laws1850")], "region", sort = F)
  dat.correct.to <- ref[ref$year==correct.to,]
  dat.correct.to$laws.correct.to <- with(dat.correct.to, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat.correct.to[,c("region", "laws.correct.to")], "region", sort = F)
  data$Laws.cor <- round(with(data, (laws.current - laws1850) - (laws.correct.to - laws1850)), 2)

  SuessR.out <- data[,c("id","year", "d13c", "Laws.cor")]
  names(SuessR.out)[3] <- "d13c.uncor"
  SuessR.out$Suess.cor  <- round(with(data, 0.014 * exp((year-1850)*0.027)
                                - 0.014*exp((correct.to-1850)*0.027)),2)
  SuessR.out$net.cor    <- SuessR.out$Suess.cor + SuessR.out$Laws.cor
  SuessR.out$d13c.cor   <- data$d13c + SuessR.out$net.cor
  print(SuessR.out)
}




# Function for calculating regional uptake constants

#' Calculate a regional CO2 uptake constant from empirical DIC data
#'
#' @param year1 The year in which samples were collected for the first set of DIC d13C observations.
#' @param year2 The year in which samples were collected for the second set of DIC d13C observations.
#' @param d13c.change The observed change in DIC d13C between the two years.
#' @details Figure out what should go in the 'DETAILS' section for these functions.
#' @return Figure out what should go in the 'VALUES' section for these functions.
#' @examples
#' year1 <- 1970
#' year2 <- 1980
#' d13c.change <- 0.3
#'
#' reg.uptake(year1 = year1, year2 = year2, d13c.change = d13c.change)

reg.uptake <- function(year1, year2, d13c.change) {

  reg.up.const <- round(d13c.change/(exp((year2-1850)*0.027) - exp((year1 - 1850)*0.027)),2)

  print(reg.up.const)

}





