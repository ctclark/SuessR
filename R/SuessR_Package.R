# SuessR Package

# Function to calculate the Laws correction. Used within the SuessR() and SuessR.custom() functions.

laws.fun <- function(e1, e2, e.1, laws.CO2, sst, r, b, p) {

  e2 + e1 - e.1 - (1/(1+laws.CO2*((((p)*(86400)*4*pi*(r*10^-6)^2))/((0.89*1.47^((sst + 273.15 - 303.15)/10))*((0.216*((4/3)*pi*r^3)^0.939)*10^-12)*(1+b)))))*((e2 - e.1)/(b+1))

}


# Function to calculate and apply Suess and Laws corrections.

#' @encoding UTF-8
#' @title Calculate Suess and Laws corrections for stable carbon isotope data from a built-in region
#' @description Generates region-specific Suess, Laws, and net (Suess + Laws) corrections for stable carbon isotope data
#'   \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data input by the user. The net correction is then used to calculate the corrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data, which are supplied
#'   in the output. This function is specifically for data from regions currently built into SuessR
#'   ("Bering", "Aleutians","Gulf of Alaska", and "Subpolar North Atlantic", as of Version 1.0.0).
#'   \if{html}{\figure{Built-in-regions.png}{options: width=600 alt="Built-in regions"}}
#'   \if{latex}{\figure{Built-in-regions.png}{options: width=6in}}
#' @param data A matrix or data frame including columns containing sample ID ('id'), year of sample collection ('year'),
#'   region of sample collection ('region'), and uncorrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C
#'   data ('d13c').
#' @param correct.to The year to which the \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data will be corrected. Defaults to 1850 (Industrial Revolution).
#' @details The SuessR() function uses the year and region of sample collection to calculate Suess, Laws, and net (Suess + Laws)
#'   corrections for \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data from marine organisms. The Suess Correction represents the change in \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C values of dissolved
#'   inorganic carbon (DIC) in the surface ocean, and is calculated using an exponential function, calibrated to the global decline
#'   DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C values. Surface waters in different regions exhibit varying rates of CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} uptake from the atmosphere as a result of
#'   differing water mass properties and residence time at the surface. Thus, a region-specific modifier has been
#'   calculated from empirical observations of changes in DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C values through time and is applied to each region. The Laws
#'   correction accounts for changes in stable carbon isotope fractionation during CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} uptake and photosynthesis by phytoplankton,
#'   and is impacted by changes in aqueous CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} concentrations, temperature, and salinity, as well as community growth rates,
#'   average cell diameter, average organic carbon content of phytoplankton cells, permeability of the phytoplankton plasmalemma to CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}},
#'   and the ratio of net diffusional loss of CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} to carbon fixation. Historic observations and reconstructions of atmospheric CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} concentrations,
#'   sea surface temperature, and sea surface salinity are used to make these calculations. See references for more details.
#'
#' This function allows users to correct their data to a specific year, using the 'correct.to' argument. This argument defaults
#'   to AD1850, which represents onset of the Industrial Revolution and the initiation of the large-scale combustion of fossil
#'   fuels that created the Suess Effect. Correcting data to 1850 allows them to be compared to data from any samples collected during
#'   or before 1850 (including archaeological samples), as well as to any other samples Suess corrected to the year 1850. Users
#'   examining contemporary datasets may wish to correct their data to another year. For example, a user examining a dataset
#'   of \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data spanning the years 1970-2010 might choose to correct all their data back to the year 1970, to correct the
#'   older samples forward to 2010, or to correct all the samples to 1990, the middle of the time series. In any of these scenarios,
#'   corrected data would be comparable to one another. When choosing a value for 'correct.to' consideration should be given to
#'   the comparability of the user's data to existing data sets, clarity of presentation of results (i.e., clear statements of the
#'   year to which the data were corrected), and reproducibility of results (i.e., presentation of both uncorrected and corrected
#'   data so future users can repeat analyses or correct the raw data to another year).
#' @return The output of this function is a data frame that includes the sample ID ('id'), year ('year'),
#'   uncorrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C ('d13c.uncor'), Laws correction ('Laws.cor'), Suess Correction ('Suess.cor'), net correction
#'   ('net.cor'), and corrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C ('d13c.cor') for each sample. The corrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data is equal to the uncorrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data
#'   plus the net correction. The units for all values are the standard 'per mil' used for \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data.
#' @examples
#' example.data <- data.frame(id = c("Sample 1", "Sample 2", "Sample 3", "Sample 4",
#'                                   "Sample 5", "Sample 6", "Sample 7", "Sample 8"),
#'                            year = c(2017, 2017, 2017, 2017, 1977, 1977, 1977, 1977),
#'                            d13c = c(-12, -12, -12, -12, -12, -12, -12, -12),
#'                            region = c("Bering Sea", "Aleutian Islands",
#'                                       "Gulf of Alaska", "Subpolar North Atlantic",
#'                                       "Bering Sea", "Aleutian Islands",
#'                                       "Gulf of Alaska", "Subpolar North Atlantic"))
#'
#' SuessR(data = example.data)
#' @references Clark, C.T., M.R. Cape, M.D. Shapley, F.J. Mueter, B.P. Finney, and N. Misarti. (In Prep) SuessR: Regional Suess
#'   and Laws corrections for \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C from marine organisms.
#'
#' Clark, C.T., L. Horstmann, A. de Vernal, A.M. Jensen, and N. Misarti. (2019) Pacific walrus diet across 4000 years of
#'   changing sea ice conditions. \emph{Quaternary Research}, 1-17.
#'
#' Misarti, N., B. Finney, H. Maschner, and M.J. Wooller. (2009) Changes in northeast Pacific marine ecosystems over the last 4500
#'   years: evidence from stable isotope analysis of bone collagen from archaeological middens. \emph{The Holocene}, 19:8. 1139-1151.
#'
#' @importFrom stats na.omit
#' @export

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
  for(i in unique(ref$region)) {
    n <- nrow(ref[ref$region == i,])
    ref$oi <- c(NA, with(ref[ref$region == i,], (CO2atm[-1] - CO2atm[-n])*ref[ref$region == i,]$Cp[1]))
  }

  #~f(CO2)ocean
  for(i in unique(ref$region)) {
    ref$fCO2 <- c(285.78, cumsum(ref[ref$region == i,]$oi[-1]) + 285.78)
  }

  #CO2aq
  ref$CO2aq <- exp(ref$lnK0)*ref$fCO2

  # Laws expression for a given year
  ref$laws.current <- with(ref, laws.fun(e1=1, e2=26.5, e.1=1, p=1.5e-5, laws.CO2=CO2aq, sst=sst, r=r, b=0.2))

  data$order <- seq(1,length(data$id),1)
  data <- merge(data, ref, c("region", "year"), sort = F)

  dat1850 <- ref[ref$year==1850,]
  dat1850$laws1850 <- with(dat1850, laws.fun(e1=1, e2=26.5, e.1=1, p=1.5e-5, laws.CO2=CO2aq, sst=sst, r=r, b=0.2))

  data <- merge(data, dat1850[,c("region", "laws1850")], "region", sort = F)
  dat.correct.to <- ref[ref$year==correct.to,]
  dat.correct.to$laws.correct.to <- with(dat.correct.to, laws.fun(e1=1, e2=26.5, e.1=1, p=1.5e-5, laws.CO2=CO2aq, sst=sst, r=r, b=0.2))
  data <- merge(data, dat.correct.to[,c("region", "laws.correct.to")], "region", sort = F)
  data$Laws.cor <- round(with(data, (laws.current - laws1850) - (laws.correct.to - laws1850)),3)

  SuessR.out <- data[,c("id","year", "d13c", "Laws.cor")]
  names(SuessR.out)[3] <- "d13c.uncor"
  SuessR.out$Suess.cor  <- round(with(data, up.con*exp((year-1850)*0.027)
                                - up.con*exp((correct.to-1850)*0.027)),3)
  SuessR.out$net.cor    <- SuessR.out$Suess.cor + SuessR.out$Laws.cor
  SuessR.out$d13c.cor   <- data$d13c + SuessR.out$net.cor
  SuessR.out <- SuessR.out[order(data$order),]
  print(SuessR.out)
}


# A function that allows users to add custom regions for the SuessR calculations

#' @encoding UTF-8
#' @title Calculate Suess and Laws corrections for stable carbon isotope data from a custom region
#' @description The SuessR.custom() function generates region-specific Suess, Laws, and net (Suess + Laws) corrections for
#'   \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data input by the user. The net correction is then used to calculate the corrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data,
#'   which are supplied in the output. This function is specifically for data from regions not currently built into SuessR
#'   (i.e., different from "Bering", "Aleutians", "Gulf of Alaska", and "Subpolar North Atlantic", as of Version 1.0.0). Or
#'   for scenarios in which users want to change the default parameters used by SuessR for the built-in regions (e.g., averge phytoplankton cell radius).
#' @param data A dataframe including sample ID, year of sample collection, uncorrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C
#'   data, and region. This function is specifically for data from regions not
#'   currently built into SuessR. Columns must be named 'id', 'year', 'd13c', and 'region'.
#' @param custom.region.data A data frame containing environmental data for the custom region from which the
#'   samples originated. Must contain columns titled 'year', 'region', 'r', 'sst', 'sss', 'CO2atm', 'up.con', and 'Cp'.
#'   See details for information on how to supply these parameters appropriately.
#' @param correct.to The year to which the \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data will be corrected. Defaults to 1850, pre-Suess effect.
#' @details The SuessR.custom() allows users to calculate and apply Suess and Laws corrections to \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data from marine
#'   organisms collected in a region not currently built into the SuessR package. In the initial release (Version 0.1.0), the
#'   built-in regions are the Bering Sea ('Bering'), the Aleutian archipelago ('Aleutians'), and the Gulf of Alaska ('Gulf of
#'   Alaska'). Because the Suess and Laws corrections require region-specific environmental data from 1850-present, users must
#'   supply these data using the 'custom.region.data' argument. Once these data have been supplied, this function calculates
#'   the corrections exactly like the SuessR() function. See the built-in 'SuessR.reference.data' object for an example
#'   template for the data, as well as the references for detailed information on how to compile and supply the appropriate
#'   data to the SuessR.custom() function. After compiling these data, please consider sharing them with the package authors
#'   (ctclark 'at' alaska.edu) to be included as built-in regions in future version of this package.
#'
#' As with the SuessR() function, SuessR.custom() allows users to correct their data to a specific year, using the 'correct.to'
#'   argument. This argument defaults to AD1850, which represents onset of the Industrial Revolution and the initiation of the
#'   large-scale combusion of fossil fuels that created the Suess Effect. Correcting data to 1850 allows them to be compared to
#'   data from any samples collected during or before 1850 (including archaeological samples), as well as to any other samples
#'   Suess corrected to the year 1850. Users examining contemporary datasets may wish to correct their data to another time year.
#'   For example, a user examining a dataset of \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data spanning the years 1970-2010 might choose to correct all their data back
#'   to the year 1970, to correct the older samples forward to 2010, or to correct all the samples to 1990, the middle of the time
#'   series. In any of these scenarios, corrected data would be comparable to one another. When choosing a value for 'correct.to'
#'   consideration should be given to the comparability of the user's data to existing data sets, clarity of presentation of
#'   results (i.e., clear statements of the year to which the data were corrected), and reproducibility of results (i.e.,
#'   presentation of both uncorrected and corrected data so future users can repeat analyses or correct the raw data to another year).
#' @return The output of this function is a data frame that includes the sample ID ('id'), year ('year'),
#'   uncorrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C ('d13c.uncor'), Laws correction ('Laws.cor'), Suess Correction ('Suess.cor'), net correction
#'   ('net.cor'), and corrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C ('d13c.cor') for each sample. The corrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data is equal to the uncorrected \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data
#'   plus the net correction. The units for all values are the standard 'per mil' used for \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C data.
#' @examples
#' example.region.data <- data.frame(year = seq(from = 1850, to = 2019, by = 1),
#'                                   region = rep("Example Region", 170),
#'                                   r = rep(5, 170),
#'                                   sst = seq(5.9, 6.6, 170),
#'                                   sss = seq(32.3, 32.7, 170),
#'                                   CO2atm = SuessR.reference.data$CO2atm,
#'                                   up.con = rep(0.014, 170),
#'                                   Cp = rep(0.41, 170))
#'
#' example.custom.data <- data.frame(id = c("Sample 1", "Sample 2",
#'                                          "Sample 3", "Sample 4"),
#'                                   year = c(1850, 1900, 1950, 2000),
#'                                   d13c = c(-12, -12, -12, -12),
#'                                   region = rep("Example Region", 4))
#'
#' SuessR.custom(data = example.custom.data,
#'               custom.region.data = example.region.data,
#'               correct.to = 1850)
#'
#' @importFrom stats na.omit
#' @export


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
  for(i in unique(ref$region)) {
  n <- nrow(ref[ref$region == i,])
  ref$oi <- c(0, with(ref[ref$region == i,], (CO2atm[-1] - CO2atm[-n])*ref[ref$region == i,]$Cp[1]))
  }


  #~f(CO2)ocean
  for(i in unique(ref$region)) {
    ref$fCO2 <- c(285.78, cumsum(ref[ref$region == i,]$oi[-1]) + 285.78)
  }


  #CO2aq
  ref$CO2aq <- exp(ref$lnK0)*ref$fCO2

  # Laws expression for a given year
  ref$laws.current <- with(ref, laws.fun(e1=1, e2=26.5, e.1=1, p=1.5e-5, laws.CO2=CO2aq, sst=sst, r=r, b=0.2))

  data$order <- seq(1,length(data$id),1)
  data <- merge(data, ref, c("region", "year"), sort = F)

  dat1850 <- ref[ref$year==1850,]
  dat1850$laws1850 <- with(dat1850, laws.fun(e1=1, e2=26.5, e.1=1, p=1.5e-5, laws.CO2=CO2aq, sst=sst, r=r, b=0.2))
  data <- merge(data, dat1850[,c("region", "laws1850")], "region", sort = F)
  dat.correct.to <- ref[ref$year==correct.to,]
  dat.correct.to$laws.correct.to <- with(dat.correct.to, laws.fun(e1=1, e2=26.5, e.1=1, p=1.5e-5, laws.CO2=CO2aq, sst=sst, r=r, b=0.2))
  data <- merge(data, dat.correct.to[,c("region", "laws.correct.to")], "region", sort = F)
  data$Laws.cor <- round(with(data, (laws.current - laws1850) - (laws.correct.to - laws1850)), 3)

  SuessR.out <- data[,c("id","year", "d13c", "Laws.cor")]
  names(SuessR.out)[3] <- "d13c.uncor"
  SuessR.out$Suess.cor  <- round(with(data, up.con*exp((year-1850)*0.027)
                                - up.con*exp((correct.to-1850)*0.027)),3)
  SuessR.out$net.cor    <- SuessR.out$Suess.cor + SuessR.out$Laws.cor
  SuessR.out$d13c.cor   <- data$d13c + SuessR.out$net.cor
  SuessR.out <- SuessR.out[order(data$order),]
  print(SuessR.out)
}




# Function for calculating regional uptake constants

#' @encoding UTF-8
#' @title Calculate a regional carbon dioxide uptake constant from empirical DIC data
#' @description The reg.uptake() function calculates the regional uptake constant required to modify the global Suess effect
#'   curve to be region-specific.
#' @param year1 The year in which samples were collected for the first set of DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C \cr observations.
#' @param year2 The year in which samples were collected for the second set of DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C \cr observations.
#' @param d13c.change The observed change in DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C between the two years.
#' @details This function calculates the regional uptake constant used to modify the global Suess effect curve to make it
#'   specific to a region. This method requires \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C values for DIC from the area of interest from two different time points,
#'   ideally separated by at least a decade. The function uses the magnitude (i.e., absolute value) of the observed change in
#'   DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C ('d13c.change'), the year of the first observation ('year1'), and the year of the second observation ('year2')
#'   to calculate the regional uptake constant. This value can then be supplied as part of the 'custom.region.data' argument
#'   (filling a column titled 'up.con') for the SuessR.custom() function.
#' @return Returns a numerical value representing the regional uptake constant. This value can be supplied as part of the
#'   'custom.region.data' (filling a column titled 'up.con') argument for the SuessR.custom() function.
#' @examples
#' year1 <- 1970
#' year2 <- 1980
#' d13c.change <- 0.3
#'
#' reg.uptake(year1 = year1, year2 = year2, d13c.change = d13c.change)
#'
#' @export

reg.uptake <- function(year1, year2, d13c.change) {

  reg.up.const <- round(d13c.change/(exp((year2-1850)*0.027) - exp((year1 - 1850)*0.027)),3)

  print(reg.up.const)

}





