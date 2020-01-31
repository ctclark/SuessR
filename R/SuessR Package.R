# SuessR Package


# Function to calculate the Laws correction. Used within the SuessR() and SuessR.custom() functions.

laws.fun <- function(e1,e2,e.1,laws.CO2, P, sst, C, b) {
  e2 + e1 - e.1 - (1/(1+((laws.CO2*P)/ (0.5*C*(1+b))))) * ((e2 - e.1)/ (b+1))
}


# Function to calculate and apply Suess and Laws corrections.

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
                   + (S * (0.027766 + -0.025888 * ((sst+273.15)/100)
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

  data <- merge(data, ref, c("region", "year"))

  dat1850 <- ref[ref$year==1850,]
  dat1850$laws1850 <- with(dat1850, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat1850[,c("region", "laws1850")], "region")
  dat.correct.to <- ref[ref$year==correct.to,]
  dat.correct.to$laws.correct.to <- with(dat.correct.to, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat.correct.to[,c("region", "laws.correct.to")], "region")
  data$Laws.cor <- round(with(data, (laws.current - laws1850) - (laws.correct.to - laws1850)),2)

  SuessR.out <- data[,c("id","year", "d13c", "Laws.cor")]
  names(SuessR.out)[3] <- "d13c.uncor"
  SuessR.out$Suess.cor  <- round(with(data, 0.014 * exp((year-1850)*0.027)
                                - 0.014*exp((correct.to-1850)*0.027)),2)
  SuessR.out$net.cor    <- SuessR.out$Suess.cor + SuessR.out$Laws.cor
  SuessR.out$d13c.cor   <- data$d13c + SuessR.out$net.cor
  SuessR.out <- SuessR.out[order(SuessR.out$id),]  # Sort by id if desired
  print(SuessR.out)
}





# Second function that allows users to add a custom region

SuessR.custom <- function(data, custom.data, correct.to = 1850) {


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

  ref <- rbind(SuessR.reference.data, custom.data)


  # Screen input data to check for (and remove) missing values, with warning
  if(any(!(data$region %in% ref$region))) stop(paste("Unrecognized region: "), data$region[!(data$region %in% ref$region)])



  #ln(K0)
  ref$lnK0 <- with(ref,
                   -58.0931 + 90.5069 * (100/(sst + 273.15))
                   + (22.294 * log((sst + 273.15)/100))
                   + (S * (0.027766 + -0.025888 * ((sst+273.15)/100)
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

  data <- merge(data, ref, c("region", "year"))

  dat1850 <- ref[ref$year==1850,]
  dat1850$laws1850 <- with(dat1850, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat1850[,c("region", "laws1850")], "region")
  dat.correct.to <- ref[ref$year==correct.to,]
  dat.correct.to$laws.correct.to <- with(dat.correct.to, laws.fun(e1=esub1, e2=esub2, e.1=esubneg1,laws.CO2=CO2aq, P=P, sst=sst, C=C, b=beta))
  data <- merge(data, dat.correct.to[,c("region", "laws.correct.to")], "region")
  data$Laws.cor <- round(with(data, (laws.current - laws1850) - (laws.correct.to - laws1850)), 2)

  SuessR.out <- data[,c("id","year", "d13c", "Laws.cor")]
  names(SuessR.out)[3] <- "d13c.uncor"
  SuessR.out$Suess.cor  <- round(with(data, 0.014 * exp((year-1850)*0.027)
                                - 0.014*exp((correct.to-1850)*0.027)),2)
  SuessR.out$net.cor    <- SuessR.out$Suess.cor + SuessR.out$Laws.cor
  SuessR.out$d13c.cor   <- data$d13c + SuessR.out$net.cor
  SuessR.out <- SuessR.out[order(SuessR.out$id),]  # Sort by id if desired
  print(SuessR.out)
}




# Function for calculating regional uptake constants

reg.uptake <- function(year1, year2, d13C.change) {

  reg.up.const <- round(d13C.change/(exp((year2-1850)*0.027) - exp((year1 - 1850)*0.027)),2)

  print(reg.up.const)

}





