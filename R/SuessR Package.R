# SuessR Package


# Function to calculate and apply Suess and Laws corrections.

SuessR <- function(data, correct.to = 1850) {

  # Screen input data to check for (and remove) missing values, with warning
  if(anyNA(data)) warning("Rows with missing values were deleted")
  data <- na.omit(data)


  # Screen input data to make sure all years provided are integers between 1850 and the current year -1
  current.year <- as.numeric(substr(date(), 21,24))-1
  if(any(!(data$year %in% 1850:current.year)))  stop(paste("Year must be an integer between 1850 and", current.year))


  # Create output data frame
  SuessR.out <- data.frame(id = data$id, year = data$year)

  ref.array <- simplify2array(by(SuessR.reference.data[,-2], SuessR.reference.data$region, as.matrix))

  #ln(K0)
  lnK0 <- -58.0931+90.5069*(100/(SuessR.reference.data$sst+273.15))+(22.294*log((SuessR.reference.data$sst+273.15)/100))+(SuessR.reference.data$S*(0.027766+-0.025888*((SuessR.reference.data$sst+273.15)/100)+0.005058*((SuessR.reference.data$sst+273.15)/100)^2))
  laws <- data.frame(lnK0 = lnK0)
  laws$region <- SuessR.reference.data$region

  #Ocean Increase
  for(j in 2:length(SuessR.reference.data$sst)) {
    laws$oi[j] <- (SuessR.reference.data$CO2atm[j] - SuessR.reference.data$CO2atm[j-1])*0.4
  }

  #~f(CO2)ocean
  for(k in 1:length(SuessR.reference.data$sst)) {
    laws$fCO2[k] <- sum(laws$oi[1:k], na.rm = T) + 284.25
  }

  #CO2aq
  for(l in 1:length(SuessR.reference.data$sst)) {
    laws$CO2aq[l] <- exp(laws$lnK0[l])*laws$fCO2[l]
  }

  laws.array <- simplify2array(by(laws[,-2], laws$region, as.matrix))

  for(i in 1:length(data$year)) {

    SuessR.out$d13c.uncor[i] <- data$d13c[data$year == data$year[i] & data$region == data$region[i]]
    SuessR.out$suess.cor[i]  <- 0.014*exp((data$year[data$year == data$year[i] & data$region == data$region[i]]-1850)*0.027) - 0.014*exp((correct.to-1850)*0.027)
    SuessR.out$laws.cor[i]   <- ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]] +
                                    ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]] -
                                    ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]] -
                            (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == data$year[i]]*
                                    ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]])/
                                   (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]*
                                    ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]*
                                 (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]])))))*
                                  ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]-
                                    ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]])/
                                   (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]+1))) -
                                   ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] +
                                     ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
                                     ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850]) -
                             (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == 1850]*
                                     ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
                                    (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                                     ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                                  (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850])))))*
                                   ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
                                     ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
                                    (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850]+1)))) -

                                   ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] +
                                     ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] -
                                     ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] -
                             (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == correct.to]*
                                     ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == correct.to])/
                                    (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == correct.to]*
                                     ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == correct.to]*
                                  (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == correct.to])))))*
                                   ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] -
                                     ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == correct.to])/
                                    (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == correct.to]+1))) -
                                   ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] +
                                     ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
                                     ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850]) -
                             (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == 1850]*
                                     ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
                                    (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                                     ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                                  (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850])))))*
                                   ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
                                     ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
                                    (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850]+1))))
    SuessR.out$net.cor[i]    <- SuessR.out$Suess.cor[i] + SuessR.out$laws.cor[i]
    SuessR.out$d13c.cor[i]   <- data$d13c[data$year == data$year[i] & data$region == data$region[i]] + SuessR.out$net.cor[i]
  }
  print(SuessR.out)
}





# Second function that allows users to add a custom region

SuessR.custom <- function(data, custom.data, correct.to = 1850) {


  # Screen input data to check for (and remove) missing values, with warning
  if(anyNA(data)) warning("Rows with missing values were deleted")
  data <- na.omit(data)


  # Screen input data to make sure all years provided are integers between 1850 and the current year -1
  current.year <- as.numeric(substr(date(), 21,24))-1
  if(any(!(data$year %in% 1850:current.year)))  stop(paste("Year must be an integer between 1850 and", current.year))


  # Create output data frame
     SuessR.out <- data.frame(id = data$id, year = data$year)

  SuessR.reference.data.custom <- as.data.frame(rbind(SuessR.reference.data, custom.data))

  ref.array <- simplify2array(by(SuessR.reference.data.custom[,-2], SuessR.reference.data.custom$region, as.matrix))

  #ln(K0)
  lnK0 <- -58.0931+90.5069*(100/(SuessR.reference.data.custom$sst+273.15))+(22.294*log((SuessR.reference.data.custom$sst+273.15)/100))+(SuessR.reference.data.custom$S*(0.027766+-0.025888*((SuessR.reference.data.custom$sst+273.15)/100)+0.005058*((SuessR.reference.data.custom$sst+273.15)/100)^2))
  laws <- data.frame(lnK0 = lnK0)
  laws$region <- SuessR.reference.data.custom$region

  #Ocean Increase
  for(j in 2:length(SuessR.reference.data.custom$sst)) {
    laws$oi[j] <- (SuessR.reference.data.custom$CO2atm[j] - SuessR.reference.data.custom$CO2atm[j-1])*0.4
  }

  #~f(CO2)ocean
  for(k in 1:length(SuessR.reference.data.custom$sst)) {
    laws$fCO2[k] <- sum(laws$oi[1:k], na.rm = T) + 284.25
  }

  #CO2aq
  for(l in 1:length(SuessR.reference.data.custom$sst)) {
    laws$CO2aq[l] <- exp(laws$lnK0[l])*laws$fCO2[l]
  }

  laws.array <- simplify2array(by(laws[,-2], laws$region, as.matrix))

  for(i in 1:length(data$year)) {

    SuessR.out$d13c.uncor[i] <- data$d13c[data$year == data$year[i] & data$region == data$region[i]]
    SuessR.out$suess.cor[i]  <- 0.014*exp((data$year[data$year == data$year[i] & data$region == data$region[i]]-1850)*0.027) - 0.014*exp((correct.to-1850)*0.027)
    SuessR.out$laws.cor[i]   <- ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]] +
                                  ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]] -
                                  ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]] -
                          (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == data$year[i]]*
                                  ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]])/
                                 (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]*
                                                ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]*
                                                (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]])))))*
                                    ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]-
                                        ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]])/
                                       (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == data$year[i]]+1))) -
                                   ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] +
                                       ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
                                       ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850]) -
                                      (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == 1850]*
                                                ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
                                               (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                                                  ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                                                  (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850])))))*
                                      ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
                                          ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
                                         (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850]+1)))) -

      ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] +
          ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] -
          ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] -
          (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == correct.to]*
                    ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == correct.to])/
                   (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == correct.to]*
                      ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == correct.to]*
                      (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == correct.to])))))*
          ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == correct.to] -
              ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == correct.to])/
             (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == correct.to]+1))) -
         ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] +
             ref.array[,"esub1",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
             ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850]) -
            (1/(1+((laws.array[,"CO2aq",as.character(data$region[i])][ref.array[,"year",data$region[i]] == 1850]*
                      ref.array[,"P",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
                     (ref.array[,"mu",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                        ref.array[,"C",data$region[i]][ref.array[,"year",data$region[i]] == 1850]*
                        (1+ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850])))))*
            ((ref.array[,"esub2",data$region[i]][ref.array[,"year",data$region[i]] == 1850] -
                ref.array[,"esubneg1",data$region[i]][ref.array[,"year",data$region[i]] == 1850])/
               (ref.array[,"beta",data$region[i]][ref.array[,"year",data$region[i]] == 1850]+1))))
    SuessR.out$net.cor[i]    <- SuessR.out$Suess.cor[i] + SuessR.out$laws.cor[i]
    SuessR.out$d13c.cor[i]   <- data$d13c[data$year == data$year[i] & data$region == data$region[i]] + SuessR.out$net.cor[i]
  }
  print(SuessR.out)
}
