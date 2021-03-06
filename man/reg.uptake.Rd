% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SuessR_Package.R
\encoding{UTF-8}
\name{reg.uptake}
\alias{reg.uptake}
\title{Calculate a regional carbon dioxide uptake constant from empirical DIC data}
\usage{
reg.uptake(year1, year2, d13c.change)
}
\arguments{
\item{year1}{The year in which samples were collected for the first set of DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C \cr observations.}

\item{year2}{The year in which samples were collected for the second set of DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C \cr observations.}

\item{d13c.change}{The observed change in DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C between the two years.}
}
\value{
Returns a numerical value representing the regional uptake constant. This value can be supplied as part of the
  'custom.region.data' (filling a column titled 'up.con') argument for the SuessR.custom() function.
}
\description{
The reg.uptake() function calculates the regional uptake constant required to modify the global Suess effect
  curve to be region-specific.
}
\details{
This function calculates the regional uptake constant used to modify the global Suess effect curve to make it
  specific to a region. This method requires \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C values for DIC from the area of interest from two different time points,
  ideally separated by at least a decade. The function uses the magnitude (i.e., absolute value) of the observed change in
  DIC \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C ('d13c.change'), the year of the first observation ('year1'), and the year of the second observation ('year2')
  to calculate the regional uptake constant. This value can then be supplied as part of the 'custom.region.data' argument
  (filling a column titled 'up.con') for the SuessR.custom() function.
}
\examples{
year1 <- 1970
year2 <- 1980
d13c.change <- 0.3

reg.uptake(year1 = year1, year2 = year2, d13c.change = d13c.change)

}
