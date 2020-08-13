#' Reference dataset for the SuessR Package
#'
#' A dataset containing the information used by SuessR to calculate regional Suess and Laws corrections
#' for \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C from marine organisms.
#'
#' @format A data frame with 680 rows and 8 variables:
#' \describe{
#'     \item{year}{calendar year, 1850-2019}
#'     \item{region}{geographic region, Aleutians, Bering, Gulf of Alaska, or Subpolar North Pacific}
#'     \item{r}{average cell radius of phytoplankton community, in microns}
#'     \item{sst}{yearly average sea surface temperature, in degrees C}
#'     \item{sss}{yearly average sea surface salinity, in practical salinity units}
#'     \item{CO2atm}{yearly average global atmospheric CO2 concentrations, in parts per million}
#'     \item{up.con}{regional uptake constant for Suess correction, unitless}
#'     \item{Cp}{Proportional change constant of oceanic/atmospheric CO\ifelse{html}{\out{<sub>2</sub>}}{\eqn{_2}} (regional), a proportion}
#'     }
#'
#' @source \url{https://www.ncdc.noaa.gov/data-access/marineocean-data/extended-reconstructed-sea-surface-temperature-ersst-v5}\cr
#'
#' \url{https://iridl.ldeo.columbia.edu/SOURCES/.CARTON-GIESE/.SODA/.v2p2p4/?Set-Language=en}\cr
#'
#' \url{https://coastwatch.noaa.gov/cw/satellite-data-products/sea-surface-salinity/miras-smos.html}\cr
#'
#' \url{https://scrippsco2.ucsd.edu/data/atmospheric_co2/icecore_merged_products.html}\cr
#'
#' @references Clark, C.T., M.R. Cape, M.D. Shapley, F.J. Mueter, B.P. Finney, and N. Misarti. (In Prep) SuessR: Regional Suess
#'   and Laws corrections for \ifelse{html}{\out{&delta;<sup>13</sup>}}{\eqn{{\delta}^{13}}}C from marine organisms.\cr
#'
#' Huang, B., Thorne, P. W., Banzon, V. F., Boyer, T., Chepurin, G., Lawrimore, J. H., … Zhang, H.-M. (2017).
#'   NOAA extended reconstructed sea surface temperature (ERSST), version 5. doi:10.7289/V5T72FNM.\cr
#'
#' Giese, B. S., & Ray, S. (2011). El Niño variability in simple ocean data assimilation (SODA), 1871-2008.
#'   Journal of Geophysical Research: Oceans, 116(2), 1–17. doi:10.1029/2010JC006695.\cr
#'
#' Keeling, C. D., Piper, S. C., Bacastow, R. B., Wahlen, M., Whorf, T. P., Heimann, M., & Meijer, H. A. (2005).
#'   Atmospheric CO2 and 13CO2 exchange with the terrestrial biosphere and oceans from 1978 to 2000: Observations
#'   and carbon cycle implications. In A history of atmospheric CO2 and its effects on plants, animals, and
#'   ecosystems (pp. 83–113). Springer.\cr
#'
#' MacFarling Meure, C., Etheridge, D., Trudinger, C., Steele, P., Langenfelds, R., Van Ommen, T., … Elkins, J. (2006).
#'   Law Dome CO2, CH4 and N2O ice core records extended to 2000 years BP. Geophysical Research Letters, 33(14),
#'   2000–2003. doi:10.1029/2006GL026152.\cr
#'
"SuessR.reference.data"
