# Created by: Craig Thamotheram
# Created on: 15/02/2022

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Cumulative cases of Covid-19 in the South African province of Gauteng.
#'
#' @docType data
#'
#' @usage data(gauteng)
#'
#' @format An object of class `"xts"`;
#' \describe{
#'   \item{cum_cases}{Cumulative cases of Covid-19 from 10th March 2020}
#' }
#'
#' @keywords datasets
#'
#' @references Downloaded from https://sacoronavirus.co.za/
#'
#' @examples
#' data(gauteng)
#' # plot daily cases
#' plot(diff(gauteng))
"gauteng"



#' Covid-19 metrics for England.
#'
#' @docType data
#'
#' @usage data(england)
#'
#' @format An object of class `"xts"`;
#' \describe{
#'   \item{cum_cases}{Cumulative cases of Covid-19}
#'   \item{cum_admissions}{Cumulative hospital admissions with Covid-19 since the start of the pandemic.}
#'   \item{cum_deaths}{Cumulative deaths within 28 days of a positive test for Covid-19 by death date since the start of the pandemic.}
#'   \item{hospital_cases}{Number of patients in hospital with confirmed Covid-19 each day. Note this is NOT cumulative.}
#' }
#'
#' @keywords datasets
#'
#' @references Downloaded from https://ukhsa-dashboard.data.gov.uk/topics/covid-19. Now at https://ukhsa-dashboard.data.gov.uk/covid-19-archive-data-download. Full documentation from UKHSA included in downloads from the website.  
#'
#' @examples
#' data(england)
#' # plot daily cases
#' plot(diff(england))
"england"




#' Cumulative cases of Covid-19 in Italy and UK, before 14 Dec 2020.
#'
#' @docType data
#'
#' @usage data(ukitaly)
#'
#' @format An object of class `"xts"`;
#' \describe{
#'   \item{UK}{Cumulative cases of Covid-19 in the UK}
#'   \item{Italy}{Cumulative cases of Covid-19 in Italy}
#' }
#'
#' @keywords datasets
#'
#' @references Downloaded from https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide
#'
#' @examples
#' data(ukitaly)
#' # plot daily cases
#' plot(diff(ukitaly))
"ukitaly"

#' Weather data of the South African province of Gauteng.
#'
#' @docType data
#'
#' @usage data(gauteng_weather_2021)
#'
#' @format An object of class `"xts"`;
#' \describe{
#'   \item{TempC}{Mean air temperature in degree Celsius at 2 meters}
#'   \item{TempC_max}{Maximum Temperature in degree Celsius at 2 meters}
#'   \item{RelHumid}{Relative humidity (%)}
#'   \item{WindSpeed}{Wind speed in m/s at 10 meters }
#' }
#'
#' @keywords datasets
#'
#' @references Downloaded from https://power.larc.nasa.gov/data-access-viewer/
#'
#' @examples
#' data(gauteng_weather_2021)
"gauteng_weather_2021"

#' Weather data of England.
#'
#' @docType data
#'
#' @usage data(england_weather_2021)
#'
#' @format An object of class `"xts"`;
#' \describe{
#'   \item{temperature_C}{Daily mean air temperature in degree Celsius (C) at 2 metres}
#'   \item{max_temp_C}{Daily maximum temperature in degree Celsius (C) at 2 metres}
#'   \item{relhum_percnt}{Relative humidity at 2 metres (%)}
#'   \item{windspd_mtrs_p_sec}{Wind speed in metres per second (m/s) at 10 metres}
#'   \item{precip_mtrs}{Preciptation (millimetres per day; mm/day)}
#' }
#'
#' @keywords datasets
#'
#' @references Downloaded from https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview and aggregated across the latitude and longitude corners for England by the authors.
#'
#' @examples
#' data(england_weather_2021)
#' 
"england_weather_2021"
