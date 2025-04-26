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
#'   \item{Cases}{Cumulative cases of Covid-19 from 10th March 2020}
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



#' Cumulative cases of Covid-19 in England.
#'
#' @docType data
#'
#' @usage data(england)
#'
#' @format An object of class `"xts"`;
#' \describe{
#'   \item{cum_cases}{Cumulative cases of Covid-19}
#'   \item{cum_admissions}{Cumulative hospital admissions due to Covid-19}
#'   \item{cum_deaths}{Cumulative deaths due to Covid-19}
#'   \item{hospital_cases}{???}
#' }
#'
#' @keywords datasets
#'
#' @references Downloaded from https://ukhsa-dashboard.data.gov.uk/topics/covid-19
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
#' # plot daily cases
#' plot(diff(gauteng_weather_2021))
"gauteng_weather_2021"

#' Weather data of England.
#'
#' @docType data
#'
#' @usage data(england_weather_2021)
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
#' data(england_weather_2021)
#' # plot daily cases
#' plot(diff(gauteng_weather_2021))
"england_weather_2021"
