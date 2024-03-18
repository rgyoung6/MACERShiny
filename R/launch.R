# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
#********************************************Main program section***********************************************
#Roxygen2 Documentation:
#' @export
#'
#' @title Launch MACERShiny
#'
#' @author Robert G. Young
#'
#' @description
#' This function launches the MACERShiny Application
#'
#' @details
#' This function launches a MACERShiny Application which allows a user to run the DBTC
#' functions and process high throughput sequencing data.
#'
#' @examples
#' \dontrun{
#' launchMACERShiny()
#' }
#'
#' @returns
#' There are no values or files returned from this function
#'
#' @references
#' <https://github.com/rgyoung6/MACER>
#' Young, R. G., Hanner, R. H. (Submitted October 2023). Title Here. Biodiversity Data Journal.
#'
#' @note
#' This is a wrapper function which launches the MACERShiny package as a shiny
#' application in the systems default browser program
#'


# wrapper for shiny::shinyApp()
launchMACERShiny <- function() {
  shiny::shinyApp(ui = shinyAppUI, server = shinyAppServer, options = list(launch.browser = TRUE))
}
