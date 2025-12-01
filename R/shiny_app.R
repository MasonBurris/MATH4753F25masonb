#' Title
#'
#' @returns Webpage
#' @importFrom shiny runApp
#' @export
#'
#' @examples
#' \dontrun{shinymle()}
shinymle <- function(){
  shiny::runApp(system.file("SHINY", package = "MATH4753F25masonb"),
                launch.browser = TRUE)
}
