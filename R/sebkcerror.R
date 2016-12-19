#' Internal function to Catch  errors and warnings
#'
#' @param expr an  R expression to evaluate
#'
#' @return value a list with 'value' and 'warning', where value may be an error caught
#' @export
#'
#' @examples
#' thiserror=str( sebkc.tryCatch( log( "a" ) ) )$value
sebkc.tryCatch = function(expr){
  
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}