

# S4 classes for motifActivity



# ---------------------------------------------------------------------------- #
# --------------------- #
#  Classes for identification and processing of regulation sites
# --------------------- #



#' An S4 class for storing \code{getMotifHits} function results
#'
#' The resulting object stores the results of \code{getMotifHits} function
#'
#' @section Slots:\describe{
#'                  \item{\code{matrix}: XYZ}
#'                 }
#'
#' @name MotifHitsMatrix-class
#' @rdname MotifHitsMatrix-class
#' @seealso \code{\link{getMotifHits}}
#' @export
setClass("MotifHitsMatrix",contains = "matrix")



#' An S4 class for storing \code{getSignal} function results
#'
#' The resulting object stores the results of \code{getSignal} function
#'
#' @section Slots:\describe{
#'                  \item{\code{matrix}: XYZ},
#'          \item{\code{signal.type}: XYZ}
#'                 }
#'
#' @name SignalMatrix-class
#' @rdname SignalMatrix-class
#' @seealso \code{\link{getSignalMatrix}}
#' @export
setClass("SignalMatrix",contains = "matrix")



# ---------------------------------------------------------------------------- #
# --------------------- #
#  Classes for modelling
# --------------------- #


#' An S4 class for storing \code{getModel} function results
#'
#' The resulting object stores the results of \code{getModel} function
#'
#' @section Slots:\describe{
#'                  \item{\code{signal}: XYZ},
#'          \item{\code{motif.hits}: XYZ}
#'                 }
#'
#' @name Model-class
#' @rdname Model-class
#' @seealso \code{\link{getModel}}
#' @export
setClass("Model",
         representation(
           norm.signal = "SignalMatrix",
           norm.motif.hits = "MotifHitsMatrix",
           coef = "data.frame",
           pval = "data.frame"
         ))






