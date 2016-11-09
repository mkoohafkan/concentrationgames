#' Suspended Sediment Vertical Gradient Formulas
#'
#' Formulas for specifying the vertical gradient of suspended sediment.
#'
#' @param z The depth at which to compute concentration.
#' @param H The maximum water depth.
#' @param b The active layer height.
#' @param w The settling velocity.
#' @param k The von Karmen constant.
#' @param us The shear velocity.
#' @param pizzuto Logical: apply the correction from Pizzuto (1984).
#' @name sscformula
NULL

#' @describeIn sscformula Rouse's formula
#' @export
rouse = function(z, H, b, w, k, us, pizzuto = FALSE) {
  R = w / (k * us)
  if (pizzuto)
    R = 0.740 +0.362*log(R)
  (((H - z) / z) * (b / (H - b))) ^ R
} 

#' @describeIn sscformula Lane and Kalinske's formula
#' @export
lane = function(z, H, b, w, k, us)
  exp(-15 * w * (z - b) / (us * H))

#' @describeIn sscformula Barenblatt's formula
#' @export 
barennblatt = function(z, H, b, w, k, us)
  ((sqrt(1 - z / H) / sqrt(1 - b / H)) * ((1 - sqrt(1 - b / H)) /
    (1 - sqrt(1 - z / H)))) ^ (w / (k * us))

#' @describeIn sscformula Tanaka and Sugimoto's formula
#' @export
tanaka = function(z, H, b, w, k, us)
  (((1 + (sqrt(1 - z / H))) / (1 + sqrt(1 - b / H))) *
    ((1 - sqrt(1 - b / H)) / (1 - sqrt(1 - z / H)))) ^ (w / (k * us))

#' Active Layer Concentration
#'
#' Get the active layer concentration at height \code{b} above the bed.
#' 
#' @param conc.avg The cross-section average concentration.
#' @param fun The vertical gradient function to use.
#' @inheritParams sscformula
#' @param ... Other arguments to pass to \code{fun}.
#' @return The concentration at the active layer boundary.
#'
#' @export 
get_cb = function(conc.avg, fun, H, b, ...) {
  conc.avg * (H- b) / integrate(fun, lower = b, upper = H, 
    H = H, b = b, ...)$value
}
