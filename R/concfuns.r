#' Suspended Sediment Vertical Gradient Formulas
#'
#' Formulas for specifying the vertical gradient of suspended sediment.
#'
#' @param z The depth at which to compute concentration.
#' @param H The maximum water depth.
#' @param crb The ratio of water depth to active layer depth \code{b/H}.
#' @param w The settling velocity.
#' @param k The von Karmen constant.
#' param us The shear velocity.
#' @name sscformula
NULL

#' @describeIn sscformula Rouse's formula
#' @export
rouse = function(z, H, crb, w, k, us)
  (((H - z) / z) * (crb / (1 - crb))) ^ (w / (k * us))

#' @describeIn sscformula Lane and Kalinske's formula
#' @export
lane = function(z, H, crb, w, k, us)
  exp(-15 * w * (z - crb * H) / (us * H))

#' @describeIn sscformula Barenblatt's formula
#' @export
barennblatt = function(z, H, crb, w, k, us)
  ((sqrt(1 - z / H) / sqrt(1 - crb)) * ((1 - sqrt(1 - crb)) / 
    (1 - sqrt(1 - z / H)))) ^ (w / (k * us))

#' @describeIn sscformula Tanaka and Sugimoto's formula
#' @export
tanaka = function(z, H, crb, w, k, us)
  (((1 + (sqrt(1 - z / H))) / (1 + sqrt(1 - crb))) * 
    ((1 - sqrt(1 - crb)) / (1 - sqrt(1 - z / H)))) ^ (w / (k * us))
