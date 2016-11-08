#'Read Concentration Data
#'
#' @param f The file containing excess concentration data. Columns 
#'   are assumed to be grain size classes, in the order of 
#'   finest --> coarsest.
#' @param classes The grain class names. If missing, the standard 20
#'   classes is assumed.
#' @return A dataframe.
#'
#' @import readr
#' @export
read_concentration = function(f, classes) {
  if (missing(classes))
    classes = c('Clay', 'VFM', 'FM', 'MM', 'CM', 'VFS', 'FS', 'MS', 
      'CS', 'VCS', 'VFG', 'FG', 'MG', 'CG', 'VCG', 'SC', 'LC', 'SB',
      'MB', 'LB')
  setNames(read_csv(f), classes)

}

#' Cross Section
#' 
#' Define a cross-section.
#'
#' @param C.width The channel width.
#' ROB.height The right bank height.
#' ROB.width The right bank floodplain width.
#' LOB.height The left bank height.
#' LOB.width The left bank floodplain width.
#' @return an object of class 'xs'.
#'
#' @export
define_xs = function(C.width, ROB.height, ROB.width, LOB.height,
  LOB.width) {

  structure(list(
    C.width = C.width,
    ROB.height= ROB.height,
    ROB.width = ROB.width,
    LOB.height = LOB.height,
    LOB.width = LOB.width
    ), class = "xs")
}

#' Distribute Mass Over Cross Section
#'
#' Distribute sediment mass over a rectangular compound cross-section.
#'
#' @param mass The total mass of sediment to be distributed over the
#' cross-section.
#' @param H The water level. The channel invert elevation is assumed
#'   to be zero.
#' @param xs The cross-section object.
#' @param conc.fun The function describing vertical concentration 
#'   gradient. If missing, a constant gradient is assumed.
#' @return The cross-section object with mass distribution.
#'
#' @export
distribute_mass = function(mass, H, xs, conc.fun) {

  ROB2LOB = xs$LOB.height - xs$ROB.height
  LOB2ROB = -ROB2LOB
  channel.only = c(0, min(xs$ROB.height, xs$LOB.height))
  channel.and.ROB = (ROB2LOB > 0) * c(max(channel.only), 
    max(xs$LOB.height, xs$ROB.height))
  channel.and.LOB = (LOB2ROB > 0) * c(max(channel.only), 
    max(xs$ROB.height, xs$LOB.height))
  channel.and.banks = c(max(xs$ROB.height, xs$LOB.height), H)

  vol.C = H * xs$C.width
  vol.ROB = (H - xs$ROB.height) * xs$ROB.width
  vol.LOB = (H - xs$LOB.height) * xs$LOB.width

  vol.total = vol.C + vol.ROB + vol.LOB

  if (missing(conc.fun))
    conc.fun = function(x)
      rep(mass / vol.total, length(x))

  section_conc = function(f, lims, ...)
    integrate(f, lims[1], lims[2], ...)$value


  mass.C.only = section_conc(conc.fun, channel.only) * xs$C.width
  mass.C.and.ROB = section_conc(conc.fun, channel.and.ROB) *
    (xs$C.width + xs$ROB.width)
  mass.C.and.LOB = section_conc(conc.fun, channel.and.LOB) *
    (xs$C.width + xs$LOB.width)
  mass.C.and.banks = section_conc(conc.fun, channel.and.banks) *
    (xs$C.width + xs$LOB.width + xs$ROB.width)

  c(xs, list(
    C.mass = mass.C.only + mass.C.and.ROB * xs$C.width / 
      (xs$C.width + xs$ROB.width) + mass.C.and.LOB * xs$C.width / 
      (xs$C.width + xs$LOB.width) + mass.C.and.banks * xs$C.width / 
      (xs$C.width + xs$ROB.width + xs$LOB.width),
    ROB.mass = mass.C.and.ROB * xs$ROB.width / (xs$C.width + 
      xs$ROB.width) + mass.C.and.banks * xs$ROB.width / 
      (xs$C.width + xs$ROB.width + xs$LOB.width),
    LOB.mass = mass.C.and.LOB * xs$LOB.width / (xs$C.width + 
      xs$LOB.width) + mass.C.and.banks * xs$LOB.width / 
      (xs$C.width + xs$ROB.width + xs$LOB.width)
  ))
}


