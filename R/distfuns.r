#'Read Concentration Data
#'
#' @param f The file containing excess concentration data. Columns 
#'   are assumed to be grain size classes, in the order of 
#'   finest --> coarsest.
#' @return A dataframe.
#'
#' @import readr
#' @export
read_concentration = function(f) {
  d = read_csv(f)
  n = ncol(d)
  names(f) = c()[1:n]
  f


}


define_xs = function(C.width, ROB.height, ROB.width, LOB.height,
  LOB.width) {

  l = list(
    C.width = C.width,
    C.height = C.height,
    ROB.width = ROB.width,
    LOB.width = LOB.width
    )
  class(l) = "xs"
  l
}

#' Deposit 
distribute_conc = function(conc, H, xs, conc.fun) {
  ROB2LOB = xs$LOB.height - xs$ROB.height
  LOB2ROB = -ROB2LOB

  channel.only = c(0, min(xs$ROB, xs$LOB))
  channel.and.ROB = (ROB2LOB > 0) * c(channel.only, xs$LOB - xs$ROB)
  channel.and.LOB = (LOB2ROB > 0) * c(channel.only, xs$ROB - xs$LOB)
  channel.and.banks = c(max(channel.and.ROB, channel.and.LOB), H)

  section_conc = function(f, lims) integrate(f, lims[1], lims[2])$value

  conc.channel.only = section_conc(conc.fun, channel.only) *
     xs$channel.width
  conc.channel.and.ROB = section_conc(conc.fun, channel.and.ROB) *
    (xs$channel.width + xs$ROB.width)
  conc.channel.and.LOB = section_conc(conc.fun, channel.and.LOB) *
    (xs$channel.width + xs$LOB.width)
  conc.channel.and.banks = section_conc(conc.fun, channel.and.banks) *
    (xs$channel.width + xs$LOB.width + xs$ROB.width)

  conc.channel = (conc.channel.only + conc.channel.and.ROB +
  conc.channel.and.LOB + conc.channel.and.banks)

  
  
}


