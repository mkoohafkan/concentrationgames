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
#' Define a rectangular compound cross section.
#'
#' @param C.width The channel width.
#' ROB.height The right bank height.
#' ROB.width The right bank floodplain width.
#' LOB.height The left bank height.
#' LOB.width The left bank floodplain width.
#' @return an object of class 'xs'.
#'
#' @export
define_xs = function(C.height, C.width, ROB.height, ROB.width, LOB.height,
  LOB.width) {

  list(
    C.height = C.height,
    C.width = C.width,
    ROB.height = ROB.height,
    ROB.width = ROB.width,
    LOB.height = LOB.height,
    LOB.width = LOB.width
    )
}

#' Distribute Mass Over Cross Section
#'
#' Distribute sediment mass over a rectangular compound cross-section.
#'
#' @param mass The total mass of sediment to be distributed over the
#' cross-section.
#' @param xs The cross-section object.
#' @param H The water level relative to the channel invert elevation.
#' @param b The active layer height. 
#' @param conc.fun The function describing vertical concentration 
#'   gradient. If missing, a constant gradient is assumed.
#' @return The cross-section object with mass distribution.
#'
#' @export
distribute_mass = function(mass, xs, H, b, conc.fun = NULL, ...) {

  # divide cross section into vertical segments
  C.lims = c(b, H)
  ROB.lims = c(xs$ROB.height, H)
  LOB.lims = c(xs$LOB.height, H)

  C.area = xs$C.width * (C.lims[2] - C.lims[1])
  ROB.area = xs$ROB.width * (ROB.lims[2] - ROB.lims[1])
  LOB.area = xs$LOB.width * (LOB.lims[2] - LOB.lims[1])

  area = C.area + ROB.area + LOB.area

  # total average concentration
  c.avg = mass / area

  # concentration function and cb
  if (is.null(conc.fun)) {
    cb = c.avg
    conc.fun = function(z, ...)
      rep(1, length(z))
    } else {
      cb = get_cb(c.avg, conc.fun, H, b, ...)
    }

  # function to compute average concentration in a segment
  segment_conc = function(f, cb, lims, ...) {
    if (lims[2] != lims[1])
      cb * integrate(f, lims[1], lims[2], ...)$value
    else
      0
  }
  C.mass = segment_conc(conc.fun, cb, C.lims, H = H, b = b, ...) *
    xs$C.width
  ROB.mass = segment_conc(conc.fun, cb, ROB.lims, H = H, b = b, ...) *
    xs$ROB.width
  LOB.mass = segment_conc(conc.fun, cb, LOB.lims, H = H, b = b, ...) *
    xs$LOB.width

  # correct masses to add up to total mass
  int.mass = C.mass + ROB.mass + LOB.mass
  mf = mass / int.mass

  # return results
  c(xs, list(
    C.mass = mf * C.mass,
    ROB.mass = mf * ROB.mass,
    LOB.mass = mf * LOB.mass
    ))

}

#' Overbank Deposion
#' 
#' Compute overbank depositon based on an exponential decay function.
#'
#' @param ob.width The overbank width.
#' @param ob.nodes 
ob_dist = function(ob.width, ob.nodes, decay) {
  
  parint = 1 - exp( - decay * ob.width)
  if(parint == 0)
    rep(1/ob.width, length(ob.nodes))
  else
    decay * exp( - decay * ob.nodes) / parint  
}

#' Mass Deposition
#'
#' Deposit mass using a decay function.
#'
#' @param xs The cross section object with mass data.
#' @param dy Spacing between cross section nodes.
#' @param decay The decay coefficient.
#' @param truncate.dist The maximum distance from the channel to deposit
#'   on the banks.
#' @return A dataframe of node positions and mass deposited.
#'
#' @export
deposit_mass = function(xs, dy, decay, truncate.dist) {

  # channel
  C.nodes = seq(0, xs$C.width, by = dy)
  m.per.node.C = rep(xs$C.mass / xs$C.width, length(C.nodes))
 
   # OB
  if (missing(truncate.dist)) {
    ROB.dist = xs$ROB.width
    LOB.dist = xs$LOB.width
  } else {
    ROB.dist = min(truncate.dist, xs$ROB.width)
    LOB.dist = min(truncate.dist, xs$LOB.width)
  }
  ROB.nodes = seq(0, ROB.dist, by = dy)
  LOB.nodes = seq(0, LOB.dist, by = dy)
  m.per.node.ROB = xs$ROB.mass * ob_dist(xs$ROB.width, ROB.nodes, 
    decay)
  m.per.node.LOB = xs$LOB.mass * ob_dist(xs$LOB.width, LOB.nodes, 
    decay)
  
  ROB.all = seq(0, xs$ROB.width, by = dy)
  LOB.all = seq(0, xs$LOB.width, by = dy)
  
  ROB.all.m = rep(0, length(ROB.all))
  ROB.all.m[1:length(m.per.node.ROB)] = m.per.node.ROB
  LOB.all.m = rep(0, length(LOB.all))
  LOB.all.m[1:length(m.per.node.LOB)] = m.per.node.LOB
  
  #print(trapz(ROB.nodes, m.per.node.ROB) - xs$ROB.mass)

  data.frame(
    label = c(
      rep("LOB", length(LOB.all)), 
      rep("C", length(C.nodes)),
      rep("ROB", length(ROB.all))
    ), 
    y = c(
      LOB.all, 
      C.nodes + xs$LOB.width, 
      ROB.all + xs$LOB.width + xs$C.width
    ), 
    z = c(
      rep(xs$LOB.height, length(LOB.all)), 
      rep(xs$C.height, length(C.nodes)), 
      rep(xs$ROB.height, length(ROB.all))
    ),
    m = c(
      rev(LOB.all.m), 
      m.per.node.C, 
      ROB.all.m),
    dy = dy
  )
}

mass_to_bedchange = function(nd, p = 0.4, d = 2650){  
  f = 1/(d*p)  
  nd["dz"] = nd$m*f
  nd
}

