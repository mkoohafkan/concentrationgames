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
distribute_mass = function(mass, xs, H, b, conc.fun, ...) {

  # divide cross section into vertical segments
  ROB2LOB = xs$LOB.height - xs$ROB.height
  LOB2ROB = -ROB2LOB
  channel.only = c(b, min(xs$ROB.height, xs$LOB.height))
  channel.and.ROB = (ROB2LOB > 0) * c(max(channel.only),
    max(xs$LOB.height, xs$ROB.height))
  channel.and.LOB = (LOB2ROB > 0) * c(max(channel.only),
    max(xs$ROB.height, xs$LOB.height))
  channel.and.banks = c(max(xs$ROB.height, xs$LOB.height), H)

  # compute segment volumes
  vol.C.only = (channel.only[2] - channel.only[1]) * xs$C.width
  vol.C.and.ROB = (channel.and.ROB[2] - channel.and.ROB[1]) *
    (xs$C.width + xs$ROB.width)
  vol.C.and.LOB = (channel.and.LOB[2] - channel.and.LOB[1]) *
    (xs$C.width + xs$LOB.width)
  vol.C.and.banks = (channel.and.banks[2] - channel.and.banks[1]) *
    (xs$C.width + xs$LOB.width + xs$ROB.width)

  vol.C = (H - b) * xs$C.width
  vol.ROB = (H - xs$ROB.height) * xs$ROB.width
  vol.LOB = (H - xs$LOB.height) * xs$LOB.width

  vol.total = vol.C + vol.ROB + vol.LOB

  # total average concentration
  c.avg = mass / vol.total
  # check concentration function and cb
  if (missing(conc.fun)) {
    cb = c.avg
    conc.fun = function(z, ...)
      rep(1, length(z))
    } else {
      cb = get_cb(c.avg, conc.fun, H, b, ...)
    }
  # function to compute average concentration in a segment
  segment_conc = function(f, cb, lims, ...) {
    if (lims[2] != lims[1])
      cb * integrate(f, lims[1], lims[2], ...)$value / (lims[2] - lims[1])
    else
      0
  }
  # compute mass in each segment
  mass.C.only = segment_conc(conc.fun, cb, channel.only, H = H, b = b,
    ...) * vol.C.only
  mass.C.and.ROB = segment_conc(conc.fun, cb, channel.and.ROB, H = H,
    b = b, ...) * vol.C.and.ROB
  mass.C.and.LOB = segment_conc(conc.fun, cb, channel.and.LOB, H = H,
    b = b, ...) * vol.C.and.LOB
  mass.C.and.banks = segment_conc(conc.fun, cb, channel.and.banks,
    H = H, b = b, ...) * vol.C.and.banks

  # aggregate segments into Channel, ROB and LOB
  frac.C.ROB = xs$C.width / (xs$C.width + xs$ROB.width)
  frac.C.LOB = xs$C.width / (xs$C.width + xs$LOB.width)
  frac.C.banks = xs$C.width / (xs$C.width + xs$ROB.width + xs$LOB.width)
  frac.ROB.C = xs$ROB.width / (xs$C.width + xs$ROB.width)
  frac.ROB.banks = xs$ROB.width / (xs$C.width + xs$ROB.width + xs$LOB.width)
  frac.LOB.C = xs$LOB.width / (xs$C.width + xs$LOB.width)
  frac.LOB.banks = xs$LOB.width / (xs$C.width + xs$ROB.width + xs$LOB.width)

  C.mass = mass.C.only + mass.C.and.ROB * frac.C.ROB +
    mass.C.and.LOB * frac.C.LOB + mass.C.and.banks * frac.C.banks
  ROB.mass = mass.C.and.ROB * frac.ROB.C +
    mass.C.and.banks * frac.ROB.banks
  LOB.mass = mass.C.and.LOB * frac.LOB.C +
    mass.C.and.banks * frac.LOB.banks

  # correct masses to add up to total mass
  int.mass = C.mass + ROB.mass + LOB.mass
  mass.factor = mass / int.mass

  # return results
  c(xs, list(
    C.mass = mass.factor * C.mass,
    ROB.mass = mass.factor * ROB.mass,
    LOB.mass = mass.factor * LOB.mass
    ))

}


ob_dist = function(ob.width, ob.nodes, decay) {
  parint = function(Y, a = decay)
    1 - exp( - a * Y)
  decayfun = function(y, a = decay)
    a * exp( - a * y)
  decayfun(ob.nodes) / parint(ob.width)
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
  
  #trapz(ROB.nodes, m.per.node.ROB)

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
  dy = unique(nd$dy)
  
  f = dy/(d*p)  
  nd["dz"] = nd$m*f
  nd
}

