% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concfuns.r
\name{sscformula}
\alias{barenblatt}
\alias{lane}
\alias{rouse}
\alias{sscformula}
\alias{tanaka}
\title{Suspended Sediment Vertical Gradient Formulas}
\usage{
rouse(z, H, b, w, k, us, pizzuto = FALSE)

lane(z, H, b, w, k, us)

barenblatt(z, H, b, w, k, us)

tanaka(z, H, b, w, k, us)
}
\arguments{
\item{z}{The depth at which to compute concentration.}

\item{H}{The maximum water depth.}

\item{b}{The active layer height.}

\item{w}{The settling velocity.}

\item{k}{The von Karmen constant.}

\item{us}{The shear velocity.}

\item{pizzuto}{Logical: apply the correction from Pizzuto (1984).}
}
\description{
Formulas for specifying the vertical gradient of suspended sediment.
}
\section{Functions}{
\itemize{
\item \code{rouse}: Rouse's formula

\item \code{lane}: Lane and Kalinske's formula

\item \code{barenblatt}: Barenblatt's formula

\item \code{tanaka}: Tanaka and Sugimoto's formula
}}

