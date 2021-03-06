% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distfuns.r
\name{read_concentration}
\alias{read_concentration}
\title{Read Concentration Data}
\usage{
read_concentration(f, classes)
}
\arguments{
\item{f}{The file containing excess concentration data. Columns 
are assumed to be grain size classes, in the order of 
finest --> coarsest.}

\item{classes}{The grain class names. If missing, the standard 20
classes is assumed.}
}
\value{
A dataframe.
}
\description{
Read Concentration Data
}

