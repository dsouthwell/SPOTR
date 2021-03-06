% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/9_Functions_test_package_4Methods_random_XYsites.R
\name{fit.occ.3method}
\alias{fit.occ.3method}
\title{Fit occupancy model with three detection methods}
\usage{
fit.occ.3method(method1, method2, method3, repeats1, repeats2, repeats3,
  s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park,
  fail.park, value, ss, two.tailed, n.park)
}
\arguments{
\item{method1}{The first detection method relevant to species ss}

\item{method2}{The second detection method relevant to species ss}

\item{method3}{The third detection method relevant to species ss}

\item{repeats1}{The number of repeat visits for the first detection method}

\item{repeats2}{The number of repeat visits for the second detection method}

\item{repeats3}{The number of repeat visits for the third detection method}

\item{s.years}{A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)}

\item{n.sites}{The number of sites monitored}

\item{xy.sites}{The XY coordinates of monitored sites}

\item{park.ID}{A vector specifying that location of each site with sub-level parks}

\item{park.level}{Set to TRUE is power is estimated within regional level management unit, FALSE otherwise}

\item{powcnt}{A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape}

\item{fail}{A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a landscape level in unmarked}

\item{pow.park}{A vector that keeps track of how many times a significant trend in occupancy is detected within each park}

\item{fail.park}{A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a park level in unmarked}

\item{value}{The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test}

\item{ss}{An index to loop through each species}

\item{two.tailed}{Set to TRUE if conducting a two-tailed test, FALSE otherwise}

\item{n.park}{The number of parks in which to estimate power}
}
\description{
This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected using 3 detection methods
A trend in occupancy is estimated and confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from each regional level management unit.
}
