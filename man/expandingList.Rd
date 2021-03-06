% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{expandingList}
\alias{expandingList}
\title{Utility function for creating and adding to a list that dynamically resizes itself as items are added.
This should hopefully make the storage of results as the simulation is running more efficient.}
\usage{
expandingList(capacity = 10)
}
\arguments{
\item{capacity}{Starting size of the initial list. Everytime capacity is reached the list is doubled in size.}
}
\description{
Utility function for creating and adding to a list that dynamically resizes itself as items are added.
This should hopefully make the storage of results as the simulation is running more efficient.
}
