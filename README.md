# nichefillr
R package implementing niche filling simulations.

<!-- badges: start -->
  [![Build Status](https://travis-ci.org/rdinnager/nichefillr.svg?branch=development)](https://travis-ci.org/rdinnager/nichefillr)
  
  [![Codecov test coverage](https://codecov.io/gh/rdinnager/nichefillr/branch/master/graph/badge.svg)](https://codecov.io/gh/rdinnager/nichefillr?branch=master)
  
  [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
  
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1166254.svg)](https://doi.org/10.5281/zenodo.1166254)
  <!-- badges: end -->

This is a very early alpha release. Currently there is little in the way of error-checking, and
testing has not been exhaustive. I will add a description and some example of uses here shortly.
In the mean-time, here is an example of what the simulation can do. It shows a set of species evolving according to a fitness landscape (denoted by grey contour lines), with competition amongst them (proportional to their closeness in niche space). The simulation can have quite complicated fitness landscapes and also allows for flexibility in the competition kernel. New species bud off from existing species according to a random birth process.

![animated gif of simulation](tester.gif)




