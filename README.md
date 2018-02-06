# nichefillr
R package implementing niche filling simulations.

This is a very early alpha release. Currently there is little in the way of error-checking, and
testing has not been exhaustive. I will add a description and some example of uses here shortly.
In the mean-time, here is an example of what the simulation can do. It shows a set of species evolving according to a fitness landscape (denoted by grey contour lines), with competition amongst them (proportional to their closeness in niche space). The simulation can have quite complicated fitness landscapes and also allows for flexibility in the competition kernel. New species bud off from existing species according to a random birth process.

![animated gif of simulation](tester.gif)