library(tidyverse)
library(raster)
library(furrr) # To parallelize!

### Open data
boliv <- raster("/nfs/agfrontiers-data/Remote Sensing/KS files/classi_bol_dry_2008_102033.tif")

### Function to calculate % neighbors of different types
percent_different_neighbors <- function(cell) {
  # pull neighboring values
  # see ?adjacent for the arguments used here
  x <- getValues(boliv)[adjacent(boliv, 
                                 cells = cell, 
                                 directions = 8, 
                                 pairs = F, 
                                 include = T)]
  # calculate the percent of them that are not equal to the focal cell (i.e., not equal to the first value of x)
  perc_diff <- sum(x!=x[1])/(length(x)-1)*100
  return(perc_diff)
}

### Optimized version of the function above, that uses pre-calculated getValues (because that only needs to be done once)
percent_different_neighbors_optimized <- function(cell, r, r_values) {
  # pull neighboring values
  # see ?adjacent for the arguments used here
  x <- r_values[adjacent(r, 
                         cells = cell, 
                         directions = 8, 
                         pairs = F, 
                         include = T)]
  # calculate the percent of them that are not equal to the focal cell (i.e., not equal to the first value of x)
  perc_diff <- sum(x!=x[1])/(length(x)-1)*100
  return(perc_diff)
}

# Get the values of the Bolivia raster to be used later.
b_values <- getValues(boliv)
ncells <- ncell(boliv)

#### TESTING ####
n <- 100 # for testing purposes

# In each test, the call is wrapped in system.time() which tells you how long the expression took to ran.
# Also, different random points are sampled from the raster each time for testing purposes.
# Later you can submit the full job to the cluster.

# Time to run original function on 100 random points: about 2 minutes
system.time(
  map_dbl(sample(1:ncells, n),  ~percent_different_neighbors(.))
)
# Time to run new optimized function on 100 random points: about 0.4 seconds
system.time(
  map_dbl(sample(1:ncells, n), ~percent_different_neighbors_optimized(., boliv, b_values))
)

# Set up parallel with furrr.
options(mc.cores = 8) # There are eight cores on a slurm node.
plan(multicore) # Specify that the job will be run in parallel across multiple cores.

# Compare non-parallel to parallel. Increase the n for this test, because there is a little overhead involved in starting the parallel jobs.
n_big <- 100000
# Time to run optimized function, not in parallel, on 100K random points: about 5 minutes
system.time(
  map_dbl(sample(1:ncells, n_big), ~percent_different_neighbors_optimized(., boliv, b_values))
)
# Time to run optimized function, in parallel across 8 cores on a single Slurm node, on 100K random points: 69 seconds or just over a minute
system.time(
  future_map_dbl(sample(1:ncells, n_big), ~ percent_different_neighbors_optimized(., boliv, b_values))
)

# Based on this, if there are 60 million points, we can estimate the time the Slurm job will take to run
# The 100K points, processed in parallel, ran in 69 seconds.
# Multiply this by the total number of points divided by the number we just ran, divided by 60*60 to get the number of hours.

69 * ncells / n_big / (60 * 60) # About 12 hours!
