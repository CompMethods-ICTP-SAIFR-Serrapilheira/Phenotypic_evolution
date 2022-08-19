# Scientific programming course ICTP-Serrapilheira
# First version: August 18th 2022
# Run the manduca simulations

#Load packages
library(MASS)
library(dplyr)
library(tidyr)
library(plyr)

#Executing the functions
source("fct/manduca.R")
source("fct/cor2cov.R")

# Generating the simulations
manduca(1, 0.1)
manduca(1, 0.5)
manduca(1, 1)

#....................................
