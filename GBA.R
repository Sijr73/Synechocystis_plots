# GBA model solver

# Clear variables
rm(list=ls(all=TRUE))

require('rstudioapi') 
require('readODS')
require('nloptr')
require('Matrix')
require('MASS')
require('lpSolve')

# Sets working directory to source file location ###############################

directory <- dirname(getActiveDocumentContext()$path)

setwd(directory) 

# Model name here ##############################################################

modelname <- "ExtendedModel"

# Reads model saved as .ods file ###############################################

source("Readmodelods.R")

# kinetics #####################################################################

source("Kinetics_rxn.R")

# Optimization on f ############################################################

source("GBA_solver.R") 



# Exporting results ############################################################

# Exporting csv file with results #######

source("Exportcsv.R")

# Plots #################################

source("GBA_Plots.R")

#phi_opt
p_opt <- opt_state[,paste("p",reaction)]
phi_opt <- p_opt/c_opt[,p]
phi_opt=data.frame(phi_opt)
# export results
opt_state=data.frame(opt_state)
opt_state$phi.Carbon=phi_opt$p.Carbon.met.
opt_state$phi.PSII=phi_opt$p.PSII
opt_state$phi.PSI=phi_opt$p.PSI
opt_state$phi.ATPase=phi_opt$p.ATPase
opt_state$phi.Cyt6bf=phi_opt$p.Cyt6bf
opt_state$phi.Rubisco=phi_opt$p.Rubisco
opt_state$phi.AA_s=phi_opt$p.AA_s
opt_state$phi.Enzyme1=phi_opt$p.Enzyme1
opt_state$phi.PSU=phi_opt$p.Photosyn.
opt_state$phi.Ribosome=phi_opt$p.Ribosome
opt_state = opt_state[opt_state$convergence == 4, ]

Extendedmodel04032024=opt_state
save(Extendedmodel04032024,file='~/Downloads/GBA R/Extendedmodel04032024.Rdata')

