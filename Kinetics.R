# Turnover time function for different kinetic rate laws
# NOTE: external concentrations x are treated as global variables, so they are 
# not an input of the functions

# first separate Km matrices for substrates and for products
KS <- K*(Mtotal<0)
KP <- K*(Mtotal>0)


# irreversible MM ##############################################################
iMM <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  tauj <- rep(0,r)
  
  for (j in 1:r) {
    
    tauj[j] <- as.numeric(prod(1 + KS[,j]/xc)/kcatf[j]) 
    
  }
  
  return(tauj)
}

# derivative with respect to c

diMM <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  ditauj <- matrix(rep(0,r*p),ncol=p)
  
  for (j in 1:r) {
    
    for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants xc
      
      ditauj[j,i2] <- -as.numeric((KS[y,j]/(c[i2]^2))*prod(1 + 
                                                             KS[-y,j]/xc[-y])/kcatf[j]) 
      
    }
    
  }
  
  return(ditauj)
  
}

# derivative with respect to x

diMMx <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  ditauj <- matrix(rep(0,r*nx),ncol=nx)
  
  for (j in 1:r) {
    
    for (n in 1:nx) {
      
      ditauj[j,n] <- -as.numeric((KS[n,j]/(x[n]^2))*prod(1 + 
                                                           KS[-n,j]/xc[-n])/kcatf[j]) 
      
    }
    
  }
  
  return(ditauj)
  
}

# irreversible MM + inhibition (here only one inhibitor per reaction) ##########

# (it's useful to first define the matrix of reciprocal KI, with zero for KI=0)

rKI <- 1/KI
rKI[rKI == Inf] <- 0

iMMi <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  tauj <- rep(0,r)
  
  for (j in 1:r) {
    
    # irreversible MM with inhibition
    tauj[j] <- as.numeric(prod(1 + xc*rKI[,j])*prod(1 + KS[,j]/xc)/kcatf[j]) 
    
  }
  
  return(tauj)
}

# derivative with respect to ci

diMMi <- function(c){
  
  # all concentrations x, c
  xc <- c(x, c) 
  
  ditauj <- matrix(rep(0,r*p),ncol=p)
  
  for (j in 1:r) {
    
    for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants
      
      ditauj[j,i2] <- as.numeric(rKI[y,j]*prod(1 + KS[,j]/xc) - 
                                   prod(1 + xc*rKI[,j])*(KS[y,j]/(c[i2]^2))*prod(1 + KS[-y,j]/xc[-y]))/kcatf[j] 
      
    }
    
  }
  
  return(ditauj)
  
}

# irreversible MM + activation (only one activator per reaction) ####################################################

iMMa <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  tauj <- rep(0,r)
  
  for (j in 1:r) {
    
    # irreversible MM with inhibition
    tauj[j] <- as.numeric(prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc)/kcatf[j]) 
    
  }
  
  return(tauj)
}

# derivative with respect to c

diMMa <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  ditauj <- matrix(rep(0,r*p),ncol=p)
  
  for (j in 1:r) {
    
    for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants
      
      ditauj[j,i2] <- -as.numeric( (KA[y,j]/(c[i2]^2))*prod(1 + KS[,j]/xc) + 
                                     (KS[y,j]/(c[i2]^2))*prod(1 + KA[,j]/xc)*prod(1 + KS[-y,j]/xc[-y]) )/kcatf[j] 
      
    }
    
  }
  
  return(ditauj)
  
}

# irreversible MM + inhibition + activation ####################################################

iMMia <- function(c){
  
  # all concentrations c
  xc <- c(x, c)
  
  tauj <- rep(0,r)
  
  for (j in 1:r) {
    
    # irreversible MM with inhibition
    tauj[j] <- as.numeric(prod(1 + xc*rKI[,j])*prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc)/kcat[j]) 
    
  }
  
  return(tauj)
}

# derivative with respect to c (assuming only one inhibitor and activator per reaction) 

diMMia <- function(c){
  
  # all concentrations c
  xc <- c(x, c)
  
  ditauj <- matrix(rep(0,r*p),ncol=r)
  
  for (j in 1:r) {
    
    for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants
      
      ditauj[j,i2] <- as.numeric( rKI[y,j]*prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc) + 
                                    
                                    prod(1 + xc*rKI[,j])*(-KA[y,j]/(c[i2]^2))*prod(1 + KS[,j]/xc) + 
                                    
                                    prod(1 + xc*rKI[,j])*prod(1 + KA[,j]/xc)*(-KS[y,j]/(c[i2]^2))*prod(1 + KS[-y,j]/xc[-y]) )/kcat[j] 
      
    }
    
  }
  
  return(ditauj)
  
}

# Reversible MM  ######################################################################################

# first separate Km matrices for substrates and for products
KS <- K*(Mtotal<0)
KP <- K*(Mtotal>0)

rMM <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  tauj <- rep(0,r)
  
  for (j in 1:r) {
    
    tauj[j] <- as.numeric( 1 / ( kcatf[j]/prod(1 + KS[,j]/xc) - kcatb[j]/prod(1 + KP[,j]/xc)  ) )
    
  }
  
  return(tauj)
}

# derivative with respect to c #################################################

drMM <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  ditauj <- matrix(rep(0,r*p),ncol=p)
  
  for (j in 1:r) {
    
    for (i in 1:p) {
      
      y <- i + nx # position of i in terms of all reactants xc
      
      # this will be multiplied my -(tau(c)^2) after loops
      ditauj[j,i] <- as.numeric( (kcatf[j]/prod(1 + KS[-y,j]/xc[-y]))*KS[y,j]/((c[i] + KS[y,j])^2) 
                                 
                                 
                                 -(kcatb[j]/prod(1 + KP[-y,j]/xc[-y]))*KP[y,j]/((c[i] + KP[y,j])^2) ) 
      
    }
    
    ditauj[j,] <-  -ditauj[j,]*as.numeric( 1 / ( kcatf[j]/prod(1 + KS[,j]/xc) - kcatb[j]/prod(1 + KP[,j]/xc)  )^2 )
    
  }
  
  return(ditauj)
  
}

# derivative with respect to x

drMMx <- function(c){
  
  # all concentrations x,c
  xc <- c(x, c) 
  
  ditauj <- matrix(rep(0,r*nx),ncol=nx)
  
  for (j in 1:r) {
    
    for (n in 1:nx) {
      
      ditauj[j,n] <- as.numeric( (kcatf[j]/prod(1 + KS[-n,j]/xc[-n]))*KS[n,j]/((x[n] + KS[n,j])^2) 
                                 
                                 
                                 -(kcatb[j]/prod(1 + KP[-n,j]/xc[-n]))*KP[n,j]/((x[n] + KP[n,j])^2) )
      
    }
    
    ditauj[j,] <-  -ditauj[j,]*as.numeric( 1 / ( kcatf[j]/prod(1 + KS[,j]/xc) - kcatb[j]/prod(1 + KP[,j]/xc)  )^2 )
    
  }
  
  return(ditauj)
  
}


#############################################################################################
# now selects the simplest rate law necessary for optimization efficiency

# irreversible cases
if (sum(kcatb)== 0) {
  
  # 1) irreversible MM kinetics
  if (sum(KI) == 0) {
    
    tau <- iMM 
    
    dtau <- diMM
    
    epsilon <- diMMx
    
  }
  
  # 2) irreversible MM kinetics + inhibition 
  if (sum(KI) > 0)  {
    
    tau <- iMMi 
    
    dtau <- diMMi
    
    epsilon <- diMMx
    
  }
  
  # 3) irreversible MM kinetics + activation 
  if (sum(KA) > 0)  {
    
    tau <- iMMa
    
    dtau <- diMMa
    
    epsilon <- diMMx
    
  }
  
}

# reversible cases
if (sum(kcatb)> 0) {
  
  tau <- rMM
  
  dtau <- drMM
  
  epsilon <- drMMx
  
}



