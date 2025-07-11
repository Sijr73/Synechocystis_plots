# Loads model saved as .ods file ###############################################

setwd(paste(directory,"/Models",sep=""))

odsfile <- paste(modelname,".ods", sep = "")

#nsheets    <- get_num_sheets_in_ods(odsfile)
sheets <- list_ods_sheets(odsfile)
nsheets    <- length(sheets)
# Position of parameters
posM          <- (1:nsheets)[sheets == "M"]
posK          <- (1:nsheets)[sheets == "K"]
posKI         <- (1:nsheets)[sheets == "KI"]
posKA         <- (1:nsheets)[sheets == "KA"]
poskcat       <- (1:nsheets)[sheets == "kcat"]
posconditions <- (1:nsheets)[sheets == "conditions"]
#posx          <- (1:nsheets)[sheets == "x"]
#posq          <- (1:nsheets)[sheets == "q"]

# Reads q0, qT, T
#q0 <- as.numeric(read_ods(odsfile, sheet= posq)[1,-1])
#qT <- as.numeric(read_ods(odsfile, sheet= posq)[2,-1])
#T  <- as.numeric(read_ods(odsfile, sheet= posq)[3,2])


# Reads x(t) and dxdt(t) functions
#eval(parse(text=noquote(read_ods(odsfile, sheet= posx)[,2])))

# Getting data from sheets

# Mass fraction matrix Mtotal including external reactants
Mtotal <- as.matrix(read_ods(odsfile, sheet= posM)[,-1])

# reaction and reactant names
reaction <- colnames(Mtotal)
rownames(Mtotal) <- unlist(read_ods(odsfile, sheet= posM)[,1])

# Michaelis constant matrix K 
if (length(posK) > 0) {
  
  K <- as.matrix(read_ods(odsfile, sheet= posK)[,-1])
  
} else K <- 0.1*(Mtotal<0)

# inhibition constant matrix KI
if (length(posKI) > 0) {
  
  KI <- as.matrix(read_ods(odsfile, sheet= posKI)[,-1])
  
} else KI <- 0*K

# activation constant matrix KA
if (length(posKA) > 0) {
  
  KA <- as.matrix(read_ods(odsfile, sheet= posKA)[,-1])
  
} else KA <- 0*K


# kcat
kcatf <- as.numeric(read_ods(odsfile, sheet= poskcat)[1,-1])
kcatb <- as.numeric(read_ods(odsfile, sheet= poskcat)[2,-1])
names(kcatf)=reaction
kcatb <- as.numeric(read_ods(odsfile, sheet= poskcat)[2,-1])
names(kcatb)=reaction
kcat_matrix <- rbind(kcatf, kcatb)
rownames(kcat_matrix)[1]="kcat_f"
rownames(kcat_matrix)[2]="kcat_b"
# Growth condition names
condition <- colnames(read_ods(odsfile, sheet= posconditions))[-1]

# Mass density rho at each condition
rho_cond <- as.numeric(read_ods(odsfile, sheet= posconditions)[1,-1])

# external concentrations at each condition
x_cond  <- as.matrix(read_ods(odsfile, sheet= posconditions)[-1,-1])
##conditionsMatrix
conditions_matrix=as.matrix(read_ods(odsfile, sheet= posconditions))
reactant <- rownames(Mtotal)

# Definitions ##################################################################

# index for external reactants
n <- 1:dim(x_cond)[1]

# internal matrix M
M <- Mtotal[-n,]

# number of external reactants
nx <- dim(x_cond)[1]

# number of growth conditions
n_conditions <- dim(x_cond)[2]

# names of internal reactants
i_reactant <- reactant[-n]

# number of internal reactants
p <- dim(M)[1]

# number of reactions
r <- dim(M)[2]

# the sum of each M column 
sM <- colSums(M)

# delete numerical artifacts
sM[abs(sM) < 1e-10] <- 0

# indexes for reactions: s (transport), e (enzymatic), and ribosome r 

e <- c(1:(r-1))[sM[1:(r-1)] == 0]  

s <- c(1:(r-1))[sM[1:(r-1)] != 0] 

# indexes: m (metabolite), a (all proteins)

m <- 1:(p-1)

# number of transport reactions
ns <- length(s)

setwd(directory)

