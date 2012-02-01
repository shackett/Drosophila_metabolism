###### Read in absorbance posterior medians for enzyme activity traits as well as flight and respirometry traits

setwd("/Users/seanhackett/Desktop/Cornell/Drosophila_metabolism/")

#if use.line == TRUE, then the point estimates of lines will be used. otherwise population estimates will be used
use.line = TRUE
#if use.flight == TRUE, then the goal is to match the rate of ATP consumption needed to match the maximum flight velocity
#if use.flight == FALSE, then the goal is to match rates of O2 consumption and CO2 production
use.flight = TRUE

if(use.line == TRUE){

line.enzyme <- read.table("MCMC_files/lnMat.tsv", sep = "\t", header = TRUE)
line.flight <- read.table("MCMC_files/FLTlnMns.tsv", sep = "\t", header = TRUE)
rownames(line.enzyme) <- as.character(line.enzyme$ln_nam)
rownames(line.flight) <- as.character(line.flight$ln_nam)
line.enzyme <- line.enzyme[,-1]; line.flight <- line.flight[,-1]; line.resp <- line.resp[,-1]
#setting negative vmax values to 0
line.enzyme[line.enzyme < 0] <- 0

if(use.flight == FALSE){
line.resp <- read.table("MCMC_files/RESPlnMns.tsv", sep = "\t", header = TRUE)
rownames(line.resp) <- as.character(line.resp$ln_nam)
	}
if(use.flight == TRUE){
flightM <- read.table("MCMC_files/FLTlnMns.tsv", sep = "\t", header = TRUE)
line.flight <- data.frame(maxVel = flightM[,colnames(flightM) %in% "maxv"])
rownames(line.flight) <- flightM$"ln_nam"
	}
}

if(use.line == FALSE){

pop.enzyme <- read.table("MCMC_files/popMat.tsv", sep = "\t", header = TRUE)
pop.flight <- read.table("MCMC_files/FLTpopMns.tsv", sep = "\t", header = TRUE)
rownames(pop.enzyme) <- as.character(pop.enzyme$pop_nam)
rownames(pop.flight) <- as.character(pop.flight$pop_nam)
pop.enzyme <- pop.enzyme[,-1]; pop.flight <- pop.flight[,-1]; pop.resp <- pop.resp[,-1]

if(use.flight == FALSE){
pop.resp <- read.table("MCMC_files/RESPpopMns.tsv", sep = "\t", header = TRUE)
rownames(pop.resp) <- as.character(pop.resp$pop_nam)
	}
if(use.flight == TRUE){
flightM <- read.table("MCMC_files/FLTpopMns.tsv", sep = "\t", header = TRUE)
pop.flight <- data.frame(maxVel = flightM[,colnames(flightM) %in% "maxv"])
rownames(pop.flight) <- flightM$"pop_nam"
	}
}

load("drosophila_stoi.R")
source("path_length_calc.R")
source("flimsy_plate_extract.R")
library(limSolve)
library(gplots)
library("colorRamps")

if(use.line == TRUE){
	if(use.flight == FALSE){
	valid.samples <- intersect(rownames(line.enzyme), rownames(line.resp))
		}else{
		valid.samples <- intersect(rownames(line.enzyme), rownames(line.flight))
			}
		
	line.enzyme <- line.enzyme[rownames(line.enzyme) %in% valid.samples,]
	if(use.flight == FALSE){
	line.resp <- line.resp[rownames(line.resp) %in% valid.samples,]
		}
		
	#relate lines to their population
	populations <- c("N", "I", "B", "T", "Z")
	line.pop <- rep(NA, times = length(valid.samples))
	pop.color <- data.frame(pops = populations, color = c("RED", "ORANGE", "BLUE", "GREEN", "CYAN"), long.name = c("Netherlands", "Ithaca", "Beijing", "Tasmania", "Zimbabwee"))
	line.color <- rep(NA, times = length(valid.samples))
	
	for(pop in 1:length(populations)){
		line.pop[grep(populations[pop], valid.samples)] <- populations[pop]
		line.color[grep(populations[pop], valid.samples)] <- as.character(pop.color[pop,2])
		}}









##### a few functions

find.name <- function(kinetic_enzyme, assay_parameters){
	assay_parameters$Model_name[rownames(assay_parameters) == kinetic_enzyme]
	}

velocity.to.power <- function(velocity.vector, power.lm){
degree = length(power.lm$coef)-1
apply(matrix(power.lm$coef, ncol = degree +1, nrow = length(velocity.vector), byrow = TRUE)*(matrix(velocity.vector, ncol = degree +1, nrow = length(velocity.vector), byrow = FALSE)^matrix((0:(length(power.lm$coef)-1)), ncol = degree +1, nrow = length(velocity.vector), byrow = TRUE)), 1, sum)}


###### Read in and calculate the paramters needed to transform A/s into moles/s using the beer-lambert law and path-length calculated from the well-geometry and volume of a microtiter-assay

assay_parameters <- read.table("kin_par.csv", sep = ",", header = TRUE, colClasses = c("character", "numeric", "character", "numeric", "numeric", "numeric", "character"))
rownames(assay_parameters) <- assay_parameters$Enzyme
assay_parameters <- assay_parameters[,-1]
assay_parameters <- cbind(assay_parameters, data.frame(path_length = NA))

for(i in 1:length(assay_parameters[,1])){

#If assay is visible then use a different type of plate is used
if(assay_parameters$molar_absorptivity[i] %in% c(16000, 19100, 11500)){
	assay_parameters$path_length[i] <- flimsy_pl(assay_parameters$assay_volume[i]*1e6)
	}else{
assay_parameters$path_length[i] <- find.pl(assay_parameters$assay_volume[i], height_scale, bottom_d)$minimum * 100	}}
	
###### Determine sample x reaction Vmax	
	
if(use.line == TRUE){
enzymeOD <- line.enzyme[,((colnames(line.enzyme) %in% rownames(assay_parameters)))]
	}else{	
enzymeOD <- pop.enzyme[,((colnames(pop.enzyme) %in% rownames(assay_parameters)))]}
enzyme_moles_per_second <- enzymeOD
enzyme_moles_per_second[!is.na(enzyme_moles_per_second)] <- NA
	
for (i in 1:length(enzyme_moles_per_second[1,])){
	enzyme <- colnames(enzyme_moles_per_second)[i]
	params <- assay_parameters[rownames(assay_parameters) == enzyme,]
	
	if(params$molar_absorptivity != "STD"){
	params$molar_absorptivity <- as.numeric(params$molar_absorptivity)
	
	enzyme_moles_per_second[,i] <- (((enzymeOD[,i]/(params$molar_absorptivity*params$path_length))*params$assay_volume*(params$assay_volume/20e-6))/params$OD_scaling)/params$fly_fraction
	}}

kinetic_enzymes = enzyme_moles_per_second[,apply(is.na(enzyme_moles_per_second), 2, sum) == 0]
kinetic_reactions <- unlist(lapply(colnames(kinetic_enzymes), find.name, assay_parameters)) 

###### Determine the rate of ATP consumption to sustain maximum velocity

#from Sun et al. 2003:
#velocity in m/s to power in W/kg-flight-muscle for Virilis
vel.pow.corr <- data.frame(advance.ratio = c(0, 0.13, 0.27, 0.40, 0.53), velocity.virilis = c(0, 0.5, 1, 1.5, 2), power.virilis = c(29, 27, 25, 30, 40), velocity.mel = rep(NA, times = 5), power.mel = rep(NA, times = 5))

#from Lehman and Dickenson 1997:
#10% efficiency of converting ATP chemical energy into work
#60W/kg for mel hovering flight

#Tran and Unden 1998
#46.5kJ/mol

#for melanogaster, 0.85 velocity, advance ratio 0.32

#scale virilis velocity = f(advance ratio) to melanogaster
vel.pow.corr$velocity.mel <- vel.pow.corr$velocity.virilis*(0.85/((2/0.53)*0.32))

#scale virilis mechanical power to melanogaster
vel.pow.corr$power.mel <- vel.pow.corr$power.virilis*(60/vel.pow.corr$power.virilis[1])


power.lm <- lm(vel.pow.corr$power.mel ~ vel.pow.corr$velocity.mel + I(vel.pow.corr$velocity.mel^2) + I(vel.pow.corr$velocity.mel^3))
degree <- length(power.lm$coef)-1

velocities <- seq(0, 1.5, by = 0.05)
fit.velocity <- velocity.to.power(velocity.vector, power.lm)

#get total power from efficiency of conversion and weight
#gives Watts: J/s per fly
if(use.flight == TRUE){
plot(fit.velocity ~ velocities, type = "l")
if(use.line == TRUE){
power <- (velocity.to.power(as.numeric(unlist(line.flight)), power.lm)/0.1)*(line.enzyme$wts*10^-6/5)
points(velocity.to.power(as.numeric(unlist(line.flight)), power.lm) ~ unlist(line.flight), pch = 8, col = "RED")
	}else{
	power <- (velocity.to.power(as.numeric(unlist(pop.flight)), power.lm)/0.1)*(pop.enzyme$wts*10^-6/5)
	points(velocity.to.power(as.numeric(unlist(pop.flight)), power.lm) ~ unlist(pop.flight), pch = 8, col = "RED")
		}

#46500 J per mole ATP (and making the flux for 5 flies again)
ATP.consumed <- (power/46500)*5

}






#arbitrary constant to makes fluxes nicer to look at
SF <- 1e11

#correct by 10ul -> L, L -> moles, hrs -> seconds, single fly -> 5 flies
#giving moles/5fly-second

if(use.flight == FALSE){
if(use.line == TRUE){
gas_exchange = (line.resp[,c(1:2)]/1e7/22.4/3600)*5
}else{gas_exchange = (pop.resp[,c(1:2)]/1e7/22.4/3600)*5}
gas_exchange <- gas_exchange * SF
}else{
	ATP.consumed <- ATP.consumed * SF
	}

kinetic_enzymes <- kinetic_enzymes * SF

joint.stoi.red <- joint.stoi
	
S = joint.stoi.red
f = rep(0, times = length(joint.stoi.red[,1]))

#respirometry - irreversible reactions
if(use.flight == FALSE){
irreversible = c("Trehalose6P synthetase_c", "Trehalose6P phosphatase_c", "Trehalase_c", "Hexokinase_c", "Branching enzyme_c", "Glycogen phosphorylase_c", "Glycogen Synthase_c", "UDP-pyrophosphorylase_c", "Lactonase_c", "6-phosphogluconate dehydrogenase_c", "Phosphofructokinase_c", "Glycerol 3-phosphate shuttle", "Glycerol 3-phosphate dehydrogenase (NAD)_c", "Phosphoenolpyruvate carboxylase_c", "Pyrophosphatase_c", "Isocitrate dehydrogenase (NADP)_m", "Isocitrate dehydrogenase (NADP)_m", "Alpha-ketoglutarate dehydrogenase I_m", "ATP synthetase", "NADH oxidation for proton transport", "FADH2 oxidation for proton transport", "Malic enzyme (NADP)_c", "Malic enzyme (NAD)_c", "C8 synthesis_c", "Pyruvate dehydrogenase_m", "ATPase_c", "fructose 1,6 bisphosphatase_c")
}else{
#flight-flux constraints
irreversible = c("Trehalose6P phosphatase_c", "Hexokinase_c", "Branching enzyme_c", "Lactonase_c", "6-phosphogluconate dehydrogenase_c", "Phosphofructokinase_c", "Glycerol 3-phosphate shuttle", "Glycerol 3-phosphate dehydrogenase (NAD)_c", "Phosphoenolpyruvate carboxylase_c", "Pyrophosphatase_c", "Isocitrate dehydrogenase (NADP)_m", "Isocitrate dehydrogenase (NADP)_m", "Alpha-ketoglutarate dehydrogenase I_m", "ATP synthetase", "NADH oxidation for proton transport", "FADH2 oxidation for proton transport", "Malic enzyme (NADP)_c", "Malic enzyme (NAD)_c", "C8 synthesis_c", "Pyruvate dehydrogenase_m", "ATPase_c", "fructose 1,6 bisphosphatase_c")
}
#sample matricies: kinetic_enzymes & gas_exchange

#iterate through all samples
nsamples <- length(kinetic_enzymes[,1])

calc.fluxes <- matrix(NA, nrow = length(S[1,]), ncol = nsamples)
rownames(calc.fluxes) <- colnames(S)
colnames(calc.fluxes) <- rownames(kinetic_enzymes)
optim.resid <- rep(NA, times = nsamples)

#line <- c(1:nsamples)[rownames(kinetic_enzymes) %in% "N15"]

for (line in c(1:nsamples)){

#reactions corresponding to each measured Vmax

G_flux <- matrix(data = NA, ncol = length(kinetic_enzymes[1,]), nrow = length(joint.stoi.red[1,])) 
H_flux <- rep(NA, times = length(kinetic_enzymes[1,]))
for (i in 1:length(kinetic_enzymes[1,])){
	G_flux[,i] <- ifelse(colnames(joint.stoi.red) == kinetic_reactions[i], -1, 0)
	H_flux[i] <- -1*kinetic_enzymes[line,i]
	#exception for bidirectional vmax measurements
	if(kinetic_reactions[i] %in% c("Malate dehydrogenase_m", "Phosphoglucomutase_c")){
		G_flux <- cbind(G_flux, ifelse(colnames(joint.stoi.red) == kinetic_reactions[i], 1, 0))
		H_flux <- c(H_flux, -1*kinetic_enzymes[line,i])
		}
	}
G_flux <- t(G_flux)

G_irr = matrix(data = NA, ncol = length(irreversible), nrow = length(joint.stoi.red[1,]))
for (i in 1:length(irreversible)){
G_irr[,i] = ifelse(colnames(joint.stoi.red) == irreversible[i], 1, 0)
	}
G_irr <- t(G_irr)
H_irr = rep(0, times = length(irreversible))

G <- rbind(G_flux, G_irr)
h <- c(H_flux, H_irr)

#for flight, remove palmitate as an energy source
if(use.flight == TRUE){
	G <- rbind(G, ifelse(colnames(joint.stoi.red) == "Palmitate usage", -1, 0), ifelse(colnames(joint.stoi.red) == "Palmitate usage", 1, 0))
	h <- c(h, 0, 0)
	}
#colnames(G) <- colnames(S)

if(use.flight == FALSE){
#reactions for CO2 and O2 exchange
exchange_rxns = c("CO2 leaving", "O2 entering")
A = matrix(data = NA, ncol = length(exchange_rxns), nrow = length(joint.stoi.red[1,]))
for (i in 1:length(exchange_rxns)){
	A[,i] <- ifelse(colnames(joint.stoi.red) == exchange_rxns[i], 1, 0)
	}
u = gas_exchange[line,c(1,2)]
	}else{
		#maximize energy generation and minimize carbon usage (sum of trehalase and GP)
		exchange_rxns = c("ATPase_c", "Trehalase_c", "Glycogen phosphorylase_c")
		A = matrix(data = NA, ncol = length(exchange_rxns), nrow = length(joint.stoi.red[1,]))
		for (i in 1:length(exchange_rxns)){
			if(exchange_rxns[i] == "Trehalase_c"){
			A[,i] <- ifelse(colnames(joint.stoi.red) == exchange_rxns[i], 2, 0)				}else{
				A[,i] <- ifelse(colnames(joint.stoi.red) == exchange_rxns[i], 1, 0)
				}
			}
		u = t(c(ATP.consumed[line], 0, 0))
		}


#also maximize ATP production
#A = cbind(A, ifelse(colnames(joint.stoi.red) == "ATPase_c", 1, 0))
#u = c(unlist(u), 10)

#QP.optim <- lsei(A = t(A), B = t(u), E = S, F = f, G = G, H = h)
#calc.fluxes[,line] <- QP.optim$X
#optim.resid[line] <- QP.optim$solutionNorm

calc.fluxes[,line] <- linp(E = S, F = f, G = G, H = h, Cost = -1*(A[,1] - A[,2] - A[,3]), ispos = FALSE)$X

}

nonzero.flux <- calc.fluxes[apply(calc.fluxes == 0, 1, sum) != length(calc.fluxes[1,]),]

if(use.line == TRUE){pdf(file = "line_flux.pdf")}else{pdf(file = "pop_flux.pdf")}

if(use.line == TRUE){
heatmap.2((nonzero.flux - apply(nonzero.flux, 1, mean))/apply(nonzero.flux, 1, sd), trace = "none", col = blue2yellow(100), ColSideColors = line.color, cexCol = 0.4, cexRow = 0.05 + 0.8/log10(length(nonzero.flux[,1])))
	}else{
heatmap.2((nonzero.flux - apply(nonzero.flux, 1, mean))/apply(nonzero.flux, 1, sd), trace = "none", col = blue2yellow(100), cexRow = 0.05 + 0.8/log10(length(nonzero.flux[,1])))
	}
#calc.fluxes[rownames(calc.fluxes) %in% c("CO2 leaving", "O2 entering"),optim.resid > 0.0001]

#(line.resp[optim.resid > 0.0001,]/1e7/22.4/3600)*SF*5
#line.enzyme[optim.resid > 0.0001,]
                               
zero.flux <- colnames(S)[apply(calc.fluxes == 0, 1, sum) == length(calc.fluxes[1,])]

if(use.line == TRUE){

plot(log((svd(nonzero.flux, nu = 1, nv = 1)$d^2/sum(svd(nonzero.flux, nu = 2, nv = 2)$d^2))), main = "Variance explained by each PC (ln)", ylab = "log(proportion of variance explained by PC)")
n.components <- 4
	
for (pc in 1:n.components){
	
plot(svd(nonzero.flux, nu = n.components, nv = n.components)$v[,pc], col = line.color, main = paste("Principle Component ", pc, sep = ""), xlab = paste("PC ", pc, sep = ""), ylab = "PC value", pch = 8, cex = 1.2)
legend("topleft", legend = pop.color[,3], text.col = as.character(pop.color[,2]))
	
	
	}}; dev.off()

#### for measured enzymes, determine v / vmax

od.measured.carried.flux <- matrix(NA, ncol = length(kinetic_reactions), nrow = length(kinetic_enzymes[,1]))

for (rxn in 1:length(kinetic_reactions)){
	od.measured.carried.flux[,rxn] <- calc.fluxes[rownames(calc.fluxes) %in% kinetic_reactions[rxn],]
	}

kinetic_enzymes[kinetic_enzymes == 0] <- NA

Vmax.fraction <- od.measured.carried.flux/kinetic_enzymes

if(use.line == TRUE){
	write.table(calc.fluxes, file = "line_fluxes.tsv", sep = "\t")
	write.table(Vmax.fraction, file = "line_vmax_fraction.tsv", sep = "\t")
	} else {
		write.table(calc.fluxes, file = "pop_fluxes.tsv", sep = "\t")
		write.table(Vmax.fraction, file = "pop_vmax_fraction.tsv", sep = "\t")
		}
