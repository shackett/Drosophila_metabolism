use.cluster = FALSE
#if use.line == TRUE, then the point estimates of lines will be used. otherwise population estimates will be used
use.line = TRUE
#if use.flight == TRUE, then the goal is to match the rate of ATP consumption needed to match the maximum flight velocity
#if use.flight == FALSE, then the goal is to match rates of O2 consumption and CO2 production
use.flight = FALSE
#if use.mcmc == TRUE, then FBA will be done on all sample within a matrix
use.mcmc = TRUE


setwd("/Users/Sean/Desktop/Cornell/Drosophila_metabolism")
load("drosophila_stoi.R")
source("path_length_calc.R")
#source("flimsy_plate_extract.R")
library(limSolve)
library(gplots)
library(colorRamps)
library(stringr)		
library(reshape2)

n_mcmc_samples <- 10000

rxNames <- c(colnames(joint.stoi)[!(colnames(joint.stoi) %in% c("Trehalose usage", "Glycogen usage", "Trehalase_c", "Phosphoglucomutase_c"))], "Glucose uptake")
joint.stoi <- joint.stoi[,!(colnames(joint.stoi) %in% c("Trehalose usage", "Glycogen usage", "Trehalase_c", "Phosphoglucomutase_c"))]
joint.stoi <- data.frame(joint.stoi, GlucoseUptake = ifelse(rownames(joint.stoi) == "D-Glucose_c", 1, 0))
colnames(joint.stoi) <- rxNames

##### a few functions

find.name <- function(kinetic_enzyme, assay_parameters){
	assay_parameters$Model_name[rownames(assay_parameters) == kinetic_enzyme]
	}

velocity.to.power <- function(velocity.vector, power.lm){
degree = length(power.lm$coef)-1
apply(matrix(power.lm$coef, ncol = degree +1, nrow = length(velocity.vector), byrow = TRUE)*(matrix(velocity.vector, ncol = degree +1, nrow = length(velocity.vector), byrow = FALSE)^matrix((0:(length(power.lm$coef)-1)), ncol = degree +1, nrow = length(velocity.vector), byrow = TRUE)), 1, sum)}


### using only a single sample or initializing MCMC samples
if(use.line == TRUE){
  line.enzyme <- read.table("MCMC_files/lnMat.tsv", sep = "\t", header = TRUE)
  line.flight <- read.table("MCMC_files/FLTlnMns.tsv", sep = "\t", header = TRUE)
  rownames(line.enzyme) <- as.character(line.enzyme$ln_nam)
  rownames(line.flight) <- as.character(line.flight$ln_nam)
  line.enzyme <- line.enzyme[,-1]
  #setting negative vmax values to 0
  line.enzyme[line.enzyme < 0] <- 0
  
  if(use.flight == FALSE){
    line.resp <- read.table("MCMC_files/RESPlnMns.tsv", sep = "\t", header = TRUE)
    rownames(line.resp) <- as.character(line.resp$ln_nam)
    line.resp <- line.resp[,-1]  
  }
  if(use.flight == TRUE){
    flightM <- read.table("MCMC_files/FLTlnMns.tsv", sep = "\t", header = TRUE)
    line.flight <- data.frame(maxVel = flightM[,colnames(flightM) %in% "maxv"])
    rownames(line.flight) <- flightM$"ln_nam"
  }
}


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
	crosspops <- c("NN", "BB")
	line.pop <- rep(NA, times = length(valid.samples))
	pop.color <- data.frame(pops = c(populations, crosspops), color = c("RED", "ORANGE", "BLUE", "GREEN", "CYAN", "PURPLE", "YELLOW"), long.name = c("Netherlands", "Ithaca", "Beijing", "Tasmania", "Zimbabwee", "Netherlands-Intrapop", "Beijing-Intrapop"))
	line.color <- rep(NA, times = length(valid.samples))
	
	for(pop in 1:length(populations)){
		line.pop[grep(populations[pop], valid.samples)] <- populations[pop]
		line.color[grep(populations[pop], valid.samples)] <- as.character(pop.color[pop,2])
		}
	for(pop in 1:length(crosspops)){
		pop_matches <- sapply(valid.samples, function(line){
			
			if(length(unique(unlist(strsplit(crosspops[pop], "")))) == 1){
				#intrapop
				length(unlist(str_match_all(line, unlist(strsplit(crosspops[pop], ""))[1]))) == 2
				}else{
					#interpop
					length(unlist(str_match_all(line, unlist(strsplit(crosspops[pop], ""))[1]))) != 0 & length(unlist(str_match_all(line, unlist(strsplit(crosspops[pop], ""))[2]))) != 0
					}
			})
				
		line.pop[pop_matches] <- crosspops[pop]
		line.color[pop_matches] <- as.character(pop.color[length(populations) + pop,2])
		}	
		
	}
		

if(use.mcmc == TRUE){
	
	#using mcmc samples
	
  # read in gas-exchange trait markov samples
	mcmc_files <- list.files("MCMC_files/MCMC_inbred")[grep("RESP", list.files("MCMC_files/MCMC_inbred"))]
	mcmc_matrix <- NULL
	for(mcmc_file in mcmc_files){
		source(paste("MCMC_files/MCMC_inbred/", mcmc_file, sep = ""))
		mcmc_matrix <- rbind(mcmc_matrix, get(ls()[grep("RESP", ls())]))
		rm(list = ls()[grep("RESP", ls())])
		}
	
  # read in flight trait markov samples
  mcmc_files <- list.files("MCMC_files/MCMC_inbred")[grep("FLT", list.files("MCMC_files/MCMC_inbred"))]
  mcmc_flight <- NULL
	for(mcmc_file in mcmc_files){
		source(paste("MCMC_files/MCMC_inbred/", mcmc_file, sep = ""))
		mcmc_flight <- rbind(mcmc_flight, get(ls()[grep("FLT", ls())]))
		rm(list = ls()[grep("FLT", ls())])
		}
  
  
	mcmc_list = list()
	
  # convert respirometry traits to a list
  
	line.resp <- read.table("MCMC_files/RESPlnMns.tsv", sep = "\t", header = TRUE)
	rownames(line.resp) <- as.character(line.resp$ln_nam)
	line.resp <- line.resp[,-1]
	
	mcmc_list$Vcot = mcmc_matrix[,grep("mu.ln\\[Vcot,", colnames(mcmc_matrix))]
	mcmc_list$Vox = mcmc_matrix[,grep("mu.ln\\[Vox,", colnames(mcmc_matrix))]
	colnames(mcmc_list$Vcot) <- rownames(line.resp); colnames(mcmc_list$Vox) <- rownames(line.resp)
	
	
  
  line.flt <- read.table("MCMC_files/FLTlnMns.tsv", sep = "\t", header = TRUE)
  rownames(line.flt) <- as.character(line.flt$ln_nam)
	line.flt <- line.flt[,-1]
	
  
  flight_formatting <- as.data.frame(t(sapply(colnames(mcmc_flight)[grep('mu.ln', colnames(mcmc_flight))], function(x){
    strsplit(regmatches(x, regexpr('\\[[a-z,0-9A-Z]+\\]', x)), split = '\\[|,|\\]')[[1]][-1]
  })))
	colnames(flight_formatting) <- c("trait", "line")
  flight_formatting$line <- rownames(line.flt)[as.numeric(flight_formatting$line)]
  flight_formatting$index <- grep('mu.ln', colnames(mcmc_flight))
  
  mcmc_flight <- mcmc_flight[,flight_formatting$index]
	mcmc_flight_melt <- melt(data.frame(flight_formatting, t(mcmc_flight)), id.vars = colnames(flight_formatting))
  mcmc_flight_cast <- acast(mcmc_flight_melt, formula = line ~ variable ~ trait, value.var = "value")
  #names(mcmc_flight_cast[1,1,])
  }


###### Read in and calculate the paramters needed to transform A/s into moles/s using the beer-lambert law and path-length calculated from the well-geometry and volume of a microtiter-assay

assay_parameters <- read.table("kin_par.csv", sep = ",", header = TRUE, colClasses = c("character", "numeric", "character", "numeric", "numeric", "numeric", "character"))
rownames(assay_parameters) <- assay_parameters$Enzyme
assay_parameters <- assay_parameters[,-1]
assay_parameters <- cbind(assay_parameters, data.frame(path_length = NA))

for(i in 1:length(assay_parameters[,1])){
	#If assay is visible then a different type of plate is used - determine the path-length in cm
	if(assay_parameters$molar_absorptivity[i] %in% c(16000, 19100, 11500)){
		assay_parameters$path_length[i] <- flimsy_pl(assay_parameters$assay_volume[i]*1e3)
		}else{
	assay_parameters$path_length[i] <- pl_optim(assay_parameters$assay_volume[i])$minimum * 100
			}
	}

####### Rxns corresponding to rows of G and stoichiometry are defined ############

joint.stoi.red <- joint.stoi
S = joint.stoi.red
f = rep(0, times = length(joint.stoi.red[,1]))

#respirometry - irreversible reactions
if(use.flight == FALSE){
irreversible = c("Trehalose6P synthetase_c", "Trehalose6P phosphatase_c", "Trehalase_c", "Hexokinase_c", "Branching enzyme_c", "Glycogen phosphorylase_c", "Glycogen Synthase_c", "UDP-pyrophosphorylase_c", "Lactonase_c", "6-phosphogluconate dehydrogenase_c", "Phosphofructokinase_c", "Glycerol 3-phosphate shuttle", "Glycerol 3-phosphate dehydrogenase (NAD)_c", "Phosphoenolpyruvate carboxylase_c", "Pyrophosphatase_c", "Isocitrate dehydrogenase (NADP)_m", "Isocitrate dehydrogenase (NADP)_m", "Alpha-ketoglutarate dehydrogenase I_m", "ATP synthetase", "NADH oxidation for proton transport", "FADH2 oxidation for proton transport", "Malic enzyme (NADP)_c", "Malic enzyme (NAD)_c", "C8 synthesis_c", "Pyruvate dehydrogenase_m", "ATPase_c", "Peroxide generation")
}else{
#flight-flux constraints
irreversible = c("Trehalose6P phosphatase_c", "Hexokinase_c", "Branching enzyme_c", "Lactonase_c", "6-phosphogluconate dehydrogenase_c", "Phosphofructokinase_c", "Glycerol 3-phosphate shuttle", "Glycerol 3-phosphate dehydrogenase (NAD)_c", "Phosphoenolpyruvate carboxylase_c", "Pyrophosphatase_c", "Isocitrate dehydrogenase (NADP)_m", "Isocitrate dehydrogenase (NADP)_m", "Alpha-ketoglutarate dehydrogenase I_m", "ATP synthetase", "NADH oxidation for proton transport", "FADH2 oxidation for proton transport", "Malic enzyme (NADP)_c", "Malic enzyme (NAD)_c", "Pyruvate dehydrogenase_m", "ATPase_c", "fructose 1,6 bisphosphatase_c")
}



	
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
	
	enzyme_moles_per_second[,i] <- ((((enzymeOD[,i]*params$assay_volume)/(params$molar_absorptivity*params$path_length))/params$OD_scaling)/params$fly_fraction)*params$Moles_of_absorbant
	
	}}

kinetic_enzymes = enzyme_moles_per_second[,apply(is.na(enzyme_moles_per_second), 2, sum) == 0]
kinetic_reactions <- unlist(lapply(colnames(kinetic_enzymes), find.name, assay_parameters)) 

kinetic_reactions <- kinetic_reactions[!(colnames(kinetic_enzymes) %in% c("PGM", "TRE"))]
kinetic_enzymes <- kinetic_enzymes[,!(colnames(kinetic_enzymes) %in% c("PGM", "TRE"))]


if(use.flight == TRUE){

###### Determine the rate of ATP consumption to sustain maximum velocity

#from Sun et al. 2003:
#velocity in m/s to power in W/kg-flight-muscle for Virilis
vel.pow.corr <- data.frame(advance.ratio = c(0, 0.13, 0.27, 0.40, 0.53), velocity.virilis = c(0, 0.5, 1, 1.5, 2), power.virilis = c(29, 27, 25, 30, 40), velocity.mel = rep(NA, times = 5), power.mel = rep(NA, times = 5))

#from Lehman and Dickenson 1997:
#10% efficiency of converting ATP chemical energy into work
#60W/kg for mel hovering flight

#Tran and Unden 1998
#46.5kJ/mole
#for melanogaster, 0.85 velocity, advance ratio 0.32, gives us the scaling between advance ratio and velocity

#scale virilis velocity = f(advance ratio) to melanogaster
vel.pow.corr$velocity.mel <- vel.pow.corr$velocity.virilis*(0.85/((2/0.53)*0.32))

#scale virilis mechanical power to melanogaster
vel.pow.corr$power.mel <- vel.pow.corr$power.virilis*(60/vel.pow.corr$power.virilis[1])


power.lm <- lm(vel.pow.corr$power.mel ~ vel.pow.corr$velocity.mel + I(vel.pow.corr$velocity.mel^2) + I(vel.pow.corr$velocity.mel^3))
degree <- length(power.lm$coef)-1

velocities <- seq(0, 1.5, by = 0.05)
fit.velocity <- velocity.to.power(velocities, power.lm)

#get total power from efficiency of conversion and weight
#gives Watts: J/s per fly
if(use.flight == TRUE){
pdf(file = "flight_vel.pdf")
plot(fit.velocity ~ velocities, type = "l", ylab = "Power (W)")
power.mech <- (velocity.to.power(as.numeric(unlist(line.flight$maxv)), power.lm)/0.1)*(line.enzyme$wts*10^-6/5)
points(velocity.to.power(as.numeric(unlist(line.flight$maxv)), power.lm) ~ unlist(line.flight$maxv), pch = 8, col = "RED")
dev.off()

#46500 J per mole ATP (and making the flux for 5 flies again)
ATP.consumed <- (power.mech/46500)*5
}
}





#arbitrary constant to makes fluxes nicer to look at
SF <- 1e11
kinetic_enzymes <- kinetic_enzymes * SF

#correct by 10ul -> L, L -> moles, hrs -> seconds, single fly -> 5 flies
#giving moles/5fly-second

if(use.mcmc == FALSE){
	
if(use.flight == FALSE){
if(use.line == TRUE){
	
gas_exchange = (line.resp[,c(1:2)]/1e7/22.4/3600)*5
}else{gas_exchange = (pop.resp[,c(1:2)]/1e7/22.4/3600)*5}
gas_exchange <- gas_exchange * SF
}

}else{
	gas_exchange = (line.resp[,c(1:2)]/1e7/22.4/3600)*5*SF
	#gas-change per mg tissue
	mcmc_list$Vcot <- (mcmc_list$Vcot/1e7/22.4/3600*5)/t(t(rep(1, times = n_mcmc_samples))) %*% t(line.enzyme$wts)*SF
	mcmc_list$Vox  <- (mcmc_list$Vox/1e7/22.4/3600*5)/t(t(rep(1, times = n_mcmc_samples))) %*% t(line.enzyme$wts)*SF
	
	}




nsamples <- length(kinetic_enzymes[,1])

if(use.flight == FALSE){

if(use.line == TRUE){
	save(kinetic_enzymes, kinetic_reactions, line.pop, line.color, pop.color, file = "color_scheme.R")
	}

if(ifelse(use.line == TRUE, !("line_resp_fba_fluxes.Rdata" %in% list.files()), !("pop_resp_fba_fluxes.Rdata" %in% list.files()))){

FBA_list = list()
for(mcmc_step in 1:n_mcmc_samples){

#sample matricies: kinetic_enzymes & gas_exchange
#iterate through all samples

calc.fluxes <- matrix(NA, nrow = length(S[1,]), ncol = nsamples)
rownames(calc.fluxes) <- colnames(S)
colnames(calc.fluxes) <- rownames(kinetic_enzymes)
optim.resid <- rep(NA, times = nsamples)

gas_exchange <- data.frame(Vcot = mcmc_list$Vcot[mcmc_step,], Vox = mcmc_list$Vox[mcmc_step,])

for (line in c(1:nsamples)){

#reactions corresponding to each measured Vmax

#G_flux <- matrix(data = NA, ncol = length(kinetic_enzymes[1,]), nrow = length(joint.stoi.red[1,])) 
#H_flux <- rep(NA, times = length(kinetic_enzymes[1,]))
#for (i in 1:length(kinetic_enzymes[1,])){
#	if(kinetic_reactions[i] %in% "Succinate dehydrogenase_m"){
#		G_flux[,i] <- ifelse(colnames(joint.stoi.red) == kinetic_reactions[i], -1, 0)
#		H_flux[i] <- -1*kinetic_enzymes[line,i]*10
#		}else{
#		G_flux[,i] <- ifelse(colnames(joint.stoi.red) == kinetic_reactions[i], -1, 0)
#		H_flux[i] <- -1*kinetic_enzymes[line,i]}
	#exception for vmax of reversible reactions
#	if(kinetic_reactions[i] %in% c("Malate dehydrogenase_m", "Phosphoglucomutase_c")){
#		G_flux <- cbind(G_flux, ifelse(colnames(joint.stoi.red) == kinetic_reactions[i], 1, 0))
#		H_flux <- c(H_flux, -1*kinetic_enzymes[line,i])
#	}}
#G_flux <- t(G_flux)

G_irr = matrix(data = NA, ncol = length(irreversible), nrow = length(joint.stoi.red[1,]))
for (i in 1:length(irreversible)){
G_irr[,i] = ifelse(colnames(joint.stoi.red) == irreversible[i], 1, 0)
	}
G_irr <- t(G_irr)
H_irr = rep(0, times = length(irreversible))

#G <- rbind(G_flux, G_irr)
#h <- c(H_flux, H_irr)
G <- G_irr
h <- H_irr


#reactions for CO2 and O2 exchange
exchange_rxns = c("CO2 leaving", "O2 entering")
A = matrix(data = NA, ncol = length(exchange_rxns), nrow = length(joint.stoi.red[1,]))
for (i in 1:length(exchange_rxns)){
	A[,i] <- ifelse(colnames(joint.stoi.red) == exchange_rxns[i], 1, 0)
	}
u = gas_exchange[line,c(1,2)]


QP.optim <- lsei(A = t(A), B = t(u), E = S, F = f, G = G, H = h)
if(QP.optim$solutionNorm < 0.01){
	calc.fluxes[,line] <- QP.optim$X
	}else{
		calc.fluxes[,line] <- NA
		}

optim.resid[line] <- QP.optim$solutionNorm

}
FBA_list[[mcmc_step]] <- calc.fluxes
}

if(use.line == TRUE){
	save(FBA_list, file = "line_resp_fba_fluxes.Rdata")
	}else{
		save(FBA_list, file = "pop_resp_fba_fluxes.Rdata")
		}
}else{
	if(use.line == TRUE){
		load("line_resp_fba_fluxes.Rdata")
		}else{
			load("pop_resp_fba_fluxes.Rdata")
			}
	}}




##### extract a median flux profile for each line #####

median.calc.fluxes <- matrix(NA, nrow = length(S[1,]), ncol = nsamples)
rownames(median.calc.fluxes) <- colnames(S)
colnames(median.calc.fluxes) <- rownames(kinetic_enzymes)

for(line in 1:nsamples){
	flux_samples <- sapply(1:n_mcmc_samples, function(x){FBA_list[[x]][,line]})	
	median.calc.fluxes[,line] <- apply(flux_samples, 1, median, na.rm = TRUE)
	}



valid.sample <- rep(TRUE, times = length(colnames(median.calc.fluxes)))
nonzero.flux <- median.calc.fluxes[apply(median.calc.fluxes[,valid.sample] == 0, 1, sum) != length(median.calc.fluxes[1, valid.sample]),]

if(use.line == TRUE){pdf(file = "line_flux.pdf")}else{pdf(file = "pop_flux.pdf")}

if(use.line == TRUE){
heatmap.2((nonzero.flux[,valid.sample] - apply(nonzero.flux[,valid.sample], 1, mean))/apply(nonzero.flux[,valid.sample], 1, sd), trace = "none", col = blue2yellow(100), ColSideColors = line.color[valid.sample], cexCol = 0.4, cexRow = 0.05 + 0.8/log10(length(nonzero.flux[,1])))

heatmap.2((nonzero.flux[,valid.sample])/apply(nonzero.flux[,valid.sample], 1, sd), trace = "none", col = blue2yellow(100), ColSideColors = line.color[valid.sample], cexCol = 0.4, cexRow = 0.05 + 0.8/log10(length(nonzero.flux[,1])))

	}else{
heatmap.2((nonzero.flux - apply(nonzero.flux, 1, mean))/apply(nonzero.flux, 1, sd), trace = "none", col = blue2yellow(100), cexRow = 0.05 + 0.8/log10(length(nonzero.flux[,1])))
	}
            
n.components <- 12                
#zero.flux <- colnames(S)[!nonzero.flux]
pc_corr <- data.frame(cotcorr = rep(NA, times = n.components), Otcorr = rep(NA, times = n.components), RQcorr = rep(NA, times = n.components))

if(use.line == TRUE){

pca_std <- nonzero.flux; pca_std <- t(scale(t(pca_std), center = TRUE, scale = TRUE))
gas_corr_samples <- gas_exchange[sapply(colnames(pca_std), function(x){c(1:length(gas_exchange[,1]))[rownames(gas_exchange) == x]}),]
gas_corr_samples <- cbind(gas_corr_samples, RQ = gas_corr_samples$Vcot/gas_corr_samples$Vox)
 
plot(((svd(pca_std, nu = 1, nv = 1)$d^2/sum(svd(pca_std, nu = 2, nv = 2)$d^2))), main = "Variance explained by each PC", ylab = "proportion of variance explained by PC")

for (pc in 1:n.components){
	
plot(svd(pca_std, nu = n.components, nv = n.components)$v[,pc], col = line.color, main = paste("Principle Component ", pc, sep = ""), xlab = paste("PC ", pc, sep = ""), ylab = "PC value", pch = 8, cex = 1.2)
legend("topleft", legend = pop.color[,3], text.col = as.character(pop.color[,2]))

pc_corr[pc,] <- apply(gas_corr_samples, 2, function(x){cor(x, svd(pca_std, nu = n.components, nv = n.components)$v[,pc])})

	}}


#use procrustes rotation to align principal components
library(shapes)
library(reshape)
library(ggplot2)
#use median as the reference and rotate all individual mcmc v samples to align against this

avgSVD <- svd(t(scale(t(nonzero.flux), center = TRUE, scale = TRUE)), nu = n.components, nv = n.components)$v
medPCloadings <- svd(t(scale(t(nonzero.flux), center = TRUE, scale = TRUE)), nu = n.components, nv = n.components)$u
colnames(medPCloadings) <- colnames(avgSVD) <- paste("PC", c(1:n.components), sep = "")
rownames(medPCloadings) <- rownames(nonzero.flux)

princ_compDF <- data.frame(sample = 0, avgSVD, pop = line.pop, size = 3, alpha = 1, stringsAsFactors = FALSE)
princ_compSave <- list()

for(mcmc_sample in 1:n_mcmc_samples){
	
	submat <- FBA_list[[mcmc_sample]][apply(FBA_list[[mcmc_sample]] != 0, 1, sum, na.rm = TRUE) != 0, apply(is.na(FBA_list[[mcmc_sample]]), 2, sum) == 0]
	
	sampleSVD <- svd(t(scale(t(submat), center = TRUE, scale = TRUE)), nu = n.components, nv = n.components)$v
	
	#rotate the sampleSVD to align against the median SVD
	rotSVD <- sampleSVD %*% procOPA(avgSVD[apply(is.na(FBA_list[[mcmc_sample]]), 2, sum) == 0,], sampleSVD, scale = FALSE, reflect = TRUE)$R
	
	#rotate the sampleSVD to align against the median SVD
	#rotSVD <- procOPA(avgSVD, sampleSVD, scale = FALSE, reflect = TRUE)$Bhat
	colnames(rotSVD) <- paste("PC", c(1:n.components), sep = "")
	rownames(rotSVD) <- colnames(submat)
  
  princ_compSave[[mcmc_sample]] <- rotSVD
  
  princ_compDF <- rbind(princ_compDF, data.frame(sample = mcmc_sample, rotSVD, pop = line.pop[apply(is.na(FBA_list[[mcmc_sample]]), 2, sum) == 0], size = 1, alpha = 0.08))
	
	}

save(princ_compDF, princ_compSave, file = "alignedPCs.Rdata")
#load("alignedPCs.Rdata") # if you dont want to wait

medianDF <- cbind(princ_compDF[princ_compDF$sample == 0,], shape = 21, fill = princ_compDF[princ_compDF$sample == 0,]$pop)
medianDF$pop <- NA 
princ_compDF <- cbind(princ_compDF[princ_compDF$sample != 0,], shape = 21, fill = princ_compDF[princ_compDF$sample != 0,]$pop)
princ_compDF <- rbind(princ_compDF, medianDF)

#medianDF <- princ_compDF[princ_compDF$sample %in% 2,]


princ_compDF_small <- princ_compDF[princ_compDF$sample <= 1000,]
princ_compDF_small$size <- princ_compDF_small$size*2

scatter_theme <- theme(text = element_text(size = 23, face = "bold", colour = "black"), title = element_text(colour = "black", size = 25, face = "bold"), panel.background = element_rect(fill = "white"), 
      legend.position = "right", panel.grid = element_blank(), axis.ticks = element_line(size = 0.2, colour = "black"), axis.text = element_text(color = "black")) 
      

princ_comp_plot <- ggplot(princ_compDF_small, aes(x = PC1, y = PC2, fill = factor(fill), colour = factor(pop), size = size, alpha = alpha, shape = shape)) + 
  scale_size_identity() + scale_alpha_identity()  + scale_x_continuous("Principal Component 1", color = "black", limits = c(-0.4, 0.4)) + scale_y_continuous("Principal Component 2", limits = c(-0.35, 0.5)) + scale_shape_identity() + 
  scale_color_brewer(palette = "Set2", na.value = "BLACK") + scale_fill_brewer(palette = "Set2", guide = 'none') + guides(colour = guide_legend(title = "Population"))
  
princ_comp_plot + geom_point() + scatter_theme +
  geom_vline(xintercept = c(-0.4, -0.2, 0, 0.2, 0.4), linetype = "longdash", size = 0.2, colour = "black") +
  geom_hline(yintercept = c(-0.4, -0.2, 0, 0.2, 0.4), linetype = "longdash", size = 0.2, colour = "black")
ggsave("principalComponents.pdf", height = 8, width = 10)




### associate flight with principal components



pc_flight_anova <- matrix(NA, ncol = length(mcmc_flight_cast[1,1,]), nrow = n_mcmc_samples)
colnames(pc_flight_anova) <- names(mcmc_flight_cast[1,1,])

mcmc_flight_cast_reorder <- mcmc_flight_cast[names(mcmc_flight_cast[,1,1]) %in% colnames(nonzero.flux),,]
mcmc_flight_cast_reorder <- mcmc_flight_cast_reorder[chmatch(colnames(nonzero.flux), names(mcmc_flight_cast_reorder[,1,1])),,]

#sort(names(mcmc_flight_cast[,i,1])[names(mcmc_flight_cast[,i,1]) %in% colnames(nonzero.flux)]) == colnames(nonzero.flux)]

for(i in 1:n_mcmc_samples){
  matchedTraits <- mcmc_flight_cast_reorder[,i,][rownames(mcmc_flight_cast_reorder[,i,]) %in% rownames(princ_compSave[[i]]),]
  if(nrow(matchedTraits) < 50){next}
  
  design_matrix <- model.matrix( ~ PC1 + PC2 + PC1*PC2 + 0, as.data.frame(princ_compSave[[i]]))
  for(flight_pheno in 1:length(mcmc_flight_cast[1,1,])){
    pc_flight_anova[i, flight_pheno] <- anova(lm(matchedTraits[,flight_pheno] ~ design_matrix))[1,5]
    }
  }

atp_flight_anova <- matrix(NA, ncol = length(mcmc_flight_cast[1,1,]), nrow = n_mcmc_samples)
colnames(atp_flight_anova) <- names(mcmc_flight_cast[1,1,])



energyGen <- data.frame(reaction = c("Hexokinase_c", "Phosphofructokinase_c", "Phosphoglycerate kinase_c", "Pyruvate kinase_c", "ATP synthetase"), ATPchange = c(-1, -1, 1, 1, 1))

for(i in 1:n_mcmc_samples){
  mcmc_flight_cast_reorder[,i,]
  flux_dist <- FBA_list[[mcmc_sample]][rownames(FBA_list[[mcmc_sample]]) %in% energyGen$reaction, colnames(FBA_list[[mcmc_sample]]) %in% names(mcmc_flight_cast_reorder[,1,1])]
    rownames(FBA_list[[i]]) %in% rownames(mcmc_flight_cast_reorder[,i,]),]
  
  design_matrix <- model.matrix( ~ PC1 + PC2 + PC1*PC2 + 0, as.data.frame(princ_compSave[[i]]))
  for(flight_pheno in 1:length(mcmc_flight_cast[1,1,])){
    pc_flight_anova[i, flight_pheno] <- anova(lm(matchedTraits[,flight_pheno] ~ design_matrix))[1,5]
    }
  }



#### for measured enzymes, determine v / vmax

od.measured.carried.flux <- matrix(NA, ncol = length(kinetic_reactions), nrow = length(kinetic_enzymes[,1]))

for(rxn in 1:length(kinetic_reactions)){
	if(sum(rownames(median.calc.fluxes) %in% kinetic_reactions[rxn]) != 0){
	od.measured.carried.flux[,rxn] <- median.calc.fluxes[rownames(median.calc.fluxes) %in% kinetic_reactions[rxn],]
	}else{
		od.measured.carried.flux[,rxn] <- 0
		}
	}	

colnames(od.measured.carried.flux) <- kinetic_reactions; rownames(od.measured.carried.flux) <- colnames(median.calc.fluxes)

kinetic_enzymes[kinetic_enzymes == 0] <- NA

Vmax.fraction <- od.measured.carried.flux/kinetic_enzymes
#measure correlations between carried flux and vmax
corr_v_vmax <- sapply(c(1:length(od.measured.carried.flux[1,])), function(x){
	if(sd(od.measured.carried.flux[,x]) != 0){
    carried_max <- cbind(od.measured.carried.flux[,x], kinetic_enzymes[,x])
    carried_max <- carried_max[apply(is.na(carried_max), 1, sum) == 0,]
    carried_max <- carried_max[carried_max[,1] != 0,]
    
		corval <- signif(cor(abs(carried_max[,1]), carried_max[,2], method = "spearman"), 3)
		pcorval <- signif(cor(abs(carried_max[,1]), carried_max[,2], method = "pearson"), 3)
    r_square <- signif(cor(abs(carried_max[,1]), carried_max[,2], method = "pearson")^2, 3)
    print(plot(abs(carried_max[,1]) ~ carried_max[,2], main = paste(colnames(od.measured.carried.flux)[x], "s", corval, "p2", r_square, sep = ": ")))
    c(corval, pcorval, r_square, length(carried_max[,1]))
    }else{
			c(NA, NA, NA, NA)
			}
	})

flux_corrTab <- cbind(t(corr_v_vmax), abs(apply(Vmax.fraction, 2, function(x){median(x[x != 0],na.rm = TRUE)})))
flux_corrTab <- cbind(c("Alcohol Dehydrogenase", "Malic Enzyme", "Malate Dehydrogenase (cytosolic)", "Fatty Acid Synthase", "Phosphoglucose Isomerase", "Phosphofructokinase", "Glucose 6-Phosphate Dehydrogenase", "6-Phosphogluconate Dehydrogenase", "Phosphoglucomutase", "Glycogen Phosphorylase", "Glycogen Synthase", "Trehalase", "Glucokinase", "Glycerol 3-Phosphate Dehydrogenase (cytosolic)", "Glycerol 3-Phosphate Dehydrogenase (mitochondrial)", "Pyruvate Dehydrogenase", "Fumarase", "Malate Dehydrogenase (mitochondrial)", "Succinate Dehydrogenase"), flux_corrTab)
colnames(flux_corrTab) <- c("Enzyme", "Spearman", "Pearson", "rSquared", "Nlines", "vmaxFrac")
flux_corrTab <- as.data.frame(flux_corrTab, stringsAsFactors = FALSE)
flux_corrTab[,-1] <- apply(flux_corrTab[,-1], c(1,2), as.numeric)
flux_corrTab <- flux_corrTab[!is.na(flux_corrTab$vmaxFrac),]
#flux_corrTab[flux_corrTab$Nlines == length(Vmax.fraction[,1]),]
flux_corrTab[flux_corrTab$Nlines != length(Vmax.fraction[,1]),c(2:4)] <- NA
flux_corrTab <- flux_corrTab[,colnames(flux_corrTab) != "Nlines"]          
flux_corrTab <- flux_corrTab[!(flux_corrTab$Enzyme %in% c("Glycogen Phosphorylase", "Glycogen Synthase", "Trehalase")),]
flux_corrTab <- flux_corrTab[order(flux_corrTab$rSquared, decreasing = TRUE),]
  
library(xtable)
rownames(flux_corrTab) <- NULL
print(xtable(flux_corrTab, display = c("s", "s", "g", "g", "g", "g"), digits = 3), include.rownames = FALSE)


#spearman correlations between Vmax and V
#PDH: 0.21, vmax frac = 0.0157, r2 = 0.0568
#G3P shuttle (GPO): 0.27, vmax frac = 0.117, r2 = 0.109
#G3P dehydorgenase (GPDH): 0.56, vmax frac = 0.00409, r2 = 0.322
#PFK: 0.23, vmax frac  = 0.0437, r2 = 0.055
#FUM: s: -0.141, r2 = 0.024, vmax frac = 0.0219
library(xtable)
flux_corrTab <- data.frame(Enzyme = c("PDH", "GPO", "GPDH", "PFK", "FUM"), spearCorr = c(0.21, 0.27, 0.558, 0.231, -0.141), r_squared = c(0.0568, 0.109, 0.322, 0.055, 0.024), "vmaxfrac" = c(0.0157, 0.117, 0.00409, 0.0437, 0.0219)) 





#measure correlations between vmax_frac and vmax
corr_vfrac_vmax <- sapply(c(1:length(Vmax.fraction[1,])), function(x){
	if(sd(Vmax.fraction[,x], na.rm = TRUE) != 0){
		cor(abs(Vmax.fraction[,x])[!is.na(Vmax.fraction[,x])], kinetic_enzymes[,x][!is.na(Vmax.fraction[,x])], method = "spearman")
		}else{
			NA
			}
	})

plot(corr_v_vmax ~ corr_vfrac_vmax)

plot(od.measured.carried.flux[,x] ~  kinetic_enzymes[,x], pch = 16, col = line.color)

vmax.frac.mat <- scale(Vmax.fraction, center = FALSE, scale = TRUE); vmax.frac.mat <- abs(vmax.frac.mat[,!is.nan(vmax.frac.mat[1,])])
vmax.frac.mat <- vmax.frac.mat[rownames(vmax.frac.mat) != "N13",]

#heatmap.2(vmax.frac.mat, trace = "n")

if(use.line == TRUE){
	write.table(median.calc.fluxes, file = "line_fluxes.tsv", sep = "\t")
	write.table(Vmax.fraction, file = "line_vmax_fraction.tsv", sep = "\t")
	} else {
		write.table(median.calc.fluxes, file = "pop_fluxes.tsv", sep = "\t")
		write.table(Vmax.fraction, file = "pop_vmax_fraction.tsv", sep = "\t")
		}
		

#Get gas-exchange for PCs

medPCloadings[rownames(medPCloadings) %in% c("CO2 leaving", "O2 entering"),][1,]/medPCloadings[rownames(medPCloadings) %in% c("CO2 leaving", "O2 entering"),][2,]




library(gplots)
n_colors = 1001
PCgrad_colors <- colorpanel(n_colors, "blue1", "antiquewhite", "gold1")
#plot(1:n_colors,1:n_colors, col = PCgrad_colors)
  
PCflux <- data.frame(reaction = names(medPCloadings[,1]), dirflux = medPCloadings[,1], invflux = medPCloadings[,1]*-1, stringsAsFactors = FALSE)

PCflux$dirColor = PCgrad_colors[round((PCflux[,2] + max(abs(PCflux[,2])))/(2*max(abs(PCflux[,2]))) * (n_colors/(n_colors+1))*n_colors) + 1]
PCflux$invColor = PCgrad_colors[round((PCflux[,3] + max(abs(PCflux[,3])))/(2*max(abs(PCflux[,3]))) * (n_colors/(n_colors+1))*n_colors) + 1]

PC_colormap <- list()
PC_colormap[[1]] <- PCflux
image(z=matrix(1:n_colors, ncol=1), x = seq(-1*max(abs(PCflux[,2])), max(abs(PCflux[,2])), by = 2*max(abs(PCflux[,2]))/(n_colors*(n_colors/(n_colors+1)))), col=PCgrad_colors, yaxt="n", xlab = "")



PCflux <- data.frame(reaction = names(medPCloadings[,2]), dirflux = medPCloadings[,2], invflux = medPCloadings[,2]*-1, stringsAsFactors = FALSE)

PCflux$dirColor = PCgrad_colors[round((PCflux[,2] + max(abs(PCflux[,2])))/(2*max(abs(PCflux[,2]))) * (n_colors/(n_colors+1))*n_colors) + 1]
PCflux$invColor = PCgrad_colors[round((PCflux[,3] + max(abs(PCflux[,3])))/(2*max(abs(PCflux[,3]))) * (n_colors/(n_colors+1))*n_colors) + 1]

PC_colormap[[2]] <- PCflux
image(z=matrix(1:n_colors, ncol=1), x = seq(-1*max(abs(PCflux[,2])), max(abs(PCflux[,2])), by = 2*max(abs(PCflux[,2]))/(n_colors*(n_colors/(n_colors+1)))), col=PCgrad_colors, yaxt="n", xlab = "")


PCflux <- data.frame(reaction = names(medPCloadings[,3]), dirflux = medPCloadings[,3], invflux = medPCloadings[,3]*-1, stringsAsFactors = FALSE)

PCflux$dirColor = PCgrad_colors[round((PCflux[,2] + max(abs(PCflux[,2])))/(2*max(abs(PCflux[,2]))) * (n_colors/(n_colors+1))*n_colors) + 1]
PCflux$invColor = PCgrad_colors[round((PCflux[,3] + max(abs(PCflux[,3])))/(2*max(abs(PCflux[,3]))) * (n_colors/(n_colors+1))*n_colors) + 1]

PC_colormap[[3]] <- PCflux
image(z=matrix(1:n_colors, ncol=1), x = seq(-1*max(abs(PCflux[,2])), max(abs(PCflux[,2])), by = 2*max(abs(PCflux[,2]))/(n_colors*(n_colors/(n_colors+1)))), col=PCgrad_colors, yaxt="n", xlab = "")


save(PC_colormap, file = "PC_loadingsColors.Rdata")







		
		
		
		
