###### Read in absorbance posterior medians for enzyme activity traits as well as flight and respirometry traits

setwd("/Users/seanhackett/Desktop/Cornell/Drosophila_metabolism/")

line.enzyme <- read.table("MCMC_files/lnMat.tsv", sep = "\t", header = TRUE)
line.flight <- read.table("MCMC_files/FLTlnMns.tsv", sep = "\t", header = TRUE)
line.resp <- read.table("MCMC_files/RESPlnMns.tsv", sep = "\t", header = TRUE)
rownames(line.enzyme) <- as.character(line.enzyme$ln_nam)
rownames(line.flight) <- as.character(line.flight$ln_nam)
rownames(line.resp) <- as.character(line.resp$ln_nam)
line.enzyme <- line.enzyme[,-1]; line.flight <- line.flight[,-1]; line.resp <- line.resp[,-1]

pop.enzyme <- read.table("MCMC_files/popMat.tsv", sep = "\t", header = TRUE)
pop.flight <- read.table("MCMC_files/FLTpopMns.tsv", sep = "\t", header = TRUE)
pop.resp <- read.table("MCMC_files/RESPpopMns.tsv", sep = "\t", header = TRUE)
rownames(pop.enzyme) <- as.character(pop.enzyme$pop_nam)
rownames(pop.flight) <- as.character(pop.flight$pop_nam)
rownames(pop.resp) <- as.character(pop.resp$pop_nam)
pop.enzyme <- pop.enzyme[,-1]; pop.flight <- pop.flight[,-1]; pop.resp <- pop.resp[,-1]

load("drosophila_stoi.R")
source("path_length_calc.R")
source("flimsy_plate_extract.R")
library(limSolve)
library(gplots)

##### a few functions

find.name <- function(kinetic_enzyme, assay_parameters){
	assay_parameters$Model_name[rownames(assay_parameters) == kinetic_enzyme]
	}



###### Read in and calculate the paramters needed to transform A/s into moles/s using the beer-lambert law and path-length calculated from the well-geometry and volume of a microtiter-assay

assay_parameters <- read.table("kin_par.csv", sep = ",", header = TRUE, colClasses = c("character", "numeric", "character", "numeric", "numeric", "character"))
rownames(assay_parameters) <- assay_parameters$Enzyme
assay_parameters <- assay_parameters[,-1]
assay_parameters <- cbind(assay_parameters, data.frame(path_length = NA))

for(i in 1:length(assay_parameters[,1])){

#If assay is visible then use a different type of plate is used
if(assay_parameters$molar_absorptivity[i] %in% c(16000, 19100, 11500)){
	assay_parameters$path_length[i] <- flimsy_pl(assay_parameters$assay_volume[i]*1e6)
	}else{
assay_parameters$path_length[i] <- find.pl(assay_parameters$assay_volume[i], height_scale, bottom_d)$minimum * 100	}}
	
enzymeOD <- pop.enzyme[,((colnames(pop.enzyme) %in% rownames(assay_parameters)))]
enzyme_moles_per_second <- enzymeOD
enzyme_moles_per_second[!is.na(enzyme_moles_per_second)] <- NA
	
for (i in 1:length(enzyme_moles_per_second[1,])){
	enzyme <- colnames(enzyme_moles_per_second)[i]
	params <- assay_parameters[rownames(assay_parameters) == enzyme,]
	
	if(params$molar_absorptivity != "STD"){
	params$molar_absorptivity <- as.numeric(params$molar_absorptivity)
	
	enzyme_moles_per_second[,i] <- ((enzymeOD[,i]/(params$molar_absorptivity*params$path_length))*params$assay_volume*(params$assay_volume/20e-6))/params$OD_scaling
	}}

kinetic_enzymes = enzyme_moles_per_second[,apply(is.na(enzyme_moles_per_second), 2, sum) == 0]
kinetic_reactions <- unlist(lapply(colnames(kinetic_enzymes), find.name, assay_parameters)) 

#kinetic_enzymes <- kinetic_enzymes[,colnames(kinetic_enzymes) != "SDH"]
#kinetic_reactions <- kinetic_reactions[colnames(kinetic_enzymes) != "SDH"]

SF <- 1e11

gas_exchange = pop.resp[,c(1:2)]/1e7/22.4/3600

kinetic_enzymes <- kinetic_enzymes * SF
gas_exchange <- gas_exchange * SF

joint.stoi.red <- joint.stoi
	
S = joint.stoi.red
f = rep(0, times = length(joint.stoi.red[,1]))

# irreversible reactions
irreversible = c("Trehalose6P synthetase_c", "Glycogen Synthase_c", "6-phosphogluconate dehydrogenase_c", "Phosphoenolpyruvate carboxylase_c", "Pyrophosphatase_c", "FADH2 oxidation for proton transport", "NADH oxidation for proton transport", "Isocitrate dehydrogenase (NADP)_m", "Isocitrate dehydrogenase (NADP)_m", "Alpha-ketoglutarate dehydrogenase I_m", "Trehalase_c", "Malic enzyme (NADP)_c", "Malic enzyme (NAD)_c", "C8 synthesis_c", "Pyruvate dehydrogenase_m", "Phosphofructokinase_c", "Trehalase_c")
#, "Phosphoglycerate kinase_c", 

#sample matricies: kinetic_enzymes & gas_exchange

#iterate through all samples

pop.fluxes <- matrix(NA, nrow = length(S[1,]), ncol = 5)
rownames(pop.fluxes) <- colnames(S)
colnames(pop.fluxes) <- rownames(kinetic_enzymes)

for (line in c(1:5)){

#reactions corresponding to each measured Vmax

G_flux <- matrix(data = NA, ncol = length(kinetic_enzymes[1,]), nrow = length(joint.stoi.red[1,])) 
H_flux <- rep(NA, times = length(kinetic_enzymes[1,]))
for (i in 1:length(kinetic_enzymes[1,])){
	G_flux[,i] <- ifelse(colnames(joint.stoi.red) == kinetic_reactions[i], -1, 0)
	H_flux[i] <- -1*kinetic_enzymes[line,i]
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


#reactions for CO2 and O2 exchange

exchange_rxns = c("CO2 leaving", "O2 entering")
A = matrix(data = NA, ncol = length(exchange_rxns), nrow = length(joint.stoi.red[1,]))
for (i in 1:length(exchange_rxns)){
	A[,i] <- ifelse(colnames(joint.stoi.red) == exchange_rxns[i], 1, 0)
	}
u = gas_exchange[line,c(1,2)]


pop.fluxes[,line] <- lsei(A = t(A), B = t(u), E = S, F = f, G = G, H = h)$X
}

nonzero.pop.flux <- pop.fluxes[apply(pop.fluxes == 0, 1, sum) != length(pop.fluxes[1,]),]
heatmap.2(nonzero.pop.flux - apply(nonzero.pop.flux, 1, mean))/apply(nonzero.pop.flux, 1, sd)




                               

















 
 

