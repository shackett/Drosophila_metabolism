###### Read in absorbance posterior medians for enzyme activity traits as well as flight and respirometry traits

setwd("/Users/seanhackett/Desktop/Cornell/Drosophila_metabolism/MCMC_files")

line.enzyme <- read.table("lnMat.tsv", sep = "\t", header = TRUE)
line.flight <- read.table("FLTlnMns.tsv", sep = "\t", header = TRUE)
line.resp <- read.table("RESPlnMns.tsv", sep = "\t", header = TRUE)

pop.enzyme <- read.table("popMat.tsv", sep = "\t", header = TRUE)
pop.flight <- read.table("FLTpopMns.tsv", sep = "\t", header = TRUE)
pop.resp <- read.table("RESPpopMns.tsv", sep = "\t", header = TRUE)
rownames(pop.enzyme) <- as.character(pop.enzyme$pop_nam)
rownames(pop.flight) <- as.character(pop.flight$pop_nam)
rownames(pop.resp) <- as.character(pop.resp$pop_nam)
pop.enzyme <- pop.enzyme[,-1]
pop.flight <- pop.flight[,-1]
pop.resp <- pop.resp[,-1]

setwd("/Users/seanhackett/Desktop/Cornell/Drosophila_metabolism/")

load("drosophila_stoi.R")
source("path_length_calc.R")

###### Read in and calculate the paramters needed to transform A/s into moles/s using the beer-lambert law and path-length calculated from the well-geometry and volume of a microtiter-assay

assay_parameters <- read.table("kin_par.csv", sep = ",", header = TRUE, colClasses = c("character", "numeric", "character", "numeric", "numeric", "character"))
rownames(assay_parameters) <- assay_parameters$Enzyme
assay_parameters <- assay_parameters[,-1]
assay_parameters <- cbind(assay_parameters, data.frame(path_length = NA))

for(i in 1:length(assay_parameters[,1])){
assay_parameters$path_length[i] <- find.pl(assay_parameters$assay_volume[i], height_scale, bottom_d)$minimum * 100	}
	
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

#maximum flux through 

library(limSolve)

joint.stoi.store -> joint.stoi

##### add all metabolites in as boundary

#joint.stoi <- cbind(joint.stoi, diag(length(joint.stoi[,1])))
#valid.solution <- rep(NA, times = length(joint.stoi[,1]))

#for (h in 1:length(joint.stoi[,1])){

joint.stoi.red <- joint.stoi#cbind(joint.stoi, diag(length(joint.stoi[,1]))[,-c(1:80)])

	
E = joint.stoi.red
F = rep(0, times = length(joint.stoi.red[,1]))

#reactions corresponding to each measured Vmax

irreversible = c("Trehalose6P synthetase_c", "Glycogen Synthase_c", "6-phosphogluconate dehydrogenase_c", "Phosphoenolpyruvate carboxylase_c", "Pyrophosphatase_c", "Malic enzyme (NADP)_c", "Malic enzyme (NAD)_c", "C8 synthesis_c", "Pyruvate dehydrogenase_m", "Isocitrate dehydrogenase (NAD)_m", "Isocitrate dehydrogenase (NADP)_m", "Alpha-ketoglutarate dehydrogenase I_m", "ATP synthetase", "ATPase_c")

G = matrix(data = NA, ncol = length(irreversible), nrow = length(joint.stoi.red[1,]))

for (i in 1:length(irreversible)){
G = 
	}

#reactions for CO2 and O2 exchange

exchange_rxns = c("CO2 leaving", "O2 entering")
A = matrix(data = NA, ncol = length(exchange_rxns), nrow = length(joint.stoi.red[1,]))
for (i in 1:length(exchange_rxns)){
	A[,i] <- ifelse(colnames(joint.stoi.red) == exchange_rxns[i], 1, 0)
	}
A <- t(A)	
B = pop.resp[1,c(1,2)]


lsei(A = A, B = t(B), E = E, F = F)
solution <- lsei(A = A, B = t(B), E = E, F = F)$X

	
	
	
	
	


#colnames(joint.stoi)[!(colnames(joint.stoi) %in% reaction_names)]

#reduce joint.stoi
#reaction_names = c("CO2 leaving", "O2 entering", "6-phosphogluconate dehydrogenase_c", "Glucose 6-phosphate dehydrogenase_c", "Lactonase_c", "Ribulose 5-phosphate epimerase_c", "Ribose 5-phosphate isomerase_c", "Transaldolase: 6,4-7,3_c", "Transketolase: 7,3-5,5_c", "Transketolase: 6,4-5,5_c", "Phosphoglucose isomerase_c", "Trehalose6P phosphatase_c", "Trehalase_c", "Trehalose6P synthetase_c", "Hexokinase_c", "Phosphoglucomutase_c", "Branching enzyme_c", "Glycogen phosphorylase_c", "Glycogen Synthase_c", "UDP-pyrophosphorylase_c") 
#"Phosphofructokinase_c")

#"Aldolase_c", "Triose phosphate isomerase_c", "Glyceraldehyde 3P dehydrogenase_c", "Phosphoglycerate kinase_c", "Phosphoglycerate mutase_c", "Enolase_c", "Glycerol 3-phosphate dehydrogenase (NAD)_c"

#joint.stoi <- joint.stoi[,colnames(joint.stoi) %in% reaction_names]
#joint.stoi <- joint.stoi[apply(joint.stoi != 0, 1, sum) != 0,]
#added_boundary <- c("H2O_c", "H+_c", "Nicotinamide adenine dinucleotide phosphate - reduced_c", "Nicotinamide adenine dinucleotide phosphate_c", "UTP_c", "UDP_c", "Phosphate_c", "ATP_c", "ADP_c", "Nicotinamide adenine dinucleotide_c", "Nicotinamide adenine dinucleotide - reduced_c ", "Phosphoenolpyruvate_c", "Glycerol 3-phosphate_c", "D-Fructose 2,6-bisphosphate_c", "D-Glucose_c", "D-Glucose 6-phosphate_c")
#for(i in 1:length(added_boundary)){
#	joint.stoi <- cbind(joint.stoi, ifelse(rownames(joint.stoi) == added_boundary[i], 1, 0))
#	}
#colnames(joint.stoi)[colnames(joint.stoi) == ""] <- added_boundary


#linp
E = joint.stoi
F = rep(0, times = length(joint.stoi[,1]))

#reactions corresponding to each measured Vmax

#reactions for CO2 and O2 exchange

exchange_rxns = c("CO2 leaving", "O2 entering")
A = matrix(data = NA, ncol = length(exchange_rxns), nrow = length(joint.stoi[1,]))
for (i in 1:length(exchange_rxns)){
	A[,i] <- ifelse(colnames(joint.stoi) == exchange_rxns[i], 1, 0)
	}
A <- t(A)	
B = pop.resp[1,c(1,2)]


solution <- lsei(A = A, B = t(B), E = E, F = F)$X

xblah = rep(1, times = length(joint.stoi[1,]))

(A %*% xblah) == B
(E %*% xblah) == F

sum(((A %*% solution) - B)^2)
E %*% solution == F

dim(A)
length(B)

dim(E)
length(F)

dim(E_coli$A)
length(E_coli$B)

dim(E_coli$G)
length(E_coli$H)






                               












find.name <- function(kinetic_enzyme, assay_parameters){
	assay_parameters$Model_name[rownames(assay_parameters) == kinetic_enzyme]
	}




 
 

