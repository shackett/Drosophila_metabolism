setwd("/Users/seanhackett/Desktop/Cornell/Drosophila_metabolism/")


#functions to find reactions in the yeast metabolomic reconstruction
load("yeast_stoi.R")
#rxns = rxn_search(stoiMat, "Acetyl-CoA synthase", is_rxn = TRUE)
#cols <- c(232, 596, 816)
#redmat <- stoiMat[apply(stoiMat[,cols], 1, is.not.zero), cols]

valid.reactions <- read.table("mel_network.csv", sep = ",", blank.lines.skip = FALSE, header = TRUE)

cytosolic.reactions <- valid.reactions[valid.reactions$Compartment == "c",][!is.na(valid.reactions[valid.reactions$Compartment == "c",]$ReactionNumber),]
cytosolic.stoi <- stoiMat[,cytosolic.reactions$ReactionNumber][apply(stoiMat[,cytosolic.reactions$ReactionNumber], 1, is.not.zero),]
rownames(cytosolic.stoi) <- lapply(rownames(cytosolic.stoi), paste, "_c", sep = "")
colnames(cytosolic.stoi) <- paste(cytosolic.reactions$Name, cytosolic.reactions$Compartment, sep = "_")

#Hrep <- 

#cytosolic.stoi[rownames(cytosolic.stoi) == "H+_c",] == 0

#cytosolic.stoi[rownames(cytosolic.stoi) == "H+_c",] <- 



mitochondrial.reactions <- valid.reactions[valid.reactions$Compartment == "m",][!is.na(valid.reactions[valid.reactions$Compartment == "m",]$ReactionNumber),]
mitochondrial.stoi <- stoiMat[, mitochondrial.reactions $ReactionNumber][apply(stoiMat[, mitochondrial.reactions $ReactionNumber], 1, is.not.zero),]
rownames(mitochondrial.stoi) <- lapply(rownames(mitochondrial.stoi), paste, "_m", sep = "")
colnames(mitochondrial.stoi) <- paste(mitochondrial.reactions$Name, mitochondrial.reactions$Compartment, sep = "_")

#Hand-added reactions

cytosolic.stoi <- add.reaction(cytosolic.stoi, c("Diphosphate_c", "H2O_c"), c(1,1), "Phosphate_c", 2, "Pyrophosphatase_c")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("H2O_c", "Coenzyme A_c", "ATP_c", "Citrate_c") , c(1,1,1,1), 
c("Oxaloacetate_c", "Acetyl-CoA_c", "ADP_c", "Phosphate_c"), c(1,1,1,1), "Citrate lyase_c")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("Carnitine_c", "Palmitoyl-CoA (n-C16:0CoA)_c"), c(1,1), c("Palmitoylcarnitine_c", "Coenzyme A_c"), c(1,1), "Palmitoylcarnitine synthesis_c")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("UDP_c", "ATP_c"), c(1,1), c("UTP_c", "ADP_c"), c(1,1), "Nucloside diphosphate kinase (UTP-ATP)_c")
cytosolic.stoi <- add.reaction(cytosolic.stoi, "H2O_c", 1, c("OH-_c", "H+_c"), c(1,1), "Water dissociation_c")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("H2O_c", "ATP_c"), c(1,1), c("Phosphate_c", "ADP_c", "H+_c"), c(1,1,1), "ATPase_c")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("H+_c", "Nicotinamide adenine dinucleotide phosphate - reduced_c", "Acetyl-ACP_c", "Malonyl-[acyl-carrier protein]_c"), c(21, 14, 1, 7), c("H2O_c", "CO2_c", "Nicotinamide adenine dinucleotide phosphate_c", "acyl carrier protein_c", "Palmitoyl-ACP (n-C16:0ACP)_c"), c(7, 7, 14, 7, 1), "palmitoyl-ACP synthesis")



mitochondrial.stoi <- add.reaction(mitochondrial.stoi, c("Palmitoylcarnitine_m", "Coenzyme A_m"), c(1,1), c("Carnitine_m", "Palmitoyl-CoA (n-C16:0CoA)_m"), c(1,1), "Palmitoylcarnitine breakdown_m")
mitochondrial.stoi <- add.reaction(mitochondrial.stoi, 
c("Palmitoyl-CoA (n-C16:0CoA)_m", "Flavin adenine dinucleotide oxidized_m", "H2O_m", "Nicotinamide adenine dinucleotide_m", "Coenzyme A_m"),
c(1, 7, 7, 7, 7), c("Flavin adenine dinucleotide reduced_m", "Nicotinamide adenine dinucleotide - reduced_m", "H+_m", "Acetyl-CoA_m"), c(7, 7, 7, 8),
"Palmitoyl-CoA catabolism_m")
mitochondrial.stoi <- add.reaction(mitochondrial.stoi, c("Succinyl-CoA_m", "Phosphate_m", "GDP_m"), c(1,1,1), c("Succinate_m", "Coenzyme A_m", "GTP_m"), c(1,1,1), "Succinyl CoA synthetase_m")
mitochondrial.stoi <- add.reaction(mitochondrial.stoi, "H2O_m", 1, c("OH-_m", "H+_m"), c(1,1), "Water dissociation_m")
mitochondrial.stoi <- add.reaction(mitochondrial.stoi, c("H2O_m", "CO2_m"), c(1,1), c("H+_m", "Bicarbonate_m"), c(1,1), "Bicarbonate equilibrium_m")

joint.stoi <- cbind(rbind(cytosolic.stoi, matrix(data = 0, ncol = length(cytosolic.stoi[1,]), nrow = length(mitochondrial.stoi[,1]))), rbind(matrix(data = 0, ncol = length(mitochondrial.stoi[1,]), nrow = length(cytosolic.stoi[,1])), mitochondrial.stoi))
rownames(joint.stoi) <- c(rownames(cytosolic.stoi), rownames(mitochondrial.stoi))
colnames(joint.stoi) <- c(colnames(cytosolic.stoi), colnames(mitochondrial.stoi))

joint.stoi <- add.reaction(joint.stoi, "Trehalose_c", 1, NULL, NULL, "Trehalose usage")
joint.stoi <- add.reaction(joint.stoi, "glycogen_c", 1, NULL, NULL, "Glycogen usage")
joint.stoi <- add.reaction(joint.stoi, "Hexadecanoate (n-C16:0)_c", 1, NULL, NULL, "Palmitate usage")
#joint.stoi <- add.reaction(joint.stoi, "alpha-D-Ribose 5-phosphate_c", 1, NULL, NULL, "Ribose biosynthesis")
joint.stoi <- add.reaction(joint.stoi, "CO2_c", 1, NULL, NULL , "CO2 leaving")
joint.stoi <- add.reaction(joint.stoi, NULL, NULL, "O2_c", 1, "O2 entering")

joint.stoi <- add.reaction(joint.stoi, c("Glycerol 3-phosphate_c", "Flavin adenine dinucleotide oxidized_m"), c(1,1), 
c("Dihydroxyacetone phosphate_c", "Flavin adenine dinucleotide reduced_m"), c(1,1), "Glycerol 3-phosphate shuttle")
joint.stoi <- add.reaction(joint.stoi, c("Flavin adenine dinucleotide reduced_m", "H+_m", "O2_m"), c(1, 8, 0.5), c("Flavin adenine dinucleotide oxidized_m", "H+_c", "H2O_m"), c(1, 6, 1), "FADH2 oxidation for proton transport")
joint.stoi <- add.reaction(joint.stoi, c("Nicotinamide adenine dinucleotide - reduced_m", "H+_m", "O2_m"), c(1, 12, 0.5), c("Nicotinamide adenine dinucleotide_m", "H+_c", "H2O_m"), c(1, 10, 1), "NADH oxidation for proton transport")
joint.stoi <- add.reaction(joint.stoi, c("ADP_m", "Phosphate_m", "H+_c"), c(1, 1, 3), c("ATP_m", "H+_m"), c(1, 3), "ATP synthetase")
joint.stoi <- add.reaction(joint.stoi, c("GTP_m", "ADP_m"), c(1,1), c("GDP_m", "ATP_m"), c(1,1), "Nucleoside diphosphate kinase (GTP-ATP)")
joint.stoi <- add.reaction(joint.stoi, c("Pyruvate_c", "H+_c"), c(1,1), c("Pyruvate_m", "H+_m"), c(1,1), "Pyruvate transporter")
joint.stoi <- add.reaction(joint.stoi, c("Phosphate_c", "H+_c"), c(1,1), c("Phosphate_m", "H+_m"), c(1,1), "Phosphate transporter")
joint.stoi <- add.reaction(joint.stoi, c("ADP_c", "ATP_m"), c(1,1), c("ATP_c", "ADP_m"), c(1,1), "ATP/ADP antiport")
joint.stoi <- add.reaction(joint.stoi, c("L-Malate_m", "Phosphate_c"), c(1,1), c("L-Malate_c", "Phosphate_m"), c(1,1), "Dicarboxylate carrier")
joint.stoi <- add.reaction(joint.stoi, c("Citrate_m", "H+_m", "L-Malate_c"), c(1,1,1), c("Citrate_c", "H+_c", "L-Malate_m"), c(1,1,1), "Tricarboxylate carrier")
joint.stoi <- add.reaction(joint.stoi, "H2O_c", 1, "H2O_m", 1, "Water transport")
joint.stoi <- add.reaction(joint.stoi, "O2_c", 1, "O2_m", 1, "Oxygen transport")
joint.stoi <- add.reaction(joint.stoi, "CO2_m", 1, "CO2_c", 1, "Carbon dioxide transport")
joint.stoi <- add.reaction(joint.stoi, "H2O_c", 1, NULL, NULL, "Water boundary")
joint.stoi <- add.reaction(joint.stoi, "Palmitoylcarnitine_c", 1, "Palmitoylcarnitine_m", 1, "palmitoylcarnitine transport")


compartment.mets <- list()
compartment.mets[["c"]] <- c(rownames(cytosolic.stoi), "O2_c")
compartment.mets[["m"]] <- c(rownames(mitochondrial.stoi), "O2_m")


#joint.stoi <- add.reaction(joint.stoi, "2-Oxoglutarate_m", 1, NULL, NULL, "ketoglutarate boundary")
#joint.stoi <- add.reaction(joint.stoi, c("L-Malate_m", "H+_m"), c(1,1), c("L-Malate_c", "H+_c"), c(1,1), "Malate transporter")
#joint.stoi <- add.reaction(joint.stoi, c("Citrate_m", "H+_m"), c(1,1), c("Citrate_c", "H+_c"), c(1,1), "Citrate transporter")

save(joint.stoi, file = "drosophila_stoi.R")

oldstoi <- joint.stoi
match.mat <- matrix(data = 0, ncol = length(oldstoi[1,]), nrow = length(joint.stoi[1,]))
colnames(match.mat) <- colnames(oldstoi)

for(i in 1:length(joint.stoi[1,])){
for(j in 1:length(oldstoi[1,])){	
	
match.mat[i,j] <- ifelse(sum((names(joint.stoi[,i][joint.stoi[,i] != 0]) %in% names(oldstoi[,j][oldstoi[,j] != 0]))) == length(joint.stoi[,i][joint.stoi[,i] != 0]), 1, 0)
	
	
	}}

c(1:length(match.mat[1,]))[apply(match.mat, 2, sum) == 0]
oldstoi[,69]


add.reaction <- function(stoi, reactants, moles_reactants, products, moles_products, reaction_name){
	new.terms <- c(reactants[!(reactants %in% rownames(stoi))], products[!(products %in% rownames(stoi))])
	if(length(new.terms) != 0){
		stoi <- rbind(stoi, matrix(data = 0, ncol = length(stoi[1,]), nrow = length(new.terms)))
		rownames(stoi)[length(stoi[,1])-length(new.terms)+c(1:length(new.terms))] <- new.terms
		print(paste("adding", new.terms))
		}
		
		rxn.vec <- rep(0, length(stoi[,1]))
		for(i in 1:length(reactants)){
			rxn.vec[rownames(stoi) %in% reactants[i]] <- -1*moles_reactants[i]
			}
		for(i in 1:length(products)){
			rxn.vec[rownames(stoi) %in% products[i]] <- moles_products[i]
			}
		output <- cbind(stoi, rxn.vec)
		colnames(output)[length(stoi[1,])+1] <- reaction_name
		output
		}

rxn_search = function(stoiMat, search_string, is_rxn = TRUE){
	#search by metabolite or reactant and return all reactions and nonzero metabolites.
	if (is_rxn == TRUE){
		colz = grep(search_string, colnames(stoiMat), ignore.case = TRUE)
		} else {
		met = grep(search_string, rownames(stoiMat), ignore.case = TRUE)
		if (length(met) == 1){
			colz = c(1:length(stoiMat[1,]))[stoiMat[met,] != 0]
			} else {
			colz = c(1:length(stoiMat[1,]))[apply(stoiMat[met,], 2, is.not.zero)]
		}}
	
	if(length(colz) == 0){
		print("no hits")
		} else {
			rxns = stoiMat[,colz]
			if(is.vector(rxns)){
				c(colz, rxns[rxns != 0])
				} else {
					output <- rbind(colz, rxns[apply(rxns, 1, is.not.zero),])
					colnames(output) = colnames(stoiMat)[colz]
					output
					}
		}
	}

is.not.zero = function(vec){
	length(vec[vec!=0]) != 0
	}
	

write.table(write.equations(joint.stoi, compartment.mets, irreversible), file = "meta.rxns.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#x <- write.equations(joint.stoi, compartment.mets, irreversible)
#rownames(x) <- NULL

write.equations <- function(joint.stoi, compartment.mets, irreversible){
#stoi.mat = n x m stoichiometry matrix, where n is the number of metabolites and m is the number of reactions
#compartment.mets = a list where names are the compartments and nested vectors are the metabolites in that compartment
#reversible = a list of irreversible reactions from mel_FBA

reactions = colnames(joint.stoi)
rxn.formula <- rep(NA, times = length(reactions))

equation.frame <- NULL

for (i in 1:length(reactions)){
	species = data.frame(class = NA, names = names(joint.stoi[,i][joint.stoi[,i] != 0]), stoi = joint.stoi[,i][joint.stoi[,i] != 0], comp = NA)
	species$names <- as.character(species$names)
	species$class <- ifelse(species$stoi > 0, "product", "reactant")
	species$stoi <- abs(species$stoi)
	
	#determine which compartment metabolites are in and therefore whether the rxn is compartmentalized or transport
	for(el in 1:length(species$names)){
		for(comp in 1:length(names(compartment.mets))){
		if(species$names[el] %in% compartment.mets[[names(compartment.mets)[comp]]]){species$comp[el] <- names(compartment.mets)[comp]
		}}}
	
	if(reactions[i] %in% irreversible){
	irr <- 1
	} else {irr <- 0}
	
	equation.frame <- rbind(equation.frame, data.frame(class = "reaction", names = reactions[i], stoi = irr, comp = NA), species)
	}
	rbind(equation.frame, data.frame(class = "file_end", names = NA, stoi = NA, comp = NA))
	}
	
	

		
	