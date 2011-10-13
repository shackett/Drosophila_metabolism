setwd("/Users/seanhackett/Desktop/Cornell/mel_metabolism/")


#functions to find reactions in the yeast metabolomic reconstruction
load("yeast_stoi.R")
rxns = rxn_search(stoiMat, "glycerol kinase", is_rxn = TRUE)
cols <- c(232, 596, 816)
redmat <- stoiMat[apply(stoiMat[,cols], 1, is.not.zero), cols]


valid.reactions <- read.table("/Users/seanhackett/Desktop/Cornell/mel_metabolism/mel_network.csv", sep = ",", blank.lines.skip = FALSE, header = TRUE)

cytosolic.reactions <- valid.reactions[valid.reactions$Compartment == "Cytosol",][!is.na(valid.reactions[valid.reactions$Compartment == "Cytosol",]$ReactionNumber),]
cytosolic.stoi <- stoiMat[,cytosolic.reactions$ReactionNumber][apply(stoiMat[,cytosolic.reactions$ReactionNumber], 1, is.not.zero),]
rownames(cytosolic.stoi) <- lapply(rownames(cytosolic.stoi), paste, "_c", sep = "")

mitochondrial.reactions <- valid.reactions[valid.reactions$Compartment == "Mitochondria",][!is.na(valid.reactions[valid.reactions$Compartment == "Mitochondria",]$ReactionNumber),]
mitochondrial.stoi <- stoiMat[, mitochondrial.reactions $ReactionNumber][apply(stoiMat[, mitochondrial.reactions $ReactionNumber], 1, is.not.zero),]
rownames(mitochondrial.stoi) <- lapply(rownames(mitochondrial.stoi), paste, "_m", sep = "")

#Hand-added reactions

added.reactions <- valid.reactions[is.na(valid.reactions$ReactionNumber),]

cytosolic.stoi <- add.reaction(cytosolic.stoi, "Trehalose_c", 1, NULL, NULL, "Trehalose usage")
cytosolic.stoi <- add.reaction(cytosolic.stoi, "glycogen_c", 1, NULL, NULL, "Glycogen usage")
cytosolic.stoi <- add.reaction(cytosolic.stoi, "Hexadecanoate (n-C16:0)_c", 1, NULL, NULL, "Palmitate usage")
cytosolic.stoi <- add.reaction(cytosolic.stoi, "alpha-D-Ribose 5-phosphate_c", 1, NULL, NULL, "Ribose biosynthesis")
cytosolic.stoi <- add.reaction(cytosolic.stoi, "CO2_c", 1, NULL, NULL , "CO2 leaving")
cytosolic.stoi <- add.reaction(cytosolic.stoi, NULL, NULL, "O2_c", 1, "O2 entering")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("Diphosphate_c", "H2O_c"), c(1,1), "Phosphate_c", 2, "Pyrophosphatase")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("H2O_c", "Coenzyme A_c", "ATP_c", "Citrate_c") , c(1,1,1,1), 
c("Oxaloacetate_c", "Acetyl-CoA_c", "ADP_c", "Phosphate_c"), c(1,1,1,1), "Citrate lyase")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("Carnitine_c", "Palmitoyl-CoA (n-C16:0CoA)_c"), c(1,1), c("Palmitoylcarnitine_c", "Coenzyme A_c"), c(1,1), "Palmitoylcarnitine synthesis")
cytosolic.stoi <- add.reaction(cytosolic.stoi, c("UDP_c", "ATP_c"), c(1,1), c("UTP_c", "ADP_c"), c(1,1), "Nucloside diphosphate kinase (UTP-ATP)")

mitochondrial.stoi <- add.reaction(mitochondrial.stoi, c("Palmitoylcarnitine_m", "Coenzyme A_m"), c(1,1), c("Carnitine_m", "Palmitoyl-CoA (n-C16:0CoA)_m"), c(1,1), "Palmitoylcarnitine breakdown")
mitochondrial.stoi <- add.reaction(mitochondrial.stoi, 
c("Palmitoyl-CoA (n-C16:0CoA)_m", "Flavin adenine dinucleotide oxidized_m", "H2O_m", "Nicotinamide adenine dinucleotide_m", "Coenzyme A_m"),
c(1, 7, 7, 7, 7), c("Flavin adenine dinucleotide reduced_m", "Nicotinamide adenine dinucleotide - reduced_m", "H+_m", "Acetyl-CoA_m"), c(7, 7, 7, 8),
"Palmitoyl-CoA catabolism")
mitochondrial.stoi <- add.reaction(mitochondrial.stoi, c("Succinyl-CoA_m", "Phosphate_m", "GDP_m"), c(1,1,1), c("Succinate_m", "Coenzyme A_m", "GTP_m"), c(1,1,1), "Succinyl CoA synthetase")


joint.stoi <- cbind(rbind(cytosolic.stoi, matrix(data = 0, ncol = length(cytosolic.stoi[1,]), nrow = length(mitochondrial.stoi[,1]))), rbind(matrix(data = 0, ncol = length(mitochondrial.stoi[1,]), nrow = length(cytosolic.stoi[,1])), mitochondrial.stoi))
rownames(joint.stoi) <- c(rownames(cytosolic.stoi), rownames(mitochondrial.stoi))
colnames(joint.stoi) <- c(colnames(cytosolic.stoi), colnames(mitochondrial.stoi))

joint.stoi <- add.reaction(joint.stoi, c("Glycerol 3-phosphate_c", "Flavin adenine dinucleotide oxidized_m"), c(1,1), 
c("Dihydroxyacetone phosphate_c", "Flavin adenine dinucleotide reduced_m"), c(1,1), "Glycerol 3-phosphate shuttle")
joint.stoi <- add.reaction(joint.stoi, c("Flavin adenine dinucleotide reduced_m", "H+_m"), c(1, 6), c("Flavin adenine dinucleotide oxidized_m", "H+_c"), c(1, 6), "FADH2 oxidation for proton transport")
joint.stoi <- add.reaction(joint.stoi, c("Nicotinamide adenine dinucleotide phosphate - reduced_m", "H+_m"), c(1, 10), c("Nicotinamide adenine dinucleotide phosphate_m", "H+_c"), c(1, 10), "NADH oxidation for proton transport")
joint.stoi <- add.reaction(joint.stoi, c("ADP_m", "Phosphate_m", "H+_c"), c(1, 1, 3), c("ATP_m", "H+_m"), c(1, 3), "ATP synthetase")
joint.stoi <- add.reaction(joint.stoi, c("GTP_m", "ADP_m"), c(1,1), c("GDP_m", "ATP_m"), c(1,1), "Nucleoside diphosphate kinase (GTP-ATP)")
joint.stoi <- add.reaction(joint.stoi, c("Pyruvate_c", "H+_c"), c(1,1), c("Pyruvate_m", "H+_m"), c(1,1), "Pyruvate transporter")
joint.stoi <- add.reaction(joint.stoi, c("Phosphate_c", "H+_c"), c(1,1), c("Phosphate_m", "H+_m"), c(1,1), "Phosphate transporter")
joint.stoi <- add.reaction(joint.stoi, "ADP_c", 1, "ATP_c", 1, "ATP/ADP antiport")
joint.stoi <- add.reaction(joint.stoi, c("L-Malate_m", "Phosphate_c"), c(1,1), c("L-Malate_c", "Phosphate_m"), c(1,1), "Dicarboxylate carrier")
joint.stoi <- add.reaction(joint.stoi, c("Citrate_m", "H+_m", "L-Malate_c"), c(1,1,1), c("Citrate_c", "H+_c", "L-Malate_m"), c(1,1,1), "Tricarboxylate carrier")

#joint.stoi <- add.reaction(joint.stoi, c("L-Malate_m", "H+_m"), c(1,1), c("L-Malate_c", "H+_c"), c(1,1), "Malate transporter")
#joint.stoi <- add.reaction(joint.stoi, c("Citrate_m", "H+_m"), c(1,1), c("Citrate_c", "H+_c"), c(1,1), "Citrate transporter")





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
	
	
	








#H20 - exchange & mito
#CO2/O2 mitochondrial exchange
#UTP/ATP trans-phosphorylation



transport
921 - trehalose
455 - glucose transport
glycogen transport

trehalose
918 - trehalose p
919 - alpha trehalose p synthase
920 - alpha trehalse
518 - HEX

Glycogen
439 - branching
454 - GP
477 - GS 
434 - UDP-glucose synthesis
? - pyrophosphatase

glycolysis
735 - PGM
731 - PGI
390 - aldolase
726 - PFK
916 - TPI
437 - glyceraldehyde3P: G3PD
732 - phosphoglycerate K
734 - phosphoglycerate M
322 - enolase
810 - pyruvate kinase

PPP
4 - G6PD
733 - 6PGlactonase
483 - 6PGD
835 - ribulose 5P epimerase
836 - ribose 5P isomerase
886 - transaldolase
906, 907 - transketolase

lower glycolytic branches

606, 607 - NAD/NADP malic enzyme
760 - pep carboxykinase
604 - malate dehydrogenase

420 - glycerol3PD - cytoNAD
421 - glycerol3PD - mitoFAD
120 - ADH
136 - aldehyde dehydrogenase
82 - acetylCoAsynthase

### Mito transport rxns

232 - citrate out
596 - malate out
816 - pyruvate in
GPO shuttle
palmiotyle CoA

### FAS

56 - acetylCoA carboxylase: acetylCoA -> malonylCoA
75 - acetylCoA transacylase
603 - malonylCoA transacylase

388 - C8 synthesis lumped
365 - c10 - Decanoyl CoA
368 - c12 - Dodecanoyl CoA
371 - c14 - Tetradecanoyl CoA
376 - c16 - palmitoyl CoA synth

387, 364, 367, 370, 375 - ACP reactions
334 - palmitoylACP hydrolysis -> Hexadecanoate
789 - acyl-CoA thioesterase -> Palmitoyl-CoA
237 - carnitine transport

palmitoylCoA +  carnitine -> palmitoyl-carnitine
palmitoylCoA transport
palmitoylcarnitine + CoA -> carnitine + palmitoylCoA

###oxidation

### Mitochondrial rxns

720 - acetylCoa production by PDH
714 - pyruvate carboxylase

TCA

242 - citrate synthase
77 - aconitase
549 - isocitrate dehydrogenase, 550 - NADP IDH
111, 112, 443 - alpha-ketoglutarate dehydrogenase
881 - succinate-CoA ligase
878 - SDH
415 - FUM
604 - malate dehydrogenase


Misc

82 - acetyl coa synthetase
180 -  ATP adenylyltransferase
183,184,186
742 - phosphate transporter
511 - bicarbonate equilibrium


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
	



		
	