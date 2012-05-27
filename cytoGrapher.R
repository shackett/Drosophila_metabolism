#visualizing flux results using Rcytoscape bridge to cytoscape
#creating a graphical model 

library("RCytoscape")
library(gplots)

setwd("/Users/seanhackett/Desktop/Cornell/Drosophila_metabolism/")
load("drosophila_stoi.R")

arm_lengths <- 4
arm_ratio <- 1
spread_angle <- 60/360*2*pi
angle_set_odd <- c(0, spread_angle, -spread_angle, 2*spread_angle, -2*spread_angle)
angle_set_even <- c(spread_angle/2, -spread_angle/2, spread_angle*3/2, -spread_angle*3/2, spread_angle*5/2, -spread_angle*5/2)
		
stoisub <- joint.stoi[apply(joint.stoi[,1:10] != 0, 1, sum) != 0,1:10]

cofactors <- c("H+_c", "H2O_c", "Phosphate_c", "Diphosphate_c", "ATP_c", "ADP_c", "UTP_c", "UDP_c")
core_mets <- rownames(stoisub)[!(rownames(stoisub) %in% cofactors)]

graph_center <- c(0,0)
met_pos <- data.frame(x = rep(NA, times = length(core_mets)), y = rep(NA, times = length(core_mets))); rownames(met_pos) <- core_mets
rxn_added <- rep(FALSE, times = length(stoisub[1,]))
rxn_nodes <- matrix(NA, ncol = 4, nrow = length(stoisub[1,])); colnames(rxn_nodes) <- c("rn_x", "rn_y", "pn_x", "pn_y")
cof_nodes <- NULL
#cof_nodes <- data.frame(cofactor = NA, xpos = X, ypos = X, stoi = X, stringAsFactors = FALSE)

met_pos[2,] <- c(1,sqrt(3))


#within each iteration determine the reactions that haven't already been layed out and are connected to at least one already defined specie

rxn_to_do <- c(1:length(stoisub[1,]))[rxn_added == FALSE & apply(matrix(stoisub[rownames(stoisub) %in% rownames(met_pos)[apply(is.na(met_pos), 1, sum) == 0],], ncol = length(stoisub[1,]))!= 0, 2, sum) != 0]

rxn_added[rxn_to_do] <- TRUE

#loop through reactions that are going to be defined and 

for(rx in rxn_to_do){
	
	rxn_stoi <- stoisub[,rx][stoisub[,rx] != 0]
	cofactor_change <- rxn_stoi[names(rxn_stoi) %in% cofactors]
	principal_change <-  rxn_stoi[names(rxn_stoi) %in% core_mets]
	 
	odd_react <- odd(length(principal_change[principal_change < 0]))
	odd_prod <- odd(length(principal_change[principal_change > 0]))
	n_react <- length(principal_change[principal_change < 0])
	n_prod <- length(principal_change[principal_change > 0])
	ndefined_react <- c(1:n_react)[apply(!is.na(met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),]), 1, sum) != 0]
	ndefined_prod <- c(1:n_prod)[apply(!is.na(met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),]), 1, sum) != 0]
	
	#if only either a subset of products or reactants is defined, but not both, the principal direction vector (going from reactants to products) is determined by the center of the graph.  Otherwise this vector is determined by the position of the defined metabolites, making adjustments to account for whether the number of principal products and reactants is odd or even
	#if there are already reactions attached to a metabolite polarize the new reaction in the opposite direction from the mean angle
	
	if(length(ndefined_react) == 0 | length(ndefined_prod) == 0){
		#use ifelse to indicate whether reactants or products were defined
		
		#redo ifelse statements to allow for more values returned
		
		reactDef = ifelse(length(ndefined_react) != 0, TRUE, FALSE)
		
		#center_pos <- apply(met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),][ndefined_react,], 2, mean)
		
		if(reactDef){changing <- principal_change < 0}else{changing <- principal_change > 0}
		center_pos <- apply(met_pos[rownames(met_pos) %in% names(principal_change[changing]),][ifelse(reactDef, ndefined_react, ndefined_prod),], 2, mean)
		
		#test_exist <- matrix(rxn_nodes[stoisub[rownames(stoisub) %in% rownames(met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),][ndefined_react,]),] != 0,], ncol = 4, byrow = FALSE)
		
		test_exist <- matrix(rxn_nodes[stoisub[rownames(stoisub) %in% rownames(met_pos[rownames(met_pos) %in% names(principal_change[changing]),][ifelse(reactDef, ndefined_react, ndefined_prod),]),] != 0,], ncol = 4, byrow = FALSE)
		
		if(sum(!is.na(test_exist[,1])) != 0){
			graph_center <- apply(matrix(test_exist[c(1:length(test_exist[,1]))[!is.na(test_exist[,1])],], ncol = 4, byrow = FALSE), 2, mean)[1:2]
			}else{
				graph_center <- c(0,0)
				}
		
		angle_det <- center_pos - graph_center; angle_det <- ifelse((center_pos - graph_center)[1] >= 0, atan(angle_det[2]/angle_det[1]), pi + atan(angle_det[2]/angle_det[1]))
			
		midpoint <- c(center_pos[1] + cos(angle_det)*(arm_lengths*3 + 2*arm_lengths*arm_ratio)/2, center_pos[2] + sin(angle_det)*(arm_lengths*3 + 2*arm_lengths*arm_ratio)/2)
			
		#define the center of mass for the side of the reaction with defined species
		if((odd_react == TRUE & reactDef) | (odd_prod == TRUE & !reactDef)){
			anglez <- ifelse((midpoint - center_pos)[1] >= 0, atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1]), pi + atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1])); anglez <- c(cos(anglez), sin(anglez))
			
			if(reactDef){cell_choice <- c(1:2)}else{cell_choice <- c(3:4)}
			rxn_nodes[rx, cell_choice] <- center_pos + arm_lengths*arm_ratio*anglez
			#define met position
			}else{
				anglez <- ifelse((midpoint - center_pos)[1] >= 0, atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1]), pi + atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1])) + spread_angle*arm_lengths*arm_ratio/(arm_lengths*1.5 + arm_lengths*arm_ratio); anglez <- c(cos(anglez), sin(anglez))
				
				if(reactDef){cell_choice <- c(1:2)}else{cell_choice <- c(3:4)}
				rxn_nodes[rx, cell_choice] <- center_pos + arm_lengths*arm_ratio*anglez
				#define met position
					}
					
			#define the center of mass for the previously undefined side of the reaction		
			anglez <- ifelse(c(rxn_nodes[rx,cell_choice] - graph_center)[1] >= 0, atan((rxn_nodes[rx,cell_choice] - graph_center)[2]/(rxn_nodes[rx,cell_choice] - graph_center)[1]), pi + atan((rxn_nodes[rx,cell_choice] - graph_center)[2]/(rxn_nodes[rx,cell_choice] - graph_center)[1])); anglez <- c(cos(anglez), sin(anglez))
			
			cell_choice <- list()
			if(reactDef){cell_choice[[1]] <- c(3:4); cell_choice[[2]] <- c(1:2)}else{cell_choice[[1]] <- c(1:2); cell_choice[[2]] <- c(3:4)}
			rxn_nodes[rx, cell_choice[[1]]] <- rxn_nodes[rx, cell_choice[[2]]] + (arm_lengths*3 + arm_lengths*arm_ratio)*anglez
			
			anglez <- ifelse((rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[1] >= 0, atan((rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[2]/(rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[1]), pi + atan((rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[2]/(rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[1]))
			
			if(odd_react == TRUE){splayed_angle = angle_set_odd[c(1:n_react)]}else{splayed_angle = angle_set_even[c(1:n_react)]}
			
			met_posn <- lapply(anglez + splayed_angle, function(x){rxn_nodes[rx, cell_choice[[1]]] + c(cos(x), sin(x))*arm_lengths*arm_ratio})
			if(reactDef){changing <- principal_change > 0}else{changing <- principal_change < 0}
		
			for(i in 1:length(met_posn)){
				met_pos[rownames(met_pos) == names(principal_change[changing]),][i,] <- met_posn[[i]]
				}
				}
				
				if(length(ndefined_react) != 0 & length(ndefined_prod) != 0){
				
					#if there are defined reactants and products then the direction of the reaction edge linking them will be determined by their position with some offset to account for whether their is an odd or even number of reactants/products
					print("woot")
					}
	
	
	
	}
	
	
	
#add in exceptions:
#some cofactors should function like primary metabolites for a subset of reactions. e.g. protons should be cofactors for most reactions but are of primary interest for ATP-synthase.  
#some primary metabolites should be split into subsets whose position is indexed with a list	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	