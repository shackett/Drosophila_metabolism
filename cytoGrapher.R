#visualizing flux results using Rcytoscape bridge to cytoscape
#creating a graphical model 

library("RCytoscape")

arm_lengths <- 5
arm_ratio <- 1
spread_angle <- 60/360*2*pi

stoisub <- S[apply(S[,1:4] != 0, 1, sum) != 0,1:4]

cofactors <- rownames(stoisub)[-c(1,3,8:11)]
core_mets <- rownames(stoisub)[!(rownames(stoisub) %in% cofactors)]

graph_center <- c(0,0)
met_pos <- data.frame(x = rep(NA, times = length(core_mets)), y = rep(NA, times = length(core_mets))); rownames(met_pos) <- core_mets
rxn_added <- rep(FALSE, times = length(stoisub[1,]))
rxn_nodes <- matrix(NA, ncol = 4, nrow = length(stoisub[,1])); colnames(rxn_nodes) <- c("rn_x", "rn_y", "pn_x", "pn_y")
cof_nodes <- NULL
#cof_nodes <- data.frame(cofactor = NA, xpos = X, ypos = X, stoi = X, stringAsFactors = FALSE)

met_pos[3,] <- c(5,5)


#within each iteration determine the reactions that haven't already been layed out and are connected to at least one already defined specie

rxn_to_do <- c(1:length(stoisub[1,]))[rxn_added == FALSE & apply(matrix(stoisub[rownames(stoisub) %in% rownames(met_pos)[apply(is.na(met_pos), 1, sum) == 0],], ncol = length(stoisub[1,]))!= 0, 2, sum) != 0]

#loop through reactions that are going to be defined and 

for(rx in 1:length(rxn_to_do)){
	
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
	
	if(length(ndefined_react) == 0 | length(ndefined_prod) == 0){
	if(length(ndefined_react) != 0){
		center_pos <- apply(met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),][ndefined_react,], 2, mean)
		
		}else{
			center_pos <- apply(met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),][ndefined_prod,], 2, mean)
			
			angle_det <- center_pos - graph_center; angle_det <- atan(angle_det[1]/angle_det[2])
			
			midpoint <- c(center_pos[1] + cos(angle_det)*(arm_lengths*3 + 2*arm_lengths*arm_ratio)/2, center_pos[2] + sin(angle_det)*(arm_lengths*3 + 2*arm_lengths*arm_ratio)/2)
			
			#define the center of mass for the products
			if(odd_prod == TRUE){
				anglez <- c(cos(atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1])), sin(atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1])))
				rxn_nodes[rxn, c(3:4)] <- center_pos + arm_lengths*arm_ratio*anglez
				#define met position
				}else{
					anglez <- atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1]) + spread_angle*arm_lengths*arm_ratio/(arm_lengths*1.5 + arm_lengths*arm_ratio); anglez <- c(cos(anglez), sin(anglez))
					rxn_nodes[rxn,c(3:4)] <- center_pos + arm_lengths*arm_ratio*anglez
					#define met position
					}
			rxn_nodes[rxn, c(1:2)] <- (arm_lengths*3 + arm_lengths*arm_ratio)*c(cos(atan((rxn_nodes[rxn,c(3:4)] - graph_center)[2]/(rxn_nodes[rxn,c(3:4)] - graph_center)[1])), sin(atan((rxn_nodes[rxn,c(3:4)] - graph_center)[2]/(rxn_nodes[rxn,c(3:4)] - graph_center)[1])))
			
			anglez <- atan(c(rxn_nodes[rxn, c(1:2)] - rxn_nodes[rxn, c(3:4)])[2]/c(rxn_nodes[rxn, c(1:2)] - rxn_nodes[rxn, c(3:4)])[1])
			splayed_angle <- ifelse(odd_react, angle_set_odd[c(1:n_react)], angle_set_even[c(1:n_react)])
			
			c(cos(anglez + splayed_angle), sin(anglez + splayed_angle))
			
			}
		angle_set_odd <- c(0, spread_angle, -spread_angle, 2*spread_angle, -2*spread_angle)
		angle_set_even <- c(spread_angle/2, -spread_angle/2, spread_angle*3/2, -spread_angle*3/2, spread_angle*5/2, -spread_angle*5/2)
	}
	}
	
	principal_change[principal_change < 1]
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	