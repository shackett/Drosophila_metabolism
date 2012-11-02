# Assay plate geometry:

############# geometry of a UV-clear plate ############

top_d = 0.00686
bottom_d = 0.00635
height = 0.01067
height_scale = ((top_d - bottom_d)/(2*height))
avg_d <- mean(top_d, bottom_d)



#calculate path length in costar plates

find.pl <- function(pl, volume, height_scale, bottom_d){

#if cylinder: abs(((avg_d/2)^2*pi*pl*10^6) - volume)

#determine the volume (in mL) as a function of path-length (in meters) for costar plates
abs((pi*(bottom_d^2/4*pl + bottom_d*height_scale*pl^2/2 + height_scale^2*pl^3/3)*10^6) - volume)

}




find.pl <- function(pl, volume, height_scale, bottom_d){
#calculate path length in 
abs((pi*((1/12)*(height_scale^2)*(pl^3) + (1/4)*height_scale*bottom_d*(pl^2) + (1/4)*(bottom_d^2)*pl)*1000) - volume)
}

pl_optim <- function(volume){optimize(find.pl, c(0, 20), volume = volume, height_scale = height_scale, bottom_d = bottom_d, tol = 10^-10)}
	

########### geometry of a flimsy - visible - plate ############

data <- read.delim("flimsy_plate_pl.txt", skip= 9, sep = "\t", header = FALSE)
data <- data[,-c(1, 2, 15)]

pl.table <- data.frame(volume = c(50, 100, 125, 150, 175, 200), pl = rep(NA, times = 6), sd.pl = rep(NA, times = 6))
pl.table$volume <- pl.table$volume / 1000

for(i in 1:6){
	pl.table$pl[i] <- mean(unlist(data[,c((i*2)-1, (i*2))]))
	pl.table$sd.pl[i] <- sd(unlist(data[,c((i*2)-1, (i*2))]))
	}

vol.fit <- summary(lm(pl.table$pl ~ pl.table$volume))$coef[,1]


flimsy_pl <- function(volume){
#returns a path-length in cm when provided a volume in mL for visible assays
(vol.fit[1] + vol.fit[2]*volume)
}

### use ####
#pathlengths in cm when provided mL
pl_optim(0.2)$minimum * 100
flimsy_pl(0.2)