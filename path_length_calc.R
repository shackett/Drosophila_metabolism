# Assay plate geometry:

top_d = 0.0069
bottom_d = 0.0064
height = 0.0107
height_scale = ((top_d - bottom_d)/height)

radius_x <- function(x, top_d, bottom_d, height){
	((x * (top_d - bottom_d)/height) + bottom_d)/2
	}

radius_x(height, top_d, bottom_d, height)

pl = 0.0107

vol = 0.002

optimize(find.pl, c(0, 20), volume = vol, height_scale = height_scale, bottom_d = bottom_d, tol = 10^-10)

find.pl <- function(pl, volume, height_scale, bottom_d){

abs((pi*((1/12)*(height_scale^2)*(pl^3) + (1/4)*height_scale*bottom_d*(pl^2) + (1/4)*(bottom_d^2)*pl)*1000) - vol)
}




(pi*(radius_x(pl, top_d, bottom_d, height)^2 + radius_x(0, top_d, bottom_d, height)^2)/2)*pl*1000
