setwd("/Users/seanhackett/Desktop/Cornell/Drosophila_metabolism/")
data <- read.delim("flimsy_plate_pl.txt", skip= 9, sep = "\t", header = FALSE)
data <- data[,-c(1, 2, 15)]

pl.table <- data.frame(volume = c(50, 100, 125, 150, 175, 200), pl = rep(NA, times = 6), sd.pl = rep(NA, times = 6))


for(i in 1:6){
	pl.table$pl[i] <- mean(unlist(data[,c((i*2)-1, (i*2))]))
	pl.table$sd.pl[i] <- sd(unlist(data[,c((i*2)-1, (i*2))]))
	}

vol.fit <- summary(lm(pl.table$pl ~ pl.table$volume))$coef[,1]

#plot(pl.table$pl ~ pl.table$volume, xlab = "path-length (cm)")
#lines((vol.fit[1] + vol.fit[2]*pl.table$volume) ~ pl.table$volume)

flimsy_pl <- function(volume){
#returns a path-length in cm when provided a volume in mL
(vol.fit[1] + vol.fit[2]*volume)
}