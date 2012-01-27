#timing of respirometry (major chunk - 1 of 3)

#input - timing: line, day, incubation start/stop for n samples
#output - time: line, day, time in minutes & fractional seconds of incubation

#time directory "C:\\Users\\Sean\\Desktop\\Respirometry\\currenttiming8.01.10csv.csv")
timing <- read.table("C:\\Users\\Sean\\Desktop\\Respirometry\\currenttiming8.01.10csv.csv", header = TRUE, sep = ",")

time <- data.frame(line = timing[,1], day = timing[,2], R1 = NA, R2 = NA, R3 = NA)
timest <- data.frame(line = timing[,1], day = timing[,2], R1 = NA, R2 = NA, R3 = NA)
nlines <- length(timing[,1])

timedif <- data.frame(R1 = c(NA,NA), R2 = c(NA,NA), R3 = c(NA,NA))

#determine time of respiration
resptime <- function(to, tf) {
if (is.na(to) == TRUE){
to <- "0.0.0"
tf <- "0.0.0"
}
timedif[1,] <- unlist(strsplit(to, "\\."))
timedif[2,] <- unlist(strsplit(tf, "\\."))
timedif11 <- as.integer(timedif[1,1])
timedif12 <- as.integer(timedif[1,2])
timedif13 <- as.integer(timedif[1,3])
timedif21 <- as.integer(timedif[2,1])
timedif22 <- as.integer(timedif[2,2])
timedif23 <- as.integer(timedif[2,3])
timedifh <- timedif21 - timedif11
timedifm <- timedif22 - timedif12
timedifs <- timedif23 - timedif13
if(timedifs < 0){
timedifm <- timedifm - 1 
timedifs <- timedifs + 60
	}
if(timedifm < 0){
timedifh <- timedifh - 1
timedifm <- timedifm + 60
	}
minutesresp <- timedifh*60 + timedifm + timedifs * (1/60)
	}

#function that outputs the start of incubation (to) with minutes and seconds as fractional hours
timestart <- function(to){
hms <- unlist(strsplit(to, "\\."))
as.integer(hms[1]) + as.integer(hms[2])/60 + as.integer(hms[3])/60^2
} 

for (i in 1:nlines){
for (j in 1:3){

if (j == 1){
to <- as.character(timing[i,3])
tf <- as.character(timing[i,4])
time[i,3] <- resptime(to,tf)
timest[i,3] <- timestart(to)
	}
if (j == 2){
to <- as.character(timing[i,5])
tf <- as.character(timing[i,6])
time[i,4] <- resptime(to,tf)
timest[i,4] <- timestart(to)	
	}
if (j == 3){
to <- as.character(timing[i,7])
tf <- as.character(timing[i,8])
time[i,5] <- resptime(to,tf)
timest[i,5] <- timestart(to)	
	}

}
}





# peak area calculation (major chunk - 2 of 3)

#function of respirometry file {partial & barometric pressure for co2 and ox, temp, flow rate} and time file {vector of incubation times - 1-3 times}
#supplied data - pcot(by file), incubtime
#derived variables - nmark (from #incubtime)
#supplied variables - co2 & ox truncating buffer, "a" transition matrix, cotblank effect
#extracted data - co2 produced per minute, o2 consumed per minute (in units = ) for n = nmarks

#file <- read.csv("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\3_31\\3_31B05_B10_data")
#file <- read.csv("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\3_31\\3_31B10_B11_data")
#file <- read.csv("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\3_13\\3_13N10_data")
#cotblank (1 of 2) file <- read.csv("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\4_1\\9timedBL_data")
#cotblank (2 of 2) file <- read.csv("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\4_2\\9blanks2_data")

#max number of peaks
nmax <- 3

#externalized variables for easy modification
#cropping parameters (on either side of peaks)
buffercext <- c(10,105)
bufferoext <- 25

#CO2 HMM transition matrix probabilities 1->1, 2->1, 1->2, 2->2
# 0.99, 0.05, 0.01, 0.5
acot <- matrix(data = c(0.99,0.5,0.01,0.5), ncol=2, nrow=2)

#O2 HMM transition matrix probabilities 1->1, 2->1, 1->2, 2->2
# 0.99, 0.05, 0.01, 0.5
aox <- matrix(data = c(0.99,0.15,0.01,0.8), ncol=2, nrow=2)

#cotblank - the mean effect of injecting a fly-less syringe, on the co2 curve (as determined by 18 replicates (2 blocks)). 
#cotareablank <- 0.0621
cotvolblk <- 6.730346e-07
#cotmolesblank <- 2.413348 * 10^-10

peakarea <- NA
peakarea <- function(file, nmark){

#%co2 (file[,3]/file[,6]) partial / barometric pressure, barometric p, flow rate
pcot <- cbind(file[,1]+1,(file[,3]/file[,6]), file[,6], file[,2])

#flow rate (mL/min -> L/s)
flowr <- mean(file[,2])*(1/60)*(1/1000)

#npad - pad output with NAs
npad <- nmax - nmark

#define truncation buffer for co2
if(k %in% cbuffermod[,1]) {bufferc <- cbuffermod[,2:3][cbuffermod[,1] == k]} else {bufferc <- buffercext}

coef <- summary(lm(pcot[,2] ~ pcot[,1]))$coef[,1]
fit <- cbind(pcot[,1], (pcot[,1]*coef[2] + coef[1]))

pcot2 <- pcot[(pcot[,2] < fit[,2]),]

border <- rep(NA, times = (length(pcot2[,1])-1))
for (i in 1:length(border)){
	border[i] <- pcot2[i+1,1] - pcot2[i,1]  
	}
lgaps <- sort(border, decreasing=T)[1:nmark]
lgapstart <- (1:length(border))[border %in% lgaps]

remove <- NULL

for (i in 1:nmark){
	remove <- c(remove, (lgapstart[i] - bufferc[1]):(lgapstart[i] + (bufferc[2] + 1))) 
	}
pcot3 <- pcot2[(!(pcot2[,1] %in% pcot2[,1][remove])),]
coef3 <- summary(lm(pcot3[,2] ~ pcot3[,1] + I(pcot3[,1]^2) + I(pcot3[,1]^3)))$coef[,1]
fit3 <- cbind(pcot3[,1], (pcot3[,1]^3*coef3[4] + pcot3[,1]^2*coef3[3] + pcot3[,1]*coef3[2] + coef3[1]))
overallfit <- cbind(pcot[,1], (pcot[,1]^3*coef3[4] + pcot[,1]^2*coef3[3] + pcot[,1]*coef3[2] + coef3[1]))

sdfit <- sqrt(sum(lm(pcot3[,2] ~ pcot3[,1] + I(pcot3[,1]^2) + I(pcot3[,1]^3))$residuals^2)/(lm(pcot3[,2] ~ pcot3[,1] + I(pcot3[,1]^2) + I(pcot3[,1]^3))$df.res))

# HMM - prob 2
t <- pcot[,1]
t1 <- t[1]
deltat <- matrix(data = NA, ncol = 3, nrow = length(pcot[,1]))
psit <- matrix(data = NA, ncol = 3, nrow = length(pcot[,1]))
deltat[,1] <- pcot[,1]
psit[,1] <- pcot[,1]

#initialization
i <- 1
piyum <- matrix(data= c(1,0))
bt <- matrix(data = c(pnorm(pcot[t1,2],overallfit[overallfit[,1] == t1,2], sdfit, lower.tail=FALSE), (1-pnorm(pcot[t1,2],overallfit[overallfit[,1] == t1,2], sdfit, lower.tail=FALSE))))
deltat[i,2:3] <- piyum * bt 
psit[i,2:3] <- c(0,0)

#recursion
for (i in 2:length(pcot[,1])){
	bt <- matrix(data = c(pnorm(pcot[i,2],overallfit[overallfit[,1] == t[i],2], sdfit, lower.tail=FALSE), 1- (pnorm(pcot[i,2],overallfit[overallfit[,1] == t[i],2], sdfit, lower.tail=FALSE))))
	deltat[i,2] <- max((deltat[(i-1),2] * bt[1] * acot[1,1]),(deltat[(i-1),3] * bt[2] * acot[2,1]))
	deltat[i,3] <- max((deltat[(i-1),3] * bt[2] * acot[2,2]),(deltat[(i-1),2] * bt[1] * acot[1,2]))
	psit[i,2] <- ifelse(max((deltat[(i-1),2] * acot[1,1] * bt[1]), (deltat[(i-1),3] * acot[1,2] * bt[2]))==(deltat[(i-1),2] * acot[1,1] * bt[1]), 0, 1)
	psit[i,3] <- ifelse(max((deltat[(i-1),3] * acot[2,2] * bt[2]), (deltat[(i-1),2] * acot[2,1] * bt[1]))==(deltat[(i-1),3] * acot[2,2] * bt[2]), 1, 0)
	}

#termination
QT <- matrix(data = NA, ncol = 2, nrow = length(psit[,1]))
QT[,1] <- pcot[,1]
lpcot <- length(pcot[,1]) 

PT <- max(deltat[lpcot,2:3])
QT[lpcot,2] <- ifelse(deltat[lpcot,2]== PT, psit[lpcot,2], psit[lpcot,3])

#path (state sequence) backtracking

for (i in (lpcot-1):1){
	if(QT[(i+1),2] == 0){
		QT[i,2] <- ifelse(psit[i,2] == 0, 0, 1)
			} 
	if(QT[(i+1),2] == 1){
		QT[i,2] <- ifelse(psit[i,3] == 1, 1, 0)
			}
	}

plot(pcot[,2], pch = 16, col = ifelse (QT[,2] == 1, "red", "blue"), main = paste (k, line, sep = "-"))
points(overallfit[,2], pch = 16, col = "black")

plot(pcot3[,2] ~ pcot3[,1], col = "green")
points(overallfit[,2] ~ overallfit[,1], pch=16)

lpeaks <- rle(QT[,2])$len[rle(QT[,2])$val == 1]
npeaks <- length(lpeaks)
cvol <- c(rep(NA, times= npeaks))
cvol <- c(rep(NA, times= npeaks))



for (i in 1:npeaks){
	
	peakpoints <- pcot[(sum(rle(QT[,2])$lengths[1:(2*i-1)])+1):(sum(rle(QT[,2])$lengths[1:(2*i)])),1]
	fitted <- matrix(overallfit[overallfit[,1] %in% peakpoints,], ncol =2)
	cvol[i] <- sum((pcot[peakpoints,2] - fitted[,2])*flowr)	
	}

markvalvol <- sort(cvol, decreasing = T)[1:nmark]
cvol <- c(cvol[cvol %in% markvalvol], rep(NA, times=npad)) - cotvolblk 

##################################################################################################
#%oxygen (file[,4]/file[,7]) partial / barometric pressure
pox <- cbind(file[,1]+1,(file[,4]/file[,7]), file[,7])

#define truncation buffer for o2
if(k %in% obuffermod[,1]) {buffero <- obuffermod[,2][obuffermod[,1] == k]} else {buffero <- bufferoext}

#truncating span points and defining background (pox3)

rem <- 1:pox[1:50,1][abs(diff(pox[1:50,2])) == max(abs(diff(pox[1:50,2])))]
pox <- pox[(!(pox[,1] %in% pox[,1][rem])),]
coef <- summary(lm(pox[,2] ~ pox[,1] + I(pox[,1]^2)))$coef[,1]
fit <- cbind(pox[,1], (pox[,1]^2*coef[3] + pox[,1]*coef[2] + coef[1]))

#if fit still remains poor
if (k %in% nuisances){
	coef <- summary(lm(pox[,2] ~ pox[,1] + I(pox[,1]^2) + I(pox[,1]^3)))$coef[,1]
	fit <- cbind(pox[,1], (pox[,1]^3*coef[4] + pox[,1]^2*coef[3] + pox[,1]*coef[2] + coef[1]))
	}

pox2 <- pox[(pox[,2] > fit[,2]),]

border <- rep(NA, times = (length(pox2[,1])-1))
for (i in 1:length(border)){
	border[i] <- pox2[i+1,1] - pox2[i,1]  
	}
lgaps <- sort(border, decreasing=T)[1:nmark]
lgapstart <- (1:length(border))[border %in% lgaps]

remove <- NULL

for (i in 1:nmark){
	remove <- c(remove, (lgapstart[i] - buffero):(lgapstart[i] + (buffero + 1))) 
	}
pox3 <- pox2[(!(pox2[,1] %in% pox2[,1][remove])),]
coef3 <- summary(lm(pox3[,2] ~ pox3[,1] + I(pox3[,1]^2) + I(pox3[,1]^3) + I(pox3[,1]^4)))$coef[,1]
fit3 <- cbind(pox3[,1], (pox3[,1]^4*coef3[5] + pox3[,1]^3*coef3[4] + pox3[,1]^2*coef3[3] + pox3[,1]*coef3[2] + coef3[1]))

overallfit <- cbind(pox[,1], (pox[,1]^4*coef3[5] + pox[,1]^3*coef3[4] + pox[,1]^2*coef3[3] + pox[,1]*coef3[2] + coef3[1]))

sdfit <- sqrt(sum(lm(pox3[,2] ~ pox3[,1] + I(pox3[,1]^2) + I(pox3[,1]^3) + I(pox3[,1]^4))$residuals^2)/(lm(pox3[,2] ~ pox3[,1] + I(pox3[,1]^2) + I(pox3[,1]^3) + I(pox3[,1]^4))$df.res))



# HMM - prob 2
t <- pox[,1]
t1 <- t[1]
deltat <- matrix(data = NA, ncol = 3, nrow = length(pox[,1]))
psit <- matrix(data = NA, ncol = 3, nrow = length(pox[,1]))
deltat[,1] <- pox[,1]
psit[,1] <- pox[,1]

#initialization
i <- 1
piyum <- matrix(data= c(1,0))
bt <- matrix(data = c(pnorm(pox[t1,2],overallfit[overallfit[,1] == t1,2], sdfit), (1-pnorm(pox[t1,2],overallfit[overallfit[,1] == t1,2], sdfit))))
deltat[i,2:3] <- piyum * bt 
psit[i,2:3] <- c(0,0)

#recursion
for (i in 2:length(pox[,1])){
	bt <- matrix(data = c(pnorm(pox[i,2],overallfit[overallfit[,1] == t[i],2], sdfit), 1- (pnorm(pox[i,2],overallfit[overallfit[,1] == t[i],2], sdfit))))
	deltat[i,2] <- max((deltat[(i-1),2] * bt[1] * aox[1,1]),(deltat[(i-1),3] * bt[2] * aox[2,1]))
	deltat[i,3] <- max((deltat[(i-1),3] * bt[2] * aox[2,2]),(deltat[(i-1),2] * bt[1] * aox[1,2]))
	psit[i,2] <- ifelse(max((deltat[(i-1),2] * aox[1,1] * bt[1]), (deltat[(i-1),3] * aox[1,2] * bt[2]))==(deltat[(i-1),2] * aox[1,1] * bt[1]), 0, 1)
	psit[i,3] <- ifelse(max((deltat[(i-1),3] * aox[2,2] * bt[2]), (deltat[(i-1),2] * aox[2,1] * bt[1]))==(deltat[(i-1),3] * aox[2,2] * bt[2]), 1, 0)
	}

#termination
QT <- matrix(data = NA, ncol = 2, nrow = length(psit[,1]))
QT[,1] <- pox[,1]
lpox <- length(pox[,1]) 

PT <- max(deltat[lpox,2:3])
QT[lpox,2] <- ifelse(deltat[lpox,2]== PT, psit[lpox,2], psit[lpox,3])

#path (state sequence) backtracking

for (i in (lpox-1):1){
	if(QT[(i+1),2] == 0){
		QT[i,2] <- ifelse(psit[i,2] == 0, 0, 1)
			} 
	if(QT[(i+1),2] == 1){
		QT[i,2] <- ifelse(psit[i,3] == 1, 1, 0)
			}
	}

plot(pox[,2], pch = 16, col = ifelse (QT[,2] == 1, "red", "blue"), main = paste(k, line, sep = "-"))
points(overallfit[,2], pch = 16, col = "black")

plot(pox3[,2] ~ pox3[,1], col = "green")
points(overallfit[,2] ~ overallfit[,1], pch=16)

lpeaks <- rle(QT[,2])$len[rle(QT[,2])$val == 1]
npeaks <- length(lpeaks)
ovol <- c(rep(NA, times= npeaks))
ovol <- c(rep(NA, times= npeaks))

for (i in 1:npeaks){
	
	peakpoints <- pox[(sum(rle(QT[,2])$lengths[1:(2*i-1)])+1):(sum(rle(QT[,2])$lengths[1:(2*i)])),1]
	fitted <- overallfit[overallfit[,1] %in% peakpoints,]	
	ovol[i] <- sum((overallfit[overallfit[,1] %in% peakpoints,2] - pox[,2][pox[,1] %in% peakpoints])*flowr/(1-(pox[,2][pox[,1] %in% peakpoints])))	
	
	}
markvalvol <- sort(abs(ovol), decreasing = T)[1:nmark]
ovol <- abs(c(ovol[abs(ovol) %in% markvalvol], rep(NA, times=npad)))

c(cvol, ovol)

}


# run peakarea function (chunk 2) on files identified by time file (chunk 1) (major chunk - 3 of 3)
nfiles <- length(time[,1])

#for file: file <- read.csv("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\3_13\\3_13N10_data")

#root directory - assumed to be constant across experiment - C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\
root <- "C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\"

#folders corresponding to days within root directory (if needed) (i.e. D1 == 3_13)
daycorr <- data.frame(index = paste("D",1:11, sep=""),data = c("3_13","3_14","3_15","3_16","3_17","3_18","7_8","7_25","3_31","4_1","4_2"))

#nuisance O2 - for the truly scary
nuisances <- c(3, 45, 59, 67, 93)

#oxygen buffer mod (k, truncation buffer)
obuffermod <- matrix(c(32, 15), ncol = 2, byrow = T)

#carbon-dioxide buffer mod (k, 5' truncation buffer, 3' truncation buffer)
cbuffermod <- matrix(c(126, 10, 85, 163, 10, 85), ncol = 3, byrow = T)

#excessive runin (removing initial points) (alternating ks & #number of removed points)
runin <- matrix(c(1, 100, 26, 200, 50, 30, 86, 15, 110, 30, 156, 600, 180, 900, 200, 200, 205, 120, 208, 1850, 257, 450, 261, 30), ncol = 2, byrow = T)

#missing lines / unsalvagably bad data
missingl <- c(27, 80, c(87:91), 93)
for (z in 1:length(time[,1])){
	if (rowSums(time[z,3:5]) == 0){
		missingl <- c(missingl, z)
		}
	}
#lines where some peaks should be used but others discarded (alternating ks, peak to remove (if two peaks are to be removed, use two entries))
partialread <- matrix(c(93, 2, 93, 3), ncol =2, byrow = T)	

vgas <- NULL
Rq <- NULL
test <- NULL

pdf("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\grandpdf.pdf", onefile=TRUE, bg="transparent", pointsize=12, width=8.5, height = 10.5)
par(mfrow = c(2,2))

for (k in 1:nfiles){
	
	line <- time[k,1]
	day <- time[k,2]
	date <- daycorr[,2][daycorr[,1] == day]
	directory <- paste(root, date, "\\", date, line, "_data", sep="")
	if(!(k %in% missingl)){
		file <- read.csv(directory)
		missingfile <- 0
		} else {missingfile <- 1}
	
	if(k %in% runin[,1]){
		file <- file[runin[,2][runin[,1] == k]:length(file[,1]),]
		file[,1] <- c(0:(length(file[,1])-1))
		}

	times <- time[k,3:5]
	nmark <- length(times[times != 0])
	times <- times[1:nmark]
	
	if (missingfile == 1){
		output <- rep(NA, times = nmax*4)
		} else{output <- peakarea(file,nmark)}
				
	cvolmin <- ifelse(is.na(output[1:3]), NA, (output[1:3]/times))
	ovolmin <- ifelse(is.na(output[4:6]), NA, (output[4:6]/times))
	output2 <- data.frame(line, day, cvolmin[1],cvolmin[2],cvolmin[3],ovolmin[1],ovolmin[2],ovolmin[3])
	names(output2) <- c("line", "day", "vco2R1", "vco2R2", "vco2R3", "vo2R1", "vo2R2", "vo2R3")
	
	vgas <- rbind(vgas, output2)
	names(vgas) <- c("line", "day", "vco2R1", "vco2R2", "vco2R3", "vo2R1", "vo2R2", "vo2R3")
	
	output3 <- data.frame(data = line, day, output2[,3:5]/output2[,6:8])
	names(output3) <- c("line", "day", "Rq1", "Rq2", "Rq3") 
	
	Rq <- rbind(Rq, output3)
	names(Rq) <- c("line", "day", "Rq1", "Rq2", "Rq3") 
	}

dev.off(which=dev.cur())

#correct for number of flies (some had four flies)
vgas[3,5] <- vgas[3,5]*5/4
vgas[3,8] <- vgas[3,8]*5/4
vgas[40,5] <- vgas[40,5]*5/4
vgas[40,8] <- vgas[40,8]*5/4
vgas[126,4] <- vgas[126,4]*5/4
vgas[126,7] <- vgas[126,7]*5/4
vgas[174,4] <- vgas[174,4]*5/4
vgas[174,7] <- vgas[174,7]*5/4

vgasfly <- cbind(vgas[,1],vgas[,2],vgas[,c(3:8)]/5)
write.csv(vgas, "C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\output8.1.10.csv") 