#polynomial fitting requires libraries 
#mda - for polyreg


#break apart discrete trajectories or trajectory switching in a file by looking at extreme instantaneous acceleration

#conditions for a good flight trajectory
#traj include at least 10 points
#traj starts at expected takeoff level, z ~= 0.
#	deviation from norm w/ dist of z at time t=1.
#traj is restricted to ascending flight 
#	zt >= z(t-1)
#flight is captured
#	velocity is greater than a value adequate to distinguish walking and flight

#categorical function specifying the degree of a polynomial fitted to the flight trajectory (because number of points will differ)
scalar <- 2.5
nxlarray <- data.frame(length = c(10 + scalar*(0:50)), degree = c(3:53))
degree <- c(3:53)
Ndegree <- function(l){
	if(l < max(nxlarray[,1])){
	corval <- as.logical(0 <= l - nxlarray[1] & l - nxlarray[1] < scalar)
	degree[cumsum(rle(corval)$length)[rle(corval)$values == TRUE]]
	} else {nxlarray[,2][nxlarray[,1] == max(nxlarray[,1])]}
	}

#take the derivatives of a vector of polynomial coefficients: where c(a,b,c) are coefficients of the order (0, 1, 2)
coefderiv <- function(dispcoef){
	velcomp <- NULL
	for (p in 1:(length(dispcoef)-1)){
		velcomp <- c(velcomp, dispcoef[p+1] * p)
	}
	velcomp
	}

#use polynomial coefficients to evaluate the value at time t
fitcoef <- function(t, coef){
	margval <- NULL
	for (o in 1:(length(coef))){
		margval <- sum(margval, coef[o] * t^(o-1))
		}
	margval
	}

#for saving path information for QC
#for first 300 files
#partialsummary <- data.frame(summarystats[c(1:300),5],summarystats[c(1:300),6],summarystats[c(1:300),7])
#names(partialsummary) <- c("filenum", "trajnum", "chunk")

#fastflight <- sort(summarystats[,10], decreasing = TRUE)[100]
#fastflights <- summarystats[,10] >= fastflight
#fastsummary <- data.frame(summarystats[fastflights,5],summarystats[fastflights,6],summarystats[fastflights,7])
#names(fastsummary) <- c("filenum", "trajnum", "chunk")

#slowflight <- sort(summarystats[,10])[100]
#slowflights <- summarystats[,10] <= slowflight
#slowsummary <- data.frame(summarystats[slowflights,5],summarystats[slowflights,6],summarystats[slowflights,7])
#names(slowsummary) <- c("filenum", "trajnum", "chunk")

inputfiles2 <- read.csv("C:\\Users\\Sean\\Desktop\\Flight_traj\\flightflielocalcsv.csv")
inputfiles2 <- cbind(inputfiles2, rep(NA, times = length(inputfiles2[,1]))) 

#minimum # of points for a flight
nmin <- 10

#experimentally determined unusable flights
wackflight <- data.frame(k = c(78, 150, 169, 342), i = c(1, 7, 6, 1), comp = rep(NA, times =2))
wackflight[,3] <- paste(wackflight[,1], wackflight[,2], sep =".")

#flight with a long track resulting in a small usable portion - truncate
truncflight <- data.frame(k = c(0), i = c(0), to = c(0), tf = c(0), comp = rep(NA, times=1))
truncflight[,5] <- paste(truncflight[,1], truncflight[,2], sep =".")



#summary statistics
#total velocity - sqrt(vx ^2 + vy ^2 + vz ^2)
#meanv, sdv

#tangential acceleration - (ax * vx + ay * vy + az * vz) / v
#meanaccelt, sdaccelt

#centripetal acceleration - sqrt((ax - accelt*vx/v)^2 + (ay - accelt*vy/v)^2 + (az - accelt*vz/v)^2)
#meanaccelc, sdaccelc

#angular acceleration - accelc/v
#meanaccela, sdaccela

#takeoff angle between points 3 and 8 - arctan(dz / sqrt(dx^2 + dy^2))
#toangle

summarystats <- NULL
ntraj <- 0
fntraj <- 0
sntraj <- 0

for (k in 1:length(inputfiles2[,1])){
directory <- paste("C:\\Users\\Sean\\Desktop\\Flight_traj\\",inputfiles2[k,3],"\\",inputfiles2[k,2], sep = "")
activesheet <- read.table(directory)

n <- length(activesheet[1,])/3

activesheet <- cbind(c(1: length(activesheet[,1])), activesheet)

for (i in 1:n){
#initial values are in mm/measurement == mm / 1/120 s	
	
	instposi <- activesheet[,c(1,(3*i-1):(3*i+1))]
	posi <- as.matrix(instposi[(apply(instposi, 1, sum)!=instposi[,1]),])
		
	#to weed out bad files - truncate them to 4 points so they're removed as short trajectories
	index <- paste(k, i, sep=".")
	if(index %in% wackflight[,3]){posi <- data.frame(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(1,2,3,4))} else { posi[,1] <- posi[,1] - posi[1,1] }
	

	#truncate flight according to truncflight
	if(index %in% truncflight[,5]){posi <- posi[posi[,1] %in% (truncflight[,3][truncflight[,5] %in% index]:truncflight[,4][truncflight[,5] %in% index]),]}

	#look for extreme acceleration (>5) indicating a failure to track a single flight

	vinst <- apply(posi, 2, diff)
	vinsttot <- sqrt(vinst[,2]^2 + vinst[,3]^2 + vinst[,4]^2)
	ainst <- apply(vinst, 2, diff)
	ainsttot <- sqrt(ainst[,2]^2+ainst[,3]^2+ainst[,4]^2)
	speedbreak <- data.frame(higha = rle(ainsttot > 5)$values, length = rle(ainsttot > 5)$lengths)

	#intervals to compare for best signiture of flight
	
	nchunk <- length(speedbreak$length[speedbreak$higha == FALSE])
	short <- rep(NA, times = nchunk)
	wrongstart <- rep(NA, times = nchunk)
	
	#data will be added for unambigous flights if short = F and wrongstart = F


	for (m in 1:nchunk){
	
	#implement a skip based upon currated data (i.e. chunk 2 is flight, 1 is garbage)

	#determine subinterval
		if(speedbreak$higha[1] == FALSE){
		if(m == 1){ 
			interval <- c(1:(speedbreak$length[1]+2))
			} else {
				interval <- c((cumsum(speedbreak$length)[((m-1)*2)]+1):(cumsum(speedbreak$length)[(m*2)-1]+2))
				}
		} else {
			interval <- c((cumsum(speedbreak$length)[(m*2)-1]+1):(cumsum(speedbreak$length)[m*2] + 2))
			}

	posisub <- posi[interval,]
	
	if(length(posisub[,1]) < 10){short[m] <- 1
	next}
	
	#categorically specify N
	N <- Ndegree(length(posisub[,1]))
	polyx <- polyreg(posisub[,1], posisub[,2], degree = N)$fitted.values
	polyy <- polyreg(posisub[,1], posisub[,3], degree = N)$fitted.values	
	polyz <- polyreg(posisub[,1], posisub[,4], degree = N)$fitted.values	
	polyzm <- cbind(posisub[,1], polyx, polyy, polyz)


	#select for suffient movement to indicate flight, > 1 mm/120s
	
	if (length(posi[,4]) < nmin){short[i] <- TRUE} else {
	short[i] <- FALSE	

		vinst <- apply(polyzm, 2, diff)
		vinsttot <- sqrt(vinst[,2]^2+vinst[,3]^2+vinst[,4]^2)
		fastmov <- vinsttot > 1

		#check for a single slow element -> next/skip
		if((length(rle(fastmov)$lengths) == 1) & (rle(fastmov)$values[1] == FALSE)){fastpoints <- 1} else{

		rlefast <- rle(fastmov)
 		fastmoving <- data.frame(count <- c(1: length(rlefast$values)), values = rlefast$values, lengths = rlefast$lengths, cumsuml = cumsum(rlefast$lengths))
		
		maxn <- max(rlefast$lengths[rlefast$values == TRUE])
		int <- fastmoving[,1][fastmoving[,2] == TRUE & fastmoving[,3] == maxn]

		if(int == 1){fastpoints <- c(1:(fastmoving[1,4]+1))} else {fastpoints <- fastmoving[(int-1),4]:(fastmoving[int,4]+1)}
		
		posisub <- matrix(posisub[posisub[,1] %in% posisub[fastpoints,1]], ncol=4)
	}}

	#select for upwards flight

	if (length(posisub[,4]) < nmin){short[i] <- TRUE} else {
	short[i] <- FALSE

		ascend <- diff(posisub[,4])>0
		ascending <- NULL
		ascending[1] <- ifelse(rle(ascend)$values[1] == TRUE, 1, 0) 
		
		#check for pure descent
		if((length(rle(ascend)$lengths) == 1 & (rle(ascend)$values[1] == FALSE))){ascending <- 1} else{
		
		if(ascending[1] == 1){
			ascending <- 1:(rle(ascend)$lengths[1]+1)		
		} else {
			ascending <- c((cumsum(rle(ascend)$lengths)[1]+1):(cumsum(rle(ascend)$lengths)[2]+1))
			}}
			
		posisub <- matrix(posisub[posisub[,1] %in% posisub[ascending,1]], ncol=4)
		posisub[,1] <- (ascending - ascending[1])
	}
	
	if(length(posisub[,1]) < 10){short[i] <- TRUE}

	 # check that initial height is between -40 & 20 mm
	
	#ifelse((posisub[1,4] > -40) & (posisub[1,4] < 20), wrongstart[i] <- FALSE, wrongstart[i] <- TRUE)

	if(short[i] == FALSE){
			
	
	#units still in mm
				
		N <- Ndegree(length(posisub[,1]))
		
		polyopx <- polyreg(posisub[,1], posisub[,2], degree = N)
		polyopy <- polyreg(posisub[,1], posisub[,3], degree = N)
		polyopz <- polyreg(posisub[,1], posisub[,4], degree = N)
		
		polyx <- polyopx$fitted.values
		polyy <- polyopy$fitted.values	
		polyz <- polyopz$fitted.values	
		polyxyz <- cbind(posisub[,1], polyx, polyy, polyz)
		
		polyxcoef <- polyopx$coef
		polyycoef <- polyopy$coef	
		polyzcoef <- polyopz$coef
		
		polyxyzcoef <- cbind(polyxcoef, polyycoef, polyzcoef)				
		
		while(sum(is.na(polyxyzcoef)) != 0){
			nacoef <- is.na(polyxyzcoef)
			nacoef2 <- data.frame(degree = c(0:(length(nacoef[,1])-1)), value = rep(NA, times = (length(nacoef[,1]))))
			for (j in 1:length(nacoef[,1]))
				{
				if(sum(nacoef[j,]) == 0){nacoef2[j,2] <- 0} else{nacoef2[j,2] <- 1}
				}
			N <- (min(nacoef2[,1][nacoef2[,2] == 1]) -1 )
			
			polyopx <- polyreg(posisub[,1], posisub[,2], degree = N)
			polyopy <- polyreg(posisub[,1], posisub[,3], degree = N)
			polyopz <- polyreg(posisub[,1], posisub[,4], degree = N)
		
			polyx <- polyopx$fitted.values
			polyy <- polyopy$fitted.values	
			polyz <- polyopz$fitted.values	
			polyxyz <- cbind(posisub[,1], polyx, polyy, polyz)
		
			polyxcoef <- polyopx$coef
			polyycoef <- polyopy$coef	
			polyzcoef <- polyopz$coef
			
			polyxyzcoef <- cbind(polyxcoef, polyycoef, polyzcoef)

			}
			
		posxyzfit <- matrix(data = NA, ncol = 4, nrow = length(posisub[,1]))
		posxyzfit[,1] <- posisub[,1]
		for (j in 1:3){
		for (o in 1:length(posisub[,1])){
			posxyzfit[o,(j+1)] <- fitcoef(posisub[o,1], polyxyzcoef[,j])
			}
			}		
		
		velxyzcoef <- apply(polyxyzcoef, 2, coefderiv)
		
		velxyzfit <- matrix(data = NA, ncol = 4, nrow = length(posisub[,1]))
		velxyzfit[,1] <- posisub[,1]
		 
		for (j in 1:3){
		for (o in 1:length(posisub[,1])){
			velxyzfit[o,(j+1)] <- fitcoef(posisub[o,1], velxyzcoef[,j])
			}
			}
		
		#convert units of mm / 120th s to m/s
		velxyzfit[,c(2:4)] <- velxyzfit[,c(2:4)] * 120/1000

		veltot <- sqrt(velxyzfit[,2] ^2 + velxyzfit[,3]^2 + velxyzfit[,4]^2)


		accelxyzcoef <- apply(velxyzcoef, 2, coefderiv)
	
		accelxyzfit <- matrix(data = NA, ncol = 4, nrow = length(posisub[,1]))
		accelxyzfit[,1] <- posisub[,1]

		for (j in 1:3){
		for (o in 1:length(posisub[,1])){
			accelxyzfit[o,(j+1)] <- fitcoef(posisub[o,1], accelxyzcoef[,j])
			}
			}		
		
		accelxyzfit[,c(2:4)] <- accelxyzfit[,c(2:4)] * (120^2)/1000

		accelt <- (velxyzfit[,2]*accelxyzfit[,2] + velxyzfit[,3]*accelxyzfit[,3] + velxyzfit[,4]*accelxyzfit[,4])/veltot

		accelc <- sqrt((accelxyzfit[,2] - (accelt/veltot)*velxyzfit[,2])^2 + (accelxyzfit[,3] - (accelt/veltot)*velxyzfit[,3])^2 + (accelxyzfit[,4] - (accelt/veltot)*velxyzfit[,4])^2)
		
		accela <- accelc/veltot 
		
		#arctan(dz / sqrt(dx^2 + dy^2))
		takeoffangle <- atan2((posxyzfit[8,4] -  posxyzfit[3,4]),sqrt((posxyzfit[8,2] -  posxyzfit[3,2])^2 + (posxyzfit[8,3] -  posxyzfit[3,3])^2))
		
		#arctan(dy/dx)
		todirbias <- atan2(posxyzfit[8,3] - posxyzfit[3,3], posxyzfit[8,2] - posxyzfit[3,2])

		#crop first and last two points
		cropbuf <- 2
		croppedpos <- (1:length(posxyzfit[,1]))[!((1:length(posxyzfit[,1])) %in% c((1:cropbuf), (length(posxyzfit[,1]) - cropbuf + 1): length(posxyzfit[,1])))]
		
		#summary statistics
		#meanvelocity, sdvelocity
		meanv <- mean(veltot[croppedpos])
		sdv <- sd(veltot[croppedpos])
		
		meanVz <- mean(velxyzfit[,4][croppedpos])
		sdVz <- sd(velxyzfit[,4][croppedpos])

		#meanaccelt, sdaccelt - tangential acceleration
		meanaccelt <- mean(abs(accelt[croppedpos]))
		sdaccelt <- sd(abs(accelt[croppedpos]))

		#meanaccelc, sdaccelc - centrifugal/radial acceleration
		meanaccelc <- mean(abs(accelc[croppedpos]))
		sdaccelc <- sd(abs(accelc[croppedpos]))

		#meanaccela, sdaccela - angular acceleration
		meanaccela <- mean(abs(accela[croppedpos]))
		sdaccela <- sd(abs(accela[croppedpos]))

		runoutput <- as.data.frame(c(inputfiles2[k,1:3], filenum = k, traj = i, chunk = m, polydeg = N, points = length(posisub[,1]), meanv = meanv, sdv = sdv, meanaccelt = meanaccelt, sdaccelt = sdaccelt, meanaccelc = meanaccelc, sdaccelc = sdaccelc, meanaccela = meanaccela, sdaccela = sdaccela, toangle = takeoffangle, meanVz = meanVz, sdVz = sdVz))
		
		#max stats
		#maximum interval - minus first 2 and last 2 points
		npoints <- length(posisub[,1]) - 4
		int <- (3*100):((npoints+2)*100)/100
		
		#max velocity
		mvelxyzfit <- matrix(data = NA, ncol = 4, nrow = length(int))
		mvelxyzfit[,1] <- int
		 
		for (j in 1:3){
		for (o in 1:length(int)){
			mvelxyzfit[o,(j+1)] <- fitcoef(mvelxyzfit[o,1], velxyzcoef[,j])
			}
			}
		
		#convert units of mm / 120th s to m/s
		mvelxyzfit[,c(2:4)] <- mvelxyzfit[,c(2:4)] * 120/1000

		mveltot <- sqrt(mvelxyzfit[,2] ^2 + mvelxyzfit[,3]^2 + mvelxyzfit[,4]^2)
		maxmvel <- max(mveltot)		
		#time of max v
		mvtime <- (mvelxyzfit[,1][mveltot == maxmvel]-2)/npoints
		
		mvely <- mvelxyzfit[,4]
		maxmvely <- max(mvely)
		mvytime <- (mvelxyzfit[,1][mvely == maxmvely]-2)/npoints

		#max acceleration
				
		maccelxyzfit <- matrix(data = NA, ncol = 4, nrow = length(int))
		maccelxyzfit[,1] <- int

		for (j in 1:3){
		for (o in 1:length(int)){
			maccelxyzfit[o,(j+1)] <- fitcoef(maccelxyzfit[o,1], accelxyzcoef[,j])
			}
			}		
		
		maccelxyzfit[,c(2:4)] <- maccelxyzfit[,c(2:4)] * (120^2)/1000

		maccelt <- (mvelxyzfit[,2]*maccelxyzfit[,2] + mvelxyzfit[,3]*maccelxyzfit[,3] + mvelxyzfit[,4]*maccelxyzfit[,4])/mveltot
		maxaccelt <- max(abs(maccelt))
		mattime <- (maccelxyzfit[,1][abs(maccelt) == maxaccelt]-2)/npoints

		maccelc <- sqrt((maccelxyzfit[,2] - (maccelt/mveltot)*mvelxyzfit[,2])^2 + (maccelxyzfit[,3] - (maccelt/mveltot)*mvelxyzfit[,3])^2 + (maccelxyzfit[,4] - (maccelt/mveltot)*mvelxyzfit[,4])^2)
		maxaccelc <- max(abs(maccelc))
		mactime <- (maccelxyzfit[,1][abs(maccelc) == maxaccelc]-2)/npoints

		maccela <- maccelc/mveltot
		maxaccela <- max(abs(maccela))
		maatime <- (maccelxyzfit[,1][abs(maccela) == maxaccela]-2)/npoints

		maxstats <- data.frame(maxmeanv = maxmvel, tmaxmeanv = mvtime, maxvy = maxmvely, tmaxvy = mvytime, maxat = maxaccelt, tmaxat = mattime, maxac = maxaccelc, tmaxac = mactime, maxaa = maxaccela, tmaxaa = maatime)		

		summarystats <- rbind(summarystats, cbind(runoutput, maxstats))
				
	#k = file	
	#i = traj
	#m = chunk
	#kpartial <- partialsummary[partialsummary[,1] == k,]
	#ipartial <- kpartial[kpartial[,2] == i,]
	#mpartial <- ipartial[ipartial[,3] == m,]

	#kfast <- fastsummary[fastsummary[,1] == k,]
	#ifast <- kfast[kfast[,2] == i,]
	#mfast <- ifast[ifast[,3] == m,]

	#kslow <- slowsummary[slowsummary[,1] == k,]
	#islow <- kslow[kslow[,2] == i,]
	#mslow <- islow[islow[,3] == m,]

	#if(k %in% mpartial[1,1] & i %in% mpartial[1,2] & m %in% mpartial[1,3]){
	#	ntraj <- ntraj + 1
	#	trajname <- paste("fly", ntraj, sep ="")
	#	write.csv(cbind(posisub, posxyzfit), paste("C:\\Users\\Sean\\Desktop\\Flight_traj\\tempwrite\\",trajname, ".csv", sep =""))
	#	}

	#if(k %in% mfast[1,1] & i %in% mfast[1,2] & m %in% mfast[1,3]){
	#	fntraj <- fntraj + 1
	#	trajname <- paste("fast", fntraj, sep ="")
	#	write.csv(cbind(posisub, posxyzfit), paste("C:\\Users\\Sean\\Desktop\\Flight_traj\\tempwrite\\",trajname, ".csv", sep =""))
	#	}

	#if(k %in% mslow[1,1] & i %in% mslow[1,2] & m %in% mslow[1,3]){
	#	sntraj <- sntraj + 1
	#	trajname <- paste("slow", sntraj, sep ="")
	#	write.csv(cbind(posisub, posxyzfit), paste("C:\\Users\\Sean\\Desktop\\Flight_traj\\tempwrite\\",trajname, ".csv", sep =""))
	#	}

	

		}
	}
}
} 
summarystats <- cbind(c(1:length(summarystats[,1])), summarystats)

#for trajectories where two flights were called remove the flight with the lower (usually meaningless) velocity 
kidup <- names(sort(summary(as.factor(paste(summarystats[,5],summarystats[,6], sep=".")), maxsum = length(summarystats[,4])), decreasing = TRUE)[sort(summary(as.factor(paste(summarystats[,5],summarystats[,6], sep=".")), maxsum = length(summarystats[,4])), decreasing = TRUE) > 1])
badchunktraj <- NULL
goodchunktraj <- NULL

for (j in 1:length(kidup)){
	tempsum <- summarystats[(paste(summarystats[,5], summarystats[,6], sep=".") %in% kidup[j]),]
	badchunktraj <- c(badchunktraj, tempsum[!tempsum[,10] == max(tempsum[,10]),][,1])
	goodchunktraj <- c(goodchunktraj, tempsum[tempsum[,10] == max(tempsum[,10]),][,1])
	}

summarystats <- summarystats[-badchunktraj,]
#summarystatsval <- summarystats


#write.table(summarystatsval, "C:\\Users\\Sean\\Desktop\\Flight_traj\\summarystats72110.csv")
summarystats

