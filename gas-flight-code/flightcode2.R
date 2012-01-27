#specify input file
#inputfiles2 is line name, directory, filename

#specify number of trajectories tracked
n <- length(activesheet[1,])/3
#specify line name and vial number
line <- "N01.1"
#data frame with line/vial, fly number, and instantaneous V
lineinstvel <- data.frame(rep(line,times=1000),flynumber=2000,instvel=1000)
names(lineinstvel) <- c("line.vial", "flynumber", "instvel") 

veln <- matrix(data=NA, ncol=4, nrow=n)
acceln <- matrix(data=NA, ncol=4, nrow=n)

for (i in 1:n){
	intposi <- editpos[,(3*i-2):(3*i)]
	posi <- as.matrix(intposi[(apply(intposi, 1, sum)!=0),])
	#change in position between continuous time points, instantaneous velocity
	
	# determine if there are more then ten points
	if(length(posi[,1])>10){
		
	vinst <- apply(posi, 2, diff)
	#total displacement
	initialpos <- matrix(rep(posi[1,1:3],length(posi[,1])), ncol=3, nrow=length(posi[,1]), byrow=TRUE)
	dispti <- posi - initialpos
	#point displacement, instantaneous v
	vinsttot <- sqrt(vinst[,1]^2+vinst[,2]^2+vinst[,3]^2)
	sdvtot <- sd(vinsttot)
	veln[i,] <- c(length(posi[,1]),(mean(vinsttot)-2*sdvtot), median(vinsttot),(mean(vinsttot)+2*sdvtot))
	
	lineinstvel[Xcount:(Xcount+length(vinsttot)-1),2] <- i
	lineinstvel[Xcount:(Xcount+length(vinsttot)-1),3] <- vinsttot
	Xcount <- (Xcount+length(vinsttot))
	#plot(density(vinsttot))
	
	#instantaneous acceleration
	ainst <- diff(vinst)
	ainsttot <- sqrt(ainst[,1]^2+ainst[,2]^2+ainst[,3]^2)
	sdatot <- sd(ainsttot)
	acceln[i,] <- c(length(posi[,1]),(mean(ainsttot)-2*sdatot), median(ainsttot),(mean(ainsttot)+2*sdatot))
	}else{
		veln[i,] <- c(length(posi[,1]),rep(NA, times=3))
		acceln[i,] <- c(length(posi[,1]),rep(NA, times=3))
		}

	}
lineinstvel <- lineinstvel[lineinstvel$flynumber != 2000,]	
validsampes <- unique(lineinstvel[,2])


boxplot(instvel~as.factor(flynumber), data = lineinstvel, notch = T)
plot(density(lineinstvel[lineinstvel$flynumber==1,3]))
lines(density(lineinstvel[lineinstvel$flynumber==2,3]), col="RED")
