#input
#vgasfly
#vgasfly <- vgasfly[,-1]
vgasfly[98,1] <- as.factor("N10")
vgasfly[99,1] <- as.factor("N10")
vgasfly[193,1] <- as.factor("ZH33")
vgasfly[194,1] <- as.factor("ZH33")
vgasfly[199,1] <- as.factor("I24")
vgasfly[200,1] <- as.factor("I24")
#timest
#time
linepop <- read.table("C:\\Users\\Sean\\Desktop\\Respirometry\\Linepop.csv", sep = ",", header = TRUE)
rownames(linepop) <- linepop[,1]

#C1: line
#C2: population
#C3: day
#C4: timestart
#C5: timeincub
#C6: co2 moles
#C7: o2 moles

vgasflybypoint <- NULL

for (r in 1:length(vgasfly[,1])){
	nr <- 3 - (is.na(vgasfly[r,3]) + is.na(vgasfly[r,4]) + is.na(vgasfly[r,5]))
	
	if(nr != 0){for (s in 1:nr){
	
		C1rs <- as.factor(vgasfly[r,1])
		C2rs <- linepop[as.character(vgasfly[r,1]),2]
		C3rs <- vgasfly[r,2]
		C4rs <- timest[r,(2+s)]
		C5rs <- time[r, (2+s)]
		C6rs <- vgasfly[r,(2+s)]
		C7rs <- vgasfly[r,(5+s)]
		vgasflybypoint <- rbind(vgasflybypoint, data.frame(line = C1rs, population = C2rs, day = C3rs, timestart = C4rs, timeincub = C5rs, Vcot = C6rs, Vox = C7rs))
		}}}
#to account for subsampling
vgasflybypoint[,6] <- vgasflybypoint[,6] * 3.2/3
vgasflybypoint[,7] <- vgasflybypoint[,7] * 3.2/3

#write.csv(vgasflybypoint, "C:\\Users\\Sean\\Desktop\\Respirometry\\vgasfly9.26.10")
#vgasflybypoint <- read.csv("H:\\respgas8.2.10.csv")

inbredgas <- vgasflybypoint[vgasflybypoint$population %in% c("N", "I", "T", "B", "Z"),]
inbredgas[,3] <- as.factor(as.character(inbredgas[,3]))
inbredgas[,2] <- as.factor(as.character(inbredgas[,2]))
intrapopgas <- vgasflybypoint[vgasflybypoint$population %in% c("N", "NN", "B", "BB"),]
intrapopgas[,3] <- as.factor(as.character(intrapopgas[,3]))
intrapopgas[,2] <- as.factor(as.character(intrapopgas[,2]))
Ngas <- vgasflybypoint[vgasflybypoint$population %in% "N",]
Ngas[,1] <- as.factor(as.character(Ngas[,1]))

pdf("C:\\Users\\Sean\\Desktop\\Respirometry\\obesity_resp_march10\\outputplot.pdf", onefile=TRUE, bg="transparent", pointsize=12, width=8.5, height = 10.5)
par(mfrow = c(2,1))

#coverage		
hist(summary(vgasflybypoint$line, maxsum = 200, descending = TRUE), main = "coverage")

#boxplot of co2 production by day
boxplot(vgasflybypoint$Vcot ~ vgasflybypoint$day )

#circadian effects
summary(lm(vgasflybypoint$Vcot ~ vgasflybypoint$timestart))
plot(vgasflybypoint$Vcot ~ vgasflybypoint$timestart, main = "Carbon dioxide production by time of day")

#incubation differences

plot(vgasflybypoint$Vcot ~ vgasflybypoint$timeincub, main = "Carbon dioxide production by length of incubation")
summary(lm(vgasflybypoint$Vcot ~ vgasflybypoint$timeincub))
#remove intrapops
summary(lm(inbredgas$Vcot ~ inbredgas$timeincub))

#Inbred Lines

#Population boxplots and ANOVA
boxplot(inbredgas$Vcot ~ inbredgas$population, notch = TRUE, main = "Carbon dioxide production by population (inbred)")
anova(lm(inbredgas$Vcot ~ inbredgas$population))
boxplot(inbredgas$Vox ~ inbredgas$population, notch = TRUE, main = "Oxygen consumption by population (inbred)")
anova(lm(inbredgas$Vox ~ inbredgas$population))
boxplot((inbredgas$cot/inbredgas$Vox) ~ inbredgas$population, notch = TRUE, main = "Respiratory Quotients by population (inbred)")
anova(lm((inbredgas$cot/inbredgas$Vox)~ inbredgas$population))

#Line boxplots and ANOVA

Vcot <- inbredgas$Vcot
linec <- inbredgas$line
linec <- factor(linec, levels = names(sort(tapply(Vcot, linec, median, na.rm=TRUE))))
boxplot(Vcot  ~ linec, main = "Carbon dioxide production by line (inbred)", xlab = "Line", ylab = "L/min")
anova(lm(Vcot ~ linec))

Vox <- inbredgas$Vox
lineo <- inbredgas$line
lineo <- factor(lineo, levels = names(sort(tapply(Vox, lineo, median, na.rm=TRUE))))
boxplot(Vox  ~ lineo, main = "Oxygen consumption by line (inbred)", xlab = "Line", ylab = "L/min" )
anova(lm(Vox ~ lineo))

Rq <- inbredgas$Vcot/inbredgas$Vox
linerq <- inbredgas$line
linerq <- factor(linerq, levels = names(sort(tapply(Rq, linerq, median, na.rm=TRUE))))
boxplot(Rq  ~ linerq, main = "Respiratory Quotient by line (inbred)", xlab = "Line", ylab = "CO2 / Ox" )
anova(lm(Rq ~ linerq))

#Single Population

Vcot <- Ngas$Vcot
linec <- Ngas$line
linec <- factor(linec, levels = names(sort(tapply(Vcot, linec, median, na.rm=TRUE))))
boxplot(Vcot  ~ linec, main = "Carbon dioxide production by line (inbred)", xlab = "Line", ylab = "L/min")
anova(lm(Vcot ~ linec))

Vox <- Ngas$Vox
lineo <- Ngas$line
lineo <- factor(lineo, levels = names(sort(tapply(Vox, lineo, median, na.rm=TRUE))))
boxplot(Vox  ~ lineo, main = "Oxygen consumption by line (inbred)", xlab = "Line", ylab = "L/min" )
anova(lm(Vox ~ lineo))

Rq <- Ngas$Vcot/Ngas$Vox
linerq <- Ngas$line
linerq <- factor(linerq, levels = names(sort(tapply(Rq, linerq, median, na.rm=TRUE))))
boxplot(Rq  ~ linerq, main = "Respiratory Quotient by line (inbred)", xlab = "Line", ylab = "CO2 / Ox" )
anova(lm(Rq ~ linerq))

#Intrapopulation vs inbred

#Population boxplots and ANOVA
boxplot(intrapopgas$Vcot ~ intrapopgas$population, notch = TRUE, main = "Carbon dioxide production by population")
anova(lm(intrapopgas$Vcot ~ intrapopgas$population))
boxplot(intrapopgas$Vox ~ intrapopgas$population, notch = TRUE, main = "Oxygen consumption by population")
anova(lm(intrapopgas$Vox ~ intrapopgas$population))
boxplot((intrapopgas$cot/intrapopgas$Vox) ~ intrapopgas$population, notch = TRUE, main = "Respiratory Quotients by population")
anova(lm((intrapopgas$cot/intrapopgas$Vox)~ intrapopgas$population))

#Line boxplots and ANOVA

Vcot <- intrapopgas$Vcot
linec <- intrapopgas$line
linec <- factor(linec, levels = names(sort(tapply(Vcot, linec, median, na.rm=TRUE))))
sortedpop <- linepop[as.character(unique(sort(linec))),2]
boxplot(Vcot  ~ linec, main = "Carbon dioxide production by line", xlab = "Line", ylab = "L/min", col=ifelse(sortedpop %in% c("B","N"), "aquamarine", "darkgreen"))
anova(lm(Vcot ~ linec))

Vox <- intrapopgas$Vox
lineo <- intrapopgas$line
lineo <- factor(lineo, levels = names(sort(tapply(Vox, lineo, median, na.rm=TRUE))))
sortedpop <- linepop[as.character(unique(sort(lineo))),2]
boxplot(Vox  ~ lineo, main = "Oxygen consumption by line", xlab = "Line", ylab = "L/min", col=ifelse(sortedpop %in% c("B","N"), "aquamarine", "darkgreen"))
anova(lm(Vox ~ lineo))

Rq <- intrapopgas$Vcot/intrapopgas$Vox
linerq <- intrapopgas$line
linerq <- factor(linerq, levels = names(sort(tapply(Rq, linerq, median, na.rm=TRUE))))
sortedpop <- linepop[as.character(unique(sort(linerq))),2]
boxplot(Rq  ~ linerq, main = "Respiratory Quotient by line", xlab = "Line", ylab = "CO2 / Ox", col=ifelse(sortedpop %in% c("B","N"), "aquamarine", "darkgreen") )
anova(lm(Rq ~ linerq))

dev.off(which=dev.cur())


