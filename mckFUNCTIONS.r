resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

#Some of these packages may be needed
#library(stringr)
#library(plyr)
#library(RcppArmadillo)
#library(Rcpp)
#library(extrafont)
#loadfonts()
#multivariate normal density
#sourceCpp("dMvn.cpp")

## set home computer indicator
if(length(grep('mck14',getwd()))==1){
  computer <- 'duke'
}
if(length(grep('mckwit',getwd()))==1){
  computer <- 'home'
}

###list all functions
in_mckFUNCTIONS <- function(){
  if(computer=='home'){
    readLines("/home/mckwit/Dropbox/mckFUNCTIONS.r")[grep("[\\s]*<-[\\s]*function",readLines("/home/mckwit/Dropbox/mckFUNCTIONS.r"),perl=T)]
  }else{
    readLines("C:/Users/mck14/Dropbox/mckFUNCTIONS.r")[grep("[\\s]*<-[\\s]*function",readLines("C:/Users/mck14/Dropbox/mckFUNCTIONS.r"),perl=T)]
  }
}

####
#running weighted average to fill NAs
####

w.avg <- function(y){ #running weighted average to fill NAs
  whr.val <- which(!is.na(y)) 
  whr.na <- which(is.na(y))
  
  if(is.na(y[1]))y[1] <- y[min(whr.val)]
  if(is.na(y[length(y)]))y[length(y)] <- y[max(whr.val)]
  
  whr.val <- which(!is.na(y)) 
  whr.na <- which(is.na(y))
  
  for(i in whr.na){
    bel.x <- whr.val[max(which(whr.val < i))]
    abv.x <- whr.val[min(which(whr.val > i))]
    bel.y <- y[bel.x]
    abv.y <- y[abv.x]
    
    tot.dist <- ((i-bel.x)+(abv.x-i))
    y[i] <- (bel.y*(1-(i-bel.x)/tot.dist)) + (abv.y*(1-(abv.x-i)/tot.dist))
    
  }
  
  y
  
}


##
#fills gaps with left most value or zero warm grow
##

stuff.gap <- function(y){ #fills gaps
  pow <- y  
  n <- dim(pow)[1]
  if(sum(n)!=0){
    for(j in 1:n){
      vals <- which(!is.na(pow[j,]))
      nas <- which(is.na(pow[j,]))
  
      if(length(nas)==0)next
      if(length(vals)==0){pow[j,] <- 0;next}
      for(i in nas){
        if(i < min(vals)){ pow[j,i] <- 0
        }else{pow[j,i] <- pow[j,i-1]
        }
      }
   }
  }
  if(sum(n)==0){
    vals <- which(!is.na(pow))
    nas <- which(is.na(pow))
       
    for(i in nas){
      if(i < min(vals)){ pow[i] <- 0
      }else{pow[i] <- pow[i-1]
      }
    }
  }
  pow
  
}

####  
# random section of matrix random version of head() or tail()
####

rand <- function(mat,n=6){ # random version of head
  if(is.null(dim(mat))){
    seed <- sample(1:len(mat),1)
    xx <- mat[seed:(seed+n)]
  }else{
    seed <- sample(1:dim(mat)[1],1)
    xx <- mat[seed:(seed+n),]
    rownames(xx) <- seed:(seed+n)
  }
  xx
} 

###sets path to dropbox on either computer

path.dir <- function(working.dir = T){ #sets path to dropbox
  
  if(length(grep('mck14',getwd()))==1){
    path <- 'C:/Users/mck14/Dropbox/'
  }
  if(length(grep('mckwit',getwd()))==1){
    path <- '/home/mckwit/Dropbox/'
  }
  path <- paste(path,sep="")
  if(working.dir==T){setwd(path)}
}

####NOT IN

'%!in%' <- function(x,y)!('%in%'(x,y)) #not in

########

#
#text to columns
#

TtoC <- function(col,sep='/'){ #text to columns
  n <- length(col)
  vect <- unlist(strsplit(col,sep))
  m <- length(vect)
  matrix(vect,ncol=(m/n),byrow=T)
}

#### takes plot gives treatment
#DUKE FOREST
#G01,G04,G07,S03,S04,S08  Ambient
#G03,G05,G09,S02,S05,S09  3 degrees
#G02,G06,G08,S01,S06,S07  5 degrees
#G10,G11,G12,S10,S11,S12  Control

#HARVARD FOREST
#G02,G04,G07,S02,S05,S07
#G01,G06,G08,S03,S06,S08
#G03,G05,G09,S01,S04,S09
#G10,G11,G12,S10,S11,S12

treatment <- function(plots,combine=TRUE,site='DF',light = FALSE,combineAandC = FALSE){ #applies WS treatments
  treat <- plots
  if(site=='DF'){
  treat[which(plots=='S03'|plots=='S3'|plots=='s03'|plots=='s3')]<- 'zero'
  treat[which(plots=='S04'|plots=='S4'|plots=='s04'|plots=='s4')]<- 'zero'
  treat[which(plots=='S08'|plots=='S8'|plots=='s08'|plots=='s8')] <- 'zero'
  treat[which(plots=='S02'|plots=='S2'|plots=='s02'|plots=='s2')] <- 'three'
  treat[which(plots=='S05'|plots=='S5'|plots=='s05'|plots=='s5')] <- 'three'
  treat[which(plots=='S09'|plots=='S9'|plots=='s09'|plots=='s9')] <- 'three'
  treat[which(plots=='S01'|plots=='S1'|plots=='s01'|plots=='s1')] <- 'five'
  treat[which(plots=='S06'|plots=='S6'|plots=='s06'|plots=='s6')] <- 'five'
  treat[which(plots=='S07'|plots=='S7'|plots=='s07'|plots=='s7')] <- 'five'
  treat[which(plots=='S10'|plots=='s10')] <- 'control'
  treat[which(plots=='S11'|plots=='s11')] <- 'control'
  treat[which(plots=='S12'|plots=='s12')] <- 'control'

  treat[which(plots=='G01'|plots=='G1'|plots=='g01'|plots=='g1')] <- 'zero'
  treat[which(plots=='G04'|plots=='G4'|plots=='g04'|plots=='g4')] <- 'zero'
  treat[which(plots=='G07'|plots=='G7'|plots=='g07'|plots=='g7')] <- 'zero'
  treat[which(plots=='G03'|plots=='G3'|plots=='g03'|plots=='g3')] <- 'three'
  treat[which(plots=='G05'|plots=='G5'|plots=='g05'|plots=='g5')] <- 'three'
  treat[which(plots=='G09'|plots=='G9'|plots=='g09'|plots=='g9')] <- 'three'
  treat[which(plots=='G02'|plots=='G2'|plots=='g02'|plots=='g2')] <- 'five'
  treat[which(plots=='G06'|plots=='G6'|plots=='g06'|plots=='g6')] <- 'five'
  treat[which(plots=='G08'|plots=='G8'|plots=='g08'|plots=='g8')] <- 'five'
  treat[which(plots=='G10'|plots=='g10')] <- 'control'
  treat[which(plots=='G11'|plots=='g11')] <- 'control'
  treat[which(plots=='G12'|plots=='g12')] <- 'control'
  }
  #HARVARD FOREST
  #G02,G04,G07,S02,S05,S07 Ambient
  #G01,G06,G08,S03,S06,S08 3 degrees
  #G03,G05,G09,S01,S04,S09 5 degrees
  #G10,G11,G12,S10,S11,S12 Control
  if(site=='HF'){
    treat[which(plots=='S02'|plots=='S2'|plots=='s02'|plots=='s2')]<- 'zero'
    treat[which(plots=='S05'|plots=='S5'|plots=='s05'|plots=='s5')]<- 'zero'
    treat[which(plots=='S07'|plots=='S7'|plots=='s07'|plots=='s7')] <- 'zero'
    treat[which(plots=='S03'|plots=='S3'|plots=='s03'|plots=='s3')] <- 'three'
    treat[which(plots=='S06'|plots=='S6'|plots=='s06'|plots=='s6')] <- 'three'
    treat[which(plots=='S08'|plots=='S8'|plots=='s08'|plots=='s8')] <- 'three'
    treat[which(plots=='S01'|plots=='S1'|plots=='s01'|plots=='s1')] <- 'five'
    treat[which(plots=='S04'|plots=='S4'|plots=='s04'|plots=='s4')] <- 'five'
    treat[which(plots=='S09'|plots=='S9'|plots=='s09'|plots=='s9')] <- 'five'
    treat[which(plots=='S10'|plots=='s10')] <- 'control'
    treat[which(plots=='S11'|plots=='s11')] <- 'control'
    treat[which(plots=='S12'|plots=='s12')] <- 'control'
    
    treat[which(plots=='G02'|plots=='G2'|plots=='g02'|plots=='g2')] <- 'zero'
    treat[which(plots=='G04'|plots=='G4'|plots=='g04'|plots=='g4')] <- 'zero'
    treat[which(plots=='G07'|plots=='G7'|plots=='g07'|plots=='g7')] <- 'zero'
    treat[which(plots=='G03'|plots=='G3'|plots=='g03'|plots=='g3')] <- 'five'
    treat[which(plots=='G05'|plots=='G5'|plots=='g05'|plots=='g5')] <- 'five'
    treat[which(plots=='G09'|plots=='G9'|plots=='g09'|plots=='g9')] <- 'five'
    treat[which(plots=='G01'|plots=='G1'|plots=='g01'|plots=='g1')] <- 'three'
    treat[which(plots=='G06'|plots=='G6'|plots=='g06'|plots=='g6')] <- 'three'
    treat[which(plots=='G08'|plots=='G8'|plots=='g08'|plots=='g8')] <- 'three'
    treat[which(plots=='G10'|plots=='g10')] <- 'control'
    treat[which(plots=='G11'|plots=='g11')] <- 'control'
    treat[which(plots=='G12'|plots=='g12')] <- 'control'
  }
  
if(combine){
  if(combineAandC){treat[which(treat=='zero'|treat=='control')]<- 'A'}
      treat[which(treat=='three'|treat=='five')] <- 'H'
}
if(light){
  treat[grep("s",plots,ignore.case=T)] = paste("Shade",treat[grep("s",plots,ignore.case=T)],sep="")
  treat[grep("g",plots,ignore.case=T)] = paste("Gap",treat[grep("g",plots,ignore.case=T)],sep="")
}
treat
}


STDplots <- function(plots){ #standardize plot names for WS
  treat <- plots
  treat[which(plots=='S03'|plots=='S3'|plots=='s03'|plots=='s3')]<- 's3'
  treat[which(plots=='S04'|plots=='S4'|plots=='s04'|plots=='s4')]<- 's4'
  treat[which(plots=='S08'|plots=='S8'|plots=='s08'|plots=='s8')] <- 's8'
  treat[which(plots=='S02'|plots=='S2'|plots=='s02'|plots=='s2')] <- 's2'
  treat[which(plots=='S05'|plots=='S5'|plots=='s05'|plots=='s5')] <- 's5'
  treat[which(plots=='S09'|plots=='S9'|plots=='s09'|plots=='s9')] <- 's9'
  treat[which(plots=='S01'|plots=='S1'|plots=='s01'|plots=='s1')] <- 's1'
  treat[which(plots=='S06'|plots=='S6'|plots=='s06'|plots=='s6')] <- 's6'
  treat[which(plots=='S07'|plots=='S7'|plots=='s07'|plots=='s7')] <- 's7'
  treat[which(plots=='S10'|plots=='s10')] <- 's10'
  treat[which(plots=='S11'|plots=='s11')] <- 's11'
  treat[which(plots=='S12'|plots=='s12')] <- 's12'
  
  treat[which(plots=='G01'|plots=='G1'|plots=='g01'|plots=='g1')] <- 'g1'
  treat[which(plots=='G04'|plots=='G4'|plots=='g04'|plots=='g4')] <- 'g4'
  treat[which(plots=='G07'|plots=='G7'|plots=='g07'|plots=='g7')] <- 'g7'
  treat[which(plots=='G03'|plots=='G3'|plots=='g03'|plots=='g3')] <- 'g3'
  treat[which(plots=='G05'|plots=='G5'|plots=='g05'|plots=='g5')] <- 'g5'
  treat[which(plots=='G09'|plots=='G9'|plots=='g09'|plots=='g9')] <- 'g9'
  treat[which(plots=='G02'|plots=='G2'|plots=='g02'|plots=='g2')] <- 'g2'
  treat[which(plots=='G06'|plots=='G6'|plots=='g06'|plots=='g6')] <- 'g6'
  treat[which(plots=='G08'|plots=='G8'|plots=='g08'|plots=='g8')] <- 'g8'
  treat[which(plots=='G10'|plots=='g10')] <- 'g10'
  treat[which(plots=='G11'|plots=='g11')] <- 'g11'
  treat[which(plots=='G12'|plots=='g12')] <- 'g12'
  
  
  treat
}

###complete treamtent based on plot seperated by '-'
all_treatment <- function(plots,site='DF'){ #applies WS treatments
  treat <- plots
  if(site=='DF'){
    treat[which(plots=='S03'|plots=='S3'|plots=='s03'|plots=='s3')]<- 'zero-chamber-shade'
    treat[which(plots=='S04'|plots=='S4'|plots=='s04'|plots=='s4')]<- 'zero-chamber-shade'
    treat[which(plots=='S08'|plots=='S8'|plots=='s08'|plots=='s8')] <- 'zero-chamber-shade'
    treat[which(plots=='S02'|plots=='S2'|plots=='s02'|plots=='s2')] <- 'hot-chamber-shade'
    treat[which(plots=='S05'|plots=='S5'|plots=='s05'|plots=='s5')] <- 'hot-chamber-shade'
    treat[which(plots=='S09'|plots=='S9'|plots=='s09'|plots=='s9')] <- 'hot-chamber-shade'
    treat[which(plots=='S01'|plots=='S1'|plots=='s01'|plots=='s1')] <- 'hot-chamber-shade'
    treat[which(plots=='S06'|plots=='S6'|plots=='s06'|plots=='s6')] <- 'hot-chamber-shade'
    treat[which(plots=='S07'|plots=='S7'|plots=='s07'|plots=='s7')] <- 'hot-chamber-shade'
    treat[which(plots=='S10'|plots=='s10')] <- 'zero-nochamber-shade'
    treat[which(plots=='S11'|plots=='s11')] <- 'zero-nochamber-shade'
    treat[which(plots=='S12'|plots=='s12')] <- 'zero-nochamber-shade'
    
    treat[which(plots=='G01'|plots=='G1'|plots=='g01'|plots=='g1')] <- 'zero-chamber-gap'
    treat[which(plots=='G04'|plots=='G4'|plots=='g04'|plots=='g4')] <- 'zero-chamber-gap'
    treat[which(plots=='G07'|plots=='G7'|plots=='g07'|plots=='g7')] <- 'zero-chamber-gap'
    treat[which(plots=='G03'|plots=='G3'|plots=='g03'|plots=='g3')] <- 'hot-chamber-gap'
    treat[which(plots=='G05'|plots=='G5'|plots=='g05'|plots=='g5')] <- 'hot-chamber-gap'
    treat[which(plots=='G09'|plots=='G9'|plots=='g09'|plots=='g9')] <- 'hot-chamber-gap'
    treat[which(plots=='G02'|plots=='G2'|plots=='g02'|plots=='g2')] <- 'hot-chamber-gap'
    treat[which(plots=='G06'|plots=='G6'|plots=='g06'|plots=='g6')] <- 'hot-chamber-gap'
    treat[which(plots=='G08'|plots=='G8'|plots=='g08'|plots=='g8')] <- 'hot-chamber-gap'
    treat[which(plots=='G10'|plots=='g10')] <- 'zero-nochamber-gap'
    treat[which(plots=='G11'|plots=='g11')] <- 'zero-nochamber-gap'
    treat[which(plots=='G12'|plots=='g12')] <- 'zero-nochamber-gap'
  }
  #HARVARD FOREST
  #G02,G04,G07,S02,S05,S07 Ambient
  #G01,G06,G08,S03,S06,S08 3 degrees
  #G03,G05,G09,S01,S04,S09 5 degrees
  #G10,G11,G12,S10,S11,S12 Control
  if(site=='HF'){
    treat[which(plots=='S02'|plots=='S2'|plots=='s02'|plots=='s2')] <- 'zero-chamber-shade'
    treat[which(plots=='S05'|plots=='S5'|plots=='s05'|plots=='s5')] <- 'zero-chamber-shade'
    treat[which(plots=='S07'|plots=='S7'|plots=='s07'|plots=='s7')] <- 'zero-chamber-shade'
    treat[which(plots=='S03'|plots=='S3'|plots=='s03'|plots=='s3')] <- 'hot-chamber-shade'
    treat[which(plots=='S06'|plots=='S6'|plots=='s06'|plots=='s6')] <- 'hot-chamber-shade'
    treat[which(plots=='S08'|plots=='S8'|plots=='s08'|plots=='s8')] <- 'hot-chamber-shade'
    treat[which(plots=='S01'|plots=='S1'|plots=='s01'|plots=='s1')] <- 'hot-chamber-shade'
    treat[which(plots=='S04'|plots=='S4'|plots=='s04'|plots=='s4')] <- 'hot-chamber-shade'
    treat[which(plots=='S09'|plots=='S9'|plots=='s09'|plots=='s9')] <- 'hot-chamber-shade'
    treat[which(plots=='S10'|plots=='s10')] <- 'zero-nochamber-shade'
    treat[which(plots=='S11'|plots=='s11')] <- 'zero-nochamber-shade'
    treat[which(plots=='S12'|plots=='s12')] <- 'zero-nochamber-shade'
    
    treat[which(plots=='G02'|plots=='G2'|plots=='g02'|plots=='g2')] <- 'zero-chamber-gap'
    treat[which(plots=='G04'|plots=='G4'|plots=='g04'|plots=='g4')] <- 'zero-chamber-gap'
    treat[which(plots=='G07'|plots=='G7'|plots=='g07'|plots=='g7')] <- 'zero-chamber-gap'
    treat[which(plots=='G03'|plots=='G3'|plots=='g03'|plots=='g3')] <- 'hot-chamber-gap'
    treat[which(plots=='G05'|plots=='G5'|plots=='g05'|plots=='g5')] <- 'hot-chamber-gap'
    treat[which(plots=='G09'|plots=='G9'|plots=='g09'|plots=='g9')] <- 'hot-chamber-gap'
    treat[which(plots=='G01'|plots=='G1'|plots=='g01'|plots=='g1')] <- 'hot-chamber-gap'
    treat[which(plots=='G06'|plots=='G6'|plots=='g06'|plots=='g6')] <- 'hot-chamber-gap'
    treat[which(plots=='G08'|plots=='G8'|plots=='g08'|plots=='g8')] <- 'hot-chamber-gap'
    treat[which(plots=='G10'|plots=='g10')] <- 'zero-nochamber-gap'
    treat[which(plots=='G11'|plots=='g11')] <- 'zero-nochamber-gap'
    treat[which(plots=='G12'|plots=='g12')] <- 'zero-nochamber-gap'
  }
  
  treat
}


#
#####rolling mean of previous 'roll' numbers if roll is 6 y[7] = mean(y[1:6])
#

rollprev <- function(y,roll=6,include=F){ #roll mean of preceeding values
  ybar <- rep(NA,length(y))
  
  ii <- 1
  if(include){ii <- 0}
  
  for( i in (roll+1):length(y)){
    
    ybar[i] <- mean(y[(i-(roll)):(i-ii)],na.rm=T)
    
  }
  for(i in 2:roll){
    ybar[i] <- mean(y[1:(i-ii)],na.rm=T)
  }
  ybar[1] <- y[1]
  return(ybar)
}


#####rolling mean ofaround numbers if roll is 6 y[7] = mean(y[4:10])
rollmean <- function(y,roll=7){ #roll mean centered on value roll must be odd
  ybar <- rep(NA,length(y))
    
  for( i in 1:length(y)){
    whr <- (i - (roll-1)/2):(i+(roll-1)/2)
    whr <- whr[whr >= 1]
    whr <- whr[whr <= length(y)]
    
    ybar[i] <- mean(y[whr],na.rm=T)
    
  }
   return(ybar)
}

rollquant <- function(y,prob = c(.025,.975),roll=7){ #roll mean centered on value roll must be odd
  ybar <- matrix(NA,length(y),length(prob))
  
  for( i in 1:length(y)){
    whr <- (i - (roll-1)/2):(i+(roll-1)/2)
    whr <- whr[whr >= 1]
    whr <- whr[whr <= length(y)]
    
    ybar[i,] <- quantile(y[whr],probs=prob,na.rm=T)
    
  }
  return(ybar)
}

#####rolling mean ofaround numbers if roll is 6 y[7] = mean(y[4:10])
rollARmean <- function(y,roll=7){ #roll mean centered on value roll must be odd
  ybar <- rep(NA,length(y))
  cof = seq(2,(roll-1/2),by=2)
  cof1 = 1/c(rev(cof),1,cof)
  for( i in 1:length(y)){
    cof=cof1
    whr <- (i - (roll-1)/2):(i+(roll-1)/2)
    cof = cof[which(whr >= 1)]
    whr <- whr[whr >= 1]
    
    cof = cof[which(whr <= length(y))]
    whr <- whr[whr <= length(y)]
    
    ybar[i] <- sum(y[whr]*cof)/sum(cof)
    
  }
  return(ybar)
}

#####rolling mean ofaround numbers if roll is 6 y[7] = mean(y[4:10])
previousARmean <- function(y,roll=4){ #roll mean centered on value roll must be odd
  ybar <- y
  cof = 1/seq(roll,1,by=-1)
  for( i in (roll):length(y)){
    whr <- (i-roll):(i-1)+1

    ybar[i] <- sum(y[i-whr + 1]*rev(cof))/sum(cof)
    
  }
  return(ybar)
}


smTOsp <- function(soilM, soil='tarrus'){ #soil moist to soil potential mckFunctions
  
  if(soil == 'tarrus'){TaD <- cbind(c(46.8,42.7,40.0,37.7,35.0,33.2,30.5,28.2,26.8,25.0,23.6,21.8,20.9,18.2,17.3,16.0,15.0,14.0,13.0,12.0,11.3),
               c(0.09,0.16,0.20,0.24,0.28,0.31,0.39,0.51,0.62,0.80,0.98,1.32,1.54,2.58,3.11,4.12,5.23,6.73,8.85,11.87,15))}
  if(soil == 'DF'){TaD <- cbind(c(50.9,45.5,39.5,35,30.5,25,20,18),
               c(.06,.16,.26,.39,.83,2.48,8.48,15))}
  #if(soil=='DF_Gap'){TaD <- cbind(c(49.4,45.9,40,35,30,25,20.5,17.9),
  #             c(.07,.13,.25,.36,.86,2.42,7.05,15))}
  if(soil=='HF'){TaD <- cbind(c(49.1,44.1,39.1,34.1,29.1,24.1,19.1,14.1,10.9,9.4),
               c(0,.06,.12,.17,.23,.29,.44,2.02,7.23,15))}
  y <- TaD[,2]
  x <- TaD[,1]/100

  fit <- loess(y~x)

  if(is.numeric(dim(soilM))){
    y <- matrix(NA,dim(soilM)[1],dim(soilM)[2])
    for( i in 1:dim(soilM)[2]){
      y[,i] <- predict(fit,soilM[,i])
    }
  }else{y <-  predict(fit,soilM)}
  
  y[which(soilM < .25 & is.na(y))] <- 15
  y[which(soilM > .4 & is.na(y))]   <- .01
  y[which(y < .01)]   <- .01

  y
}

#DFnames <- c('Sand','Silt','Clay')
#Tarrus0_2 <- c(
#Tarrus2_7 <- c(
#Tarrus7_33 <- c(
  
#DFSha0_10 <- c(13.92,64.04,22.04)
#DFSha10_20 <- c(14.5,52.13,33.38)
#DFShaAvg <- c(14.21,58.085,27.71) #3
#TaD <- cbind(c(50.9,45.5,39.5,35,30.5,25,20,18),
#             c(.06,.16,.26,.39,.83,2.48,8.48,15))

#DFGap0_10 <- c(19.92,60.04,20.04) #2.4
#DFGap10_20 <- c(9.92,62.04,28.04)
#DFShaAvg <- c(14.92,61.04,24.04)
#TaD <- cbind(c(49.4,45.9,40,35,30,25,20.5,17.9),
#             c(.07,.13,.25,.36,.86,2.42,7.05,15))

#HFGuess <- c(65,25,10)
#TaD <- cbind(c(49.1,44.1,39.1,34.1,29.1,24.1,19.1,14.1,10.9,9.4),
#             c(0,.06,.12,.17,.23,.29,.44,2.02,7.23,15))

#####


#####
readENV <- function(where = computer,daily=FALSE){  #readENV data warming site from DropB
  
  if(where == 'duke')
  {
    path = 'C:/Users/mck14/Dropbox/warming/ENVDATA/processedData/'
  }else{
    path =  '/home/mckwit/Dropbox/warming/ENVDATA/processedData/' 
  }
  
  if(site == 'Duke'){
    SITE = 'Duke'
  }else{
    SITE = 'Harvard'
  } 
  
  if(daily){

    SM <- read.csv(paste(path,'daily-',SITE,'-SM.csv',sep=''))
    ST <- read.csv(paste(path,'daily-',SITE,'-ST.csv',sep=''))
    AT <- read.csv(paste(path,'daily-',SITE,'-AT.csv',sep=''))
    Q  <- read.csv(paste(path,'daily-',SITE,'-Q.csv',sep=''))
    RH <- read.csv(paste(path,'daily-',SITE,'-Rh.csv',sep=''))

  }else{
    
    SM <- read.csv(paste(path,'gapfill-data-',SITE,'-SM.csv',sep=''))
    ST <- read.csv(paste(path,'gapfill-data-',SITE,'-ST.csv',sep=''))
    AT <- read.csv(paste(path,'gapfill-data-',SITE,'-AT.csv',sep=''))
    Q  <- read.csv(paste(path,'gapfill-data-',SITE,'-Q.csv',sep=''))
    RH <- read.csv(paste(path,'gapfill-data-',SITE,'-Rh.csv',sep=''))
  
  }
  
  list(SM= SM,ST=ST,AT=AT,Q=Q,RH=RH)
}

################################################

###
#Fuzzy MATCH 
fuzWHICH <- function(DFwant, DFref){  #Fuzzy match
  if(length(DFwant) == 1){
    whichTIME <- which((abs(DFwant - DFref)) == min(abs(DFref - DFwant))) 
    out = whichTIME[1]
  }else{
    out = numeric()
    for(i in 1:length(DFwant)){
    out <- c(out,which((abs(DFwant[i] - DFref)) == min(abs(DFref - DFwant[i])))[1]) 
    out
    }
  }
  return(out)
}
##
###
nearWHICH <- function(data,number){ #find closest point
  ind <- which(abs(data-number)==min(abs(data-number),na.rm=T))
  if(length(ind)==0){ind <- NA}
  ind[1]
}

####Function to take positive area under curve######

Pos.A.U.Curve <- function(datetime,magnitude){ #positive area under curve
	magnitude[which(magnitude < 0)] <- 0
	seconds <- (as.duration(interval(time[1], time)))
	sum(diff(seconds)*rollmean(magnitude,2))
}

#####Function converts Type T Thermocouple between mVolt and Celsius and vis versa
TempC.mV <- function(input,TempC.or.mV){ #convert Type T between mVolt and C
  
  mV <- c(0,0.039,0.078,0.117,0.156,0.195,0.234,0.273,0.312,0.352,0.391,0.431,0.47,0.51,0.549,0.589,0.629
          ,0.669,0.709,0.749,0.79,0.83,0.87,0.911,0.951,0.992,1.033,1.074,1.114,1.155,1.196,1.238,1.279
          ,1.32,1.362,1.403,1.445,1.486,1.528,1.57,1.612,1.654,1.696,1.738,1.78,1.823,1.865,1.908,1.95
          ,1.993,2.036,2.079,2.122,2.165,2.208,2.251,2.294,2.338,2.381,2.425,2.468,2.512,2.556,2.6,2.643
          ,2.687,2.732,2.776,2.82,2.864,2.909,2.953,2.998,3.043,3.087,3.132,3.177,3.222,3.267,3.312,3.358
          ,3.403,3.448,3.494,3.539,3.585,3.631,3.677,3.722,3.768,3.814,3.86,3.907,3.953,3.999,4.046,4.092
          ,4.138,4.185,4.232,4.279,4.325,4.372,4.419,4.466,4.513,4.561,4.608,4.655,4.702,4.75)
  
  C <- 0:110
  TypeT <- cbind(C,mV)
	
  if(TempC.or.mV == 'mV'){
		TmV <- smooth.spline(TypeT[,'C']~TypeT[,'mV'])
	}
	if(TempC.or.mV == 'TempC'){
		TmV <- smooth.spline(TypeT[,'mV']~TypeT[,'C'])
	}
	predict(TmV,input)$y
 }

#####Function changes temperature from campbell logger into voltage(mV) and rectifies

Logger.Temp.Convert <- function(panelT,thermoT){ #campbell T into mV and rectify 
	thermoV <- TempC.mV(thermoT,'TempC') - TempC.mV(panelT,'TempC')
	if(mean(thermoV,na.rm=T) < 0 ){ thermoV <- -thermoV}
	thermoV
}

##### Day frac to minutes or hours

fracTOtime <- function(DayFraction, roundit = FALSE,hourORmin = 'hour' ){ #Day frac to minutes or hours
	if(hourORmin == 'min'){
		if(roundit == T){Min <- round(DayFraction*1440)}
		if(roundit == F){Min <- DayFraction*1440}
	}
	if(hourORmin == 'hour'){
		if(roundit == T){Min <- round(DayFraction*24)}
		if(roundit == F){Min <- DayFraction*24}
	}
		Min
} 

timeToFrac <- function(hhmmss,sep = ":"){
  
  tmp = matrix(unlist(strsplit(hhmmss,split = sep)),ncol=3,byrow=T)
  as.numeric(tmp[,1])/24 + as.numeric(tmp[,2])/1440 + as.numeric(tmp[,1])/86400

}

### Make a date from a fraction

makeDATE <- function(Years,Months,Days,Hours = 0,Mins = 0,Secs = 0){ ### Make a date from a fraction
	TheDate <- as.Date(paste(Years,Months,Days,sep='-'))
	TheDate + hours(Hours) + minutes(Mins) + seconds(Secs)
	
	}

##CR1000 timestamp to date or julian date base can be changed
TIMESTAMPtoDATE <- function(times,base=as.Date('2008-12-31'),julian=TRUE){ ##CR1000 timestamp to date or jd
  require(lubridate)
  zero <- ymd_hms('1999-12-31 01:01:01')
  dtimes <-  rep(zero,length(times))
  dtimes[!is.na(times)] <- mdy_hm(times[!is.na(times)])
  dtimes[is.na(times)] <-NA
  if(julian){
    jd <- julian(dtimes,origin=base)
    jd
  }else{
    dtimes
  }
}

######
#strptime('11/1/2012 0:00', format="%Y/%m/%d %H:%M:%S", tz = Sys.timezone())

### F to celcius

tempCONVERT<- function(inTEMP,InputUnit='C'){ #F to Celcius
if(InputUnit == 'F'){ outtemp <- (inTEMP - 32)/1.8}
if(InputUnit == 'C'){ outtemp <- (inTEMP * 1.8) + 32}
outtemp
}

##convert rh to vpd output in KPa

rhTOvpd<- function (rh,T){ #rH to vpd KPa
  
  if(sum(dim(rh))==0){
    es<- 6.112 * exp((17.67*T)/(T+243.5))
  }else{
    es <- matrix(NA,dim(rh)[1],dim(rh)[2])
    for(j in 1:dim(rh)[2]){
    es[,j] <- as.vector(6.112 * exp((17.67*T[,j])/(T[,j]+243.5))) #Saturated Vapor Pressure using Celsius
    }
  }
  
e <- (rh * es)/100
VPD = es - e

VPD/10
}	



rhTOvpd2<- function (rh,T){  #celsius Rh  %  VPD KPa
  
  if(sum(dim(rh))==0){
    SVP<- 610.7 * 10^(7.5*T/(237.3+T))
  }else{
    SVP <- matrix(NA,dim(rh)[1],dim(rh)[2])
    for(j in 1:dim(rh)[2]){
     SVP[,j] <- as.vector(610.7 * 10^(7.5*T/(237.3+T))) #Saturated Vapor Pressure using Celsius
    }
  }
  
  VPD <- ((100-rh)/100)*SVP
  VPD = VPD/1000
  
  VPD
}  


#####
#### this function returns the index of max ormin 'which.max' number needs to be middle of odd roll num
localMS <- function(inDATA,period,MorM ='MAX'){ #find local min or max
	require(zoo)
	if(period%%2 == 0)stop('Period must be odd')
	useDATA <- c(inDATA,rep(mean(inDATA,na.rm=T),period))
	if(MorM == 'MAX'){xz <- as.zoo(useDATA)
		rxz <- rollapply(xz, period, function(useDATA) which.max(useDATA)==quantile(1:period,.5))
		}
	if(MorM == 'MIN'){xz <- as.zoo(useDATA)
		rxz <- rollapply(xz, period, function(useDATA) which.min(useDATA)==quantile(1:period,.5))
		}
	
	index(rxz)[coredata(rxz)] 
}


###
degreeHOURS <- function(JD2009,month,day,temps,lowThresh = 4,jdSTART = c(365,731,1096),jdLENGTH = 150){ #calc degree hours
		
	require(zoo)
	forFIT <- na.fill(temps,18)
	
	for(i in 1:12){
		plotMEAN <- apply(temps[which(month==i),],2,mean,na.rm=T)
	
		for(j in 1:length(temps[1,])){
			forFIT[which(is.na(temps[,j]) & month==i),j] <- plotMEAN[j]
		}
		
	}
			
	for(i in 1:length(temps[1,])){
		fit <- smooth.spline(JD2009,forFIT[,i])
		temps[is.na(temps[,i]),i] <- predict(fit,JD2009[is.na(temps[,i])])$y
	}
	
	dH <- as.matrix(temps - lowThresh)
	dH[which(dH < 0,arr.ind=T)] <- 0
	
	forFIT[forFIT < 9999999] <- 0
	
	for( i in  1:3){
			START <- which(JD2009==jdSTART[i])[1]
			FINISH <- tail(which(JD2009==(jdSTART[i]+jdLENGTH)),1)
		forFIT[START:FINISH,] <- apply(dH[START:FINISH,],2,cumsum)
	
	}
	
	forFIT
	
}

chillHOURS <- function(JD2009,month,day,temps,base = 7,jdSTART = c(304,669,1035),jdLENGTH = 150){ #calc chill hours
		
	require(zoo)
	forFIT <- na.fill(temps,18)
	
	for(i in 1:12){
		plotMEAN <- apply(temps[which(month==i),],2,mean,na.rm=T)
	
		for(j in 1:length(temps[1,])){
			forFIT[which(is.na(temps[,j]) & month==i),j] <- plotMEAN[j]
		}
		
	}
			
	for(i in 1:length(temps[1,])){
		fit <- smooth.spline(JD2009,forFIT[,i])
		temps[is.na(temps[,i]),i] <- predict(fit,JD2009[is.na(temps[,i])])$y
	}
	
	cH <- as.matrix(base - temp )
	cH[which(temp > base)] <- 0
	cH[which(temp < 0)] <- 0
	
	forFIT[forFIT < 9999999] <- 0
	
	for( i in  1:3){
			START <- which(JD2009==jdSTART[i])[1]
			FINISH <- tail(which(JD2009==(jdSTART[i]+jdLENGTH)),1)
		forFIT[START:FINISH,] <- apply(dH[START:FINISH,],2,cumsum)
	
	}
	
	forFIT
	
}

################

findBREAKSrepeat <- function(REPdata,value){ #finds which location has last occurance of value
	sequen <- which(floor(REPdata)  <= value)
	c(sequen[which(diff(sequen) > 1)], tail(sequen,1))
	}
	
	cleanMEAN <- function(data.set){
	for(i in which(is.na(data.set))){
	
		if(i-1 <= 0){data.set[i] <- median(data.set[i:(i+7)],na.rm=T);next}
		if(i+1 >= length(data.set)){data.set[i] <- median(data.set[i:(i-7)],na.rm=T);next}
		data.set[i] <- median(data.set[(i-7):(i+7)],na.rm=T)
	
	}
	data.set
}

###############

rowBINDER <- function(matnow,row2add,rowName){ #binds rows (may be obsolete)

		if(length(matnow) == 0){
			matnow <- row2add
			if(!is.matrix(row2add)){
				matnow <- matrix(matnow,nrow=1)
			if(!is.null(names(row2add)))colnames(matnow) <- names(row2add)
			}
			rownames(matnow) <- rowName
			return(matnow)
		}
	addSTART <- dim(matnow)[1]+1
	matnow <- rbind(matnow,row2add)
	rownames(matnow)[addSTART:nrow(matnow)] <- rowName
	matnow
	}

#################


# this function takes a warming site env matrix 
#grabs the daily mean and sums across an interval 6 is one week
cumLAGenv <- function(tempMAT,interval=6){  #weekly mean temp WS
  
  days <- unique(tempMAT[,'JD2009'])
  dayMEAN <- matrix(NA,length(days),dim(tempMAT[,-c(1,2,3,4,5,7)])[2])
  for(i in 1:length(days)){
    ind <- which(tempMAT[,'JD2009']==days[i])
    dayMEAN[i,] <- apply(tempMAT[ind,-c(1,2,3,4,5,7)],2,mean,na.rm=T)
  }
  
  dayTOT <- matrix(0,dim(dayMEAN)[1],dim(dayMEAN)[2])
  
  for(i in (1+interval):length(days)){
    dayTOT[i,-1] <- apply(dayMEAN[(i-interval):i,-1],2,mean,na.rm=T)
  }
  dayTOT[,1] <- dayMEAN[,1]
  colnames(dayTOT) <- colnames(dayMEAN) <- colnames(tempMAT[,-c(1,2,3,4,5,7)])
  list(dailyMEAN = dayMEAN,intervalSUM <- dayTOT,COLnames <- colnames(dayTOT) )
}

######## turns in z-scores, use 2 stds to compare between covariates and factors

CENSTAND <- function(xXx,numSTDS = 1){ #z-scores mckFUNCTIONS
  dimxXx = length(dim(rep(1,10)))
  if(dimxXx == 0){ cen <- xXx - mean(xXx,na.rm=T) 
    if(numSTDS != 0){ cen <- cen/(numSTDS*sd(xXx,na.rm=T)) }
  }
  if(dimxXx > 1){
    xmean <- apply(xXx,2,mean,na.rm=T) 
    xsd <- numSTDS*apply(xXx,2,sd,na.rm=T)
    cen <- (xXx - matrix(xmean,nrow(xXx),ncol(xXx),byrow=T) )
    if(numSTDS != 0){
      cen <- cen/(matrix(xsd,nrow(xXx),ncol(xXx),byrow=T) )
    }
  }
  cen
}  


########factor to a dummy variable

facTOdummy <- function(fact, firstZERO = T){ #makes dummy variables
  if(firstZERO == T){
  is <- sort(unique(fact))
  numIS <- length(is)
  dumbMAT <- matrix(0,length(fact),(numIS-1))
  
  for(i in 2:numIS){
    dumbMAT[which(fact == is[i]),i-1] <- 1
  }
  colnames(dumbMAT) <- is[-1]
  }
  
  if(firstZERO == F){
    is <- sort(unique(fact))
    numIS <- length(is)
    dumbMAT <- matrix(0,length(fact),(numIS))
    
    for(i in 1:numIS){
      dumbMAT[which(fact == is[i]),i] <- 1
    }
    colnames(dumbMAT) <- is
  }
  
  dumbMAT
}

#################
###iteratively subtracts random effects from y
ranFy <- function(yy,bb,rr){  #iteratively subtracts random effects from y                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  
  num <- dim(bb)[2]
  start <- yy
  for(i in 1:num){
    start <- start - bb[,i]*rr[,i]
  }
  start
}

############### apply tapply
tapply.matrix <- function(mat,margin,index,funct,...){ #tapply on matrix mckFUNCTIONS
  apply(mat,margin,function(x){ tapply(x,index,funct,na.rm=T)})  
}


#####find mode
my.mode <- function(x) { #find mode
  d <- density(x)
  d$x[which.max(d$y)]
}

######check design matrix, name of intercept column
check.design <- function(x,intName='intercept'){  # design matrix diagnostics
  
  p <- ncol(x)
  
  if(is.null(colnames(x))){
    colnames(x) <- paste('x',c(1:p),sep='_')
    xrange      <- apply(x,2,range)
    wi          <- which(xrange[1,] == 1 & xrange[2,] == 1)
    if(length(wi) > 0)colnames(x)[wi] <- 'intercept'
  }
  
  wi <- which(colnames(x) == 'intercept')
  
  xname <- colnames(x)
  
  VIF <- rep(NA,p)
  names(VIF) <- xname
  for(k in 1:p){
    
    if(xname[k] == intName)next
    
    notk <- xname[xname != xname[k]]
    
    ykk <- x[,xname[k]]
    xkk <- x[,notk]
    
    tkk <- summary(lm(ykk ~ xkk))$adj.r.squared
    VIF[k] <- 1/(1 - tkk)
  }
  
  VIF <- VIF[-wi] 
  
  corx <- cor(x[,-wi])
  rankx <- qr(x)$rank
  
  list(VIF = round(VIF,1), correlation = round(corx,2), rank = rankx, p = p)
}

print.pars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,  #print gibbs and plot
                        sigOnly = F,burnin=1,xlimits = NULL){  
  
  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
      btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }
    
    wq   <- apply(btmp,2,quantile,c(.025,.975))  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    if(length(wq) > 0){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
  }
  
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
    if(burnin > (nrow(xgb) + 100))stop("burnin too large")
    xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)
  
  out <- t(rbind(apply(xgb,2,mean),apply(xgb,2,quantile,c(.025,.975))))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('mean','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('mean','0.025','0.975','true value')
  }
  
  armat <- matrix(0,nc,10)  #for AR model
  
  # for(j in 1:nc)armat[j,] <- ar(xgb[,j],aic=F,order.max = 10)$ar[1:10]
  
  # if(!is.null(colnames(xgb)))rownames(armat) <- colnames(xgb)
  #  colnames(armat) <- paste('AR',c(1:10),sep='-')
  
  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf),mar=c(1.5,2,1.5,1))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc),mar=c(1.5,2,1.5,1))
  
  
  if(CPLOT){
    for(j in 1:nc){
      plot(xgb[,j],type='l',axes=FALSE)
      abline(h=out[j,],lty=2)
      if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
      abline(h=0,col='red',lwd=2)
      title(colnames(xgb)[j])
      axis(side=2,at=quantile(xgb[,j],probs=c(.5)),
           labels=round(quantile(xgb[,j],probs=c(.5)),3))
      axis(side=1)
    }
  }
  xlims <- xlimits
  if(DPLOT){
    for(j in 1:nc){
      xj <- density(xgb[,j])
      if(is.null(xlimits))xlims <- range(xj$x)
      plot(xj$x,xj$y,type='l',xlim=xlims)
      abline(v=out[j,],lty=2)
      if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
      abline(v=0,col='red',lwd=2)
      title(colnames(xgb)[j])
    }
  }
  list(summary = signif(out,4))
  
}


#
# Helper function to draw a line for each series on existing plot
#
# Series:matrix[individuals,measures], rownames:ids, colnames:ordered names,provide: numeric equivalent in "NameColUsed"
# NameColUsed: Names of columns to include can be 2009,2010 if colnames are X2009,X2010,...
# boldSeries: Name of rows to emphasize grep(names)
# col:main color of lines; lwd:main line weight
# hcol:bold line col; hlwd:bold lind width
#

ManyLinePlots <- function(series,NameColUsed,boldSeries = NA, col="black", hcol="black", lwd=1, hlwd=2, showPoints=FALSE) { #plot many lines on one plot
  startPoint = NameColUsed[1]
  endPoint = rev(NameColUsed)[1]
  columns <- grep(paste(NameColUsed,collapse="|"),colnames(series))
  for (i in 1:length(series[,1])) {
    lines(NameColUsed, series[i,columns], col=col, lwd=lwd)
  }
  
  if (!is.na(boldSeries)) {
    boldIndex <- grep(paste(boldSeries,collapse="|"),rownames(series))
    lines(NameColUsed, series[boldIndex,columns], col=hcol, lwd=hlwd)
  }
  
  if (showPoints) {
    points(NameColUsed, series[boldIndex,columns], col=hcol, pch=4)
  }
} 

#
# Newspaper style many line plot
#

news.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit = "",colUsed = "",embolden = 1){ # Newspaper style many line plot
  if(colUsed == "")colUsed <- minX:maxX
  par(mar=c(4, 4, 3, 1), oma=c(0,0,0,0), xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.6,0), las=1, lend=1)
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), family="Helvetica", cex.axis=0.8, cex.lab=0.8)
  grid(NA, NULL, col="black", lty="dotted", lwd=0.3)
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, col="dark grey", lwd=0.7, hlwd=3, hcol="#244A5D")
  
}

#
# Feltron style many line plot
#

feltron.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit ="",colUsed = "",embolden = 1){ # Feltron style many line plot
  
  if(sum(colUsed == "") == 1)colUsed <- minX:maxX
  par(bg="#36394A", mar=c(5, 4, 3, 2), oma=c(0,0,0,0), xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.8,0.3,0.5), col.lab="white", col.axis="white", col.main="white", font.main=1, cex.main=0.8, cex.axis=0.8, cex.lab=0.8, family="Helvetica", lend=1, tck=0)
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit)))) #asp=1/2, 
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, col="white", lwd=0.35, hcol="#E3DF0C", hlwd=3)
  
}
#
# Feltron no-axis style many line plot
#

felt_noaxis.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit ="",colUsed = "",embolden = 1){ # Feltron no-axis style many line plot
  
  if(colUsed == "")colUsed <- minX:maxX
  par(bg="#36394A", mar=c(5, 4, 3, 2), oma=c(0,0,0,0), xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.8,0.3,0.5), col.lab="white", col.axis="white", col.main="white", font.main=1, cex.main=0.8, cex.axis=0.8, cex.lab=0.8, family="Helvetica", lend=1, tck=0, las=1)
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), xaxt="n", yaxt="n")
  axis(1, tick=FALSE, col.axis="white")
  axis(2, tick=FALSE, col.axis="white")
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, col="white", lwd=0.35, hcol="#E3DF0C", hlwd=3)
  
}

#
# FiveThirtyEight style many line plot
#

five38.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit ="",colUsed = "",embolden = 25){ # FiveThirtyEight style many line plot
  opar <- par()
  if(colUsed == "")colUsed <- minX:maxX
  par(mar=c(3, 4, 3, 2), oma=c(0,0,0,0), bg="#F0F0F0", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
  grid(NULL, NULL, col="#DEDEDE", lty="solid", lwd=0.9)
  axis(1, tick=FALSE, cex.axis=0.9)
  axis(2, tick=FALSE, cex.axis=0.9)
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, col="dark grey", lwd=1, hlwd=3, hcol="#008ED4")
  on.exit(par(opar))
}
#
# The Economist style many line plot
#

# Red corner rectangle
econ.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit ="",colUsed = "",embolden = 25){ # The Economist style many line plot
  opar <- par()
  if(colUsed == "")colUsed <- minX:maxX
  par(xpd=NA, oma=c(0,0,0,0), mar=c(0,0,0,0), bg="#DCE6EC", xpd=FALSE, xaxs="i", yaxs="i", lend=1)
  plot(0, 0, type = "n", bty = "n", xaxt="n", yaxt="n", xlim=c(0,100), ylim=c(0,100))
  rect(0,100,2,94, col="red", border=NA)
  
  # Actual chart
  par(mar=c(4, 3, 3, 2), oma=c(0,0,0,0), xpd=FALSE, xaxs="r", yaxs="i", mgp=c(1.8,.2,0), cex.axis=0.7, cex.lab=0.7, col.lab="black", col.axis="black", col.main="black", tck=0.02, yaxp=c(minVal, maxVal, 2), new=TRUE)
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), family="Helvetica", asp=1/2)
  grid(NA, NULL, col="white", lty="solid", lwd=1.5)
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, lwd=1.25, hlwd=2.5, col="#33A5A2", hcol="#244A5D")
  par(opar)
}

#
# Tukey style many line plot
#
tukey.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit ="",colUsed = "",embolden = 25){ # Tukey style many line plot
  opar <- par()
  if(colUsed == "")colUsed <- minX:maxX
  par(las=1, tck=0.02, mgp=c(2.8,0.3,0.5), cex.lab=0.85, cex.axis=0.8, cex.main=0.9)
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), asp=1/2)
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, col="#cccccc", hlwd=1.2, showPoints=TRUE)
  par(opar)
}

#
# Bright on Dark style many line plot
#
BonD.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit ="",colUsed = "",embolden = 25){ # Bright on Dark style many line plot
  opar <- par()
  if(colUsed == "")colUsed <- minX:maxX
  par(bg="black", las=1, tck=0, mgp=c(2.8,0.3,0), cex.lab=0.85, cex.axis=0.8, cex.main=0.9, col.axis="white", col.main="white", col.lab="white")
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), asp=1/2)
  grid(NULL, NULL, lty="solid", col="white", lwd=0.5)
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, col="#f30baa", lwd=1.2, hcol="green", hlwd=3)
  par(opar)
}


#non-rectangular 
n.rect.h <- function(I,Amax,Rd,Q,Theta){ #non-rectangular hyperbola
  (1/(2*Theta))*(Q*I+(Amax+Rd)-sqrt(((Q*I+(Amax+Rd))^2)-4*Theta*Q*I*(Amax+Rd)))+Rd
}

e.rise.m <- function(I,Amax,Rd,Q){ #exponetial rise to max
  (Amax-Rd)* (1-exp(-Q*I)) + Rd
}

mm.curve <- function(I,Amax,Rd,K){ #m-m saturation
  (Amax+Rd)*I/ (K + I) - Rd
}



#split date(month day year) into three columns produces matrix [ndates X out.order]
split.date <- function(date,sep="/",in.order = c('m','d','y'),out.order = c('m','d','y')){#split date into three columns produces matrix
  
  m <- which(out.order == 'm')
  d <- which(out.order == 'd')
  y <- which(out.order == 'y')
  
  whr <- match(out.order,in.order)
  
  out <- matrix(a.n(unlist(strsplit(a.c(date),split=sep))),nrow=length(date),ncol=3,byrow=T)[,whr]
  if(length(out) == 3){
    names(out) <- out.order
  }else{  
    colnames(out) <- out.order
  }
  out
  
}

#Julian Date using chron, input is date as character
julian.date <- function(date,sep="/",origin = c(month = 11, day = 9, year = 1980)){ #Julian Date using chron
  library(chron)
  out.order <- c('m','d','y')
  hold.date <- split.date(date)
  if(length(date) > 1){
    julian(hold.date[,'m'],hold.date[,'d'],hold.date[,'y'],origin = origin) #turn into a JD
  }else{
    julian(hold.date['m'],hold.date['d'],hold.date['y'],origin = origin) #turn into a JD
  }
}



###for WarmingScript.r

add.plot.info <- function(datatable){###for WarmingScript.r
  num <- dim(datatable)[2]
  datatable[,num+1] <- 'NA'
  colnames(datatable)[num+1] <- 'Temp'
  datatable[,num+2] <- 'NA'
  colnames(datatable)[num+2] <- 'Light'
  datatable[,num+3] <- 'NA'
  colnames(datatable)[num+3] <- 'ID'
  #datatable[,num+4] <- as.character(strptime( paste(datatable[,'Date'],datatable[,'HHMMSS'],sep=" "),format="%m/%d/%Y %T", tz="EST5EDT")) #turn into a date
  #colnames(datatable)[num+4] <- 'DTime'
  
  datetime <- chron(a.c(data.table[,'Date']),a.c(data.table[,'HHMMSS'])) #turn into a date
  
  #datatable[,'FTime'] <- as.numeric(julian(datetime,origin = as.POSIXct("2012-01-01",tz='UTC'),digits=9))-(6/24)
  
  Control <- c('g10','g11','g12','s10','s11','s12')
  Ambient <- c('g01','g04','g07','s03','s04','s08')
  Three <- c('g03','g05','g09','s02','s05','s09')
  Five <- c('g02','g06','g08','s01','s06','s07')
  
  Plot <- rep('NA',length(datatable[,'Plot']))
  
  
  Plot[which(datatable[,'Plot']  == 'g01')] <- 'g01'
  Plot[which(datatable[,'Plot']  == 'g02')] <- 'g02'
  Plot[which(datatable[,'Plot']  == 'g03')] <- 'g03'
  Plot[which(datatable[,'Plot']  == 'g04')] <- 'g04'
  Plot[which(datatable[,'Plot']  == 'g05')] <- 'g05'
  Plot[which(datatable[,'Plot']  == 'g06')] <- 'g06'
  Plot[which(datatable[,'Plot']  == 'g07')] <- 'g07'
  Plot[which(datatable[,'Plot']  == 'g08')] <- 'g08'
  Plot[which(datatable[,'Plot']  == 'g09')] <- 'g09'
  Plot[which(datatable[,'Plot']  == 's01')] <- 's01'
  Plot[which(datatable[,'Plot']  == 's02')] <- 's02'
  Plot[which(datatable[,'Plot']  == 's03')] <- 's03'
  Plot[which(datatable[,'Plot']  == 's04')] <- 's04'
  Plot[which(datatable[,'Plot']  == 's05')] <- 's05'
  Plot[which(datatable[,'Plot']  == 's06')] <- 's06'
  Plot[which(datatable[,'Plot']  == 's07')] <- 's07'
  Plot[which(datatable[,'Plot']  == 's08')] <- 's08'
  Plot[which(datatable[,'Plot']  == 's09')] <- 's09'
  
  Plot[which(datatable[,'Plot']  == 'g1')] <- 'g01'
  Plot[which(datatable[,'Plot']  == 'g2')] <- 'g02'
  Plot[which(datatable[,'Plot']  == 'g3')] <- 'g03'
  Plot[which(datatable[,'Plot']  == 'g4')] <- 'g04'
  Plot[which(datatable[,'Plot']  == 'g5')] <- 'g05'
  Plot[which(datatable[,'Plot']  == 'g6')] <- 'g06'
  Plot[which(datatable[,'Plot']  == 'g7')] <- 'g07'
  Plot[which(datatable[,'Plot']  == 'g8')] <- 'g08'
  Plot[which(datatable[,'Plot']  == 'g9')] <- 'g09'
  Plot[which(datatable[,'Plot']  == 's1')] <- 's01'
  Plot[which(datatable[,'Plot']  == 's2')] <- 's02'
  Plot[which(datatable[,'Plot']  == 's3')] <- 's03'
  Plot[which(datatable[,'Plot']  == 's4')] <- 's04'
  Plot[which(datatable[,'Plot']  == 's5')] <- 's05'
  Plot[which(datatable[,'Plot']  == 's6')] <- 's06'
  Plot[which(datatable[,'Plot']  == 's7')] <- 's07'
  Plot[which(datatable[,'Plot']  == 's8')] <- 's08'
  Plot[which(datatable[,'Plot']  == 's9')] <- 's09'
  
  Plot[which(datatable[,'Plot']  == 'g10')] <- 'g10'
  Plot[which(datatable[,'Plot']  == 'g11')] <- 'g11'
  Plot[which(datatable[,'Plot']  == 'g12')] <- 'g12'
  Plot[which(datatable[,'Plot']  == 's10')] <- 's10'
  Plot[which(datatable[,'Plot']  == 's11')] <- 's11'
  Plot[which(datatable[,'Plot']  == 's12')] <- 's12'
  
  datatable[,'Plot'] <- Plot
  
  
  datatable[datatable[,'Plot'] %in%  Ambient,'Temp'] <- 'Ambient'
  datatable[datatable[,'Plot'] %in%  Control,'Temp'] <- 'Control'
  datatable[datatable[,'Plot'] %in%  Three,'Temp'] <- 'Three'
  datatable[datatable[,'Plot'] %in%  Five,'Temp'] <- 'Five'
  
  if(sum(grep('s',datatable[,'Plot'])) > 0){ 
    datatable[grep('s',datatable[,'Plot']),'Light'] <- 'Shade'}
  if(sum(grep('g',datatable[,'Plot'])) > 0){ 
    datatable[grep('g',datatable[,'Plot']),'Light'] <- 'Gap'}
  
  datatable[,'ID'] <- paste(datatable[,'Plot'],datatable[,'Temp'],datatable[,'Row'],datatable[,'Column'],datatable[,'species'],datatable[,'Date'],sep="/")
  datatable
}

#random multivariate normal
myrmvnorm <- function (nn, mu, sigma){ #random multivariate normal
  
  sigsvd <- svd(sigma)
  retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  
  retval <- matrix(rnorm(nn * ncol(sigma)), nn) %*% retval
  retval + mu
}


tmvnorm.my <- function(avec,muvec,smat,lo=rep(-Inf,length(avec)),hi=rep(Inf,length(avec)),
         whichSample=c(1:length(avec)),times=1){   
  
  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 
  # whichSample indicates which variables to sample
  
  if(length(lo) == 1)lo <- rep(lo,length(avec))
  if(length(hi) == 1)hi <- rep(hi,length(avec))
  
  for(j in 1:times){
    for(k in whichSample){
      
      tmp <- conditionalMVN(avec[k],muvec[k,],smat,k)
      muk <- tmp$mu
      sgk <- tmp$vr
      
      if(length(muk) == 0)next
      
      avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
    }
  }
  avec
}



histo <- function(data, horiz = F,nbrks = 20,binData = data, dataRange = c(min(binData,na.rm=T),max(binData,na.rm=T)),
                  colours = rgb(160/256,82/256,45/256,.75),b.colours=rgb(160/256,82/256,45/256,1),xlimits = NULL, ylimits = NULL,addin = F,
                  box='n',xaxe = "n",yaxe="n",xlabel="",ylabel="",type='counts'){  #horizontal andv vert histo presence and adj.presence
  #tmp = hist(data,plot=F, breaks = nbrks)
  breaks <- seq(dataRange[1],dataRange[2],length =nbrks)
  bins <- cut(binData,breaks = breaks)
  
  #total counts of each bin 
  normal = tapply(bins,bins,length)
  #normal[is.na(normal)] = 0
  counts = tapply(data,bins,sum,na.rm=T)
  value = tapply(data,bins,mean,na.rm=T)
  percent = (counts/normal)
  adj.percent = percent*value
  
  if(type == 'percent') counts = percent
  if(type == 'adj.percent') counts = adj.percent
  
  counts[is.na(counts)] = 0
  if(horiz == FALSE){
    a=breaks[1:(length(breaks) - 1)]
    b=0
    c=breaks[2:length(breaks)]
    d=counts
    ylimits = c(0,max(d))
    xlimits = c(min(a),max(c))
  }
  if(horiz == TRUE){
    b=breaks[1:(length(breaks) - 1)]
    a=0
    d=breaks[2:length(breaks)]
    c=counts
    xlimits = c(0,max(c))
    ylimits = c(min(b),max(d))
  }
  if(!addin){
    plot(NULL,xlim= xlimits,ylim = ylimits,type='n',bty=box,xaxt = xaxe,yaxt=yaxe,xlab=xlabel,ylab=ylabel) 
  }
  rect(a,b,c,d,border=b.colours,col=colours)
}


colorRange <- function(percent = .5,colours = c('#AF2F03','white','#027A40'),trans = 1){
  #cols = rgb2hsv(col2rgb(colours))
  cols = col2rgb(colours)/256
  ncolours = length(colours)
  whr = percent*(ncol(cols)-1)+1
  whr.lo = floor(whr)
  whr.hi = ceiling(whr)
  return(
      rgb( #hsv
      (ncolours - 1)*(percent - (whr.lo-1)/(ncolours - 1)) * (cols[1,whr.hi] - cols[1,whr.lo]) + cols[1,whr.lo],
      (ncolours - 1)*(percent - (whr.lo-1)/(ncolours - 1)) * (cols[2,whr.hi] - cols[2,whr.lo]) + cols[2,whr.lo],
      (ncolours - 1)*(percent - (whr.lo-1)/(ncolours - 1)) * (cols[3,whr.hi] - cols[3,whr.lo]) + cols[3,whr.lo],
      trans
      )
  )
}

bivarColor = function(v1 = .5,v2 = .5, cols = c("#64acbe","#627f8c","#574249",
                                                "#b0d5df","#ad9ea5","#985356",
                                                "#e8e8e8","#e4acac","#c85a5a"), trans = 0){
 
  red = col2rgb(cols)[1,]/256
  green = col2rgb(cols)[2,]/256
  blue = col2rgb(cols)[3,]/256
  red[red > 1] = 1;green[green > 1] = 1;blue[blue > 1] = 1
  red[red < 0] = 0;green[green < 0] = 0;blue[blue < 0] = 0
  x = rep(seq(0,1,length = sqrt(length(cols))),sqrt(length(cols)))
  y = rep(seq(1,0,length = sqrt(length(cols))),each = sqrt(length(cols)))
  red = predict(lm(red   ~ x + y),data.frame(x = v1,y = v2))
  green = predict(lm(green ~ x + y),data.frame(x = v1,y = v2))
  blue = predict(lm(blue   ~ x + y),data.frame(x = v1,y = v2))
  red[red > 1] = 1;green[green > 1] = 1;blue[blue > 1] = 1
  red[red < 0] = 0;green[green < 0] = 0;blue[blue < 0] = 0
  rgb(red,green,blue,trans)

}

bivarSplineColor = function(v1 = .5,v2 = .5,cols = c("#64acbe","#627f8c","#574249",
                                                     "#b0d5df","#ad9ea5","#985356",
                                                     "#e8e8e8","#e4acac","#c85a5a"), trans = 0){
  require(akima)
  #fuzWHICH() needs in mckFunctions or helper.github
  red = col2rgb(cols)[1,]/256
  green = col2rgb(cols)[2,]/256
  blue = col2rgb(cols)[3,]/256
  x = rep(seq(0,1,length = sqrt(length(cols))),sqrt(length(cols)))
  y = rep(seq(1,0,length = sqrt(length(cols))),each = sqrt(length(cols)))
  surface = interp(x= y,y= x,z = red)
  xs = fuzWHICH(v1,surface$y)
  ys = fuzWHICH(v2,surface$x)
  red = surface$z[cbind(ys,xs)]
  surface = interp(x= y,y= x,z = green)
  xs = fuzWHICH(v1,surface$y)
  ys = fuzWHICH(v2,surface$x)
  green = surface$z[cbind(ys,xs)]
  surface = interp(x= y,y= x,z = blue)
  xs = fuzWHICH(v1,surface$y)
  ys = fuzWHICH(v2,surface$x)
  blue = surface$z[cbind(ys,xs)]
  rgb(red,green,blue,trans)
  #rgb(predict(lm(red   ~ x + y),data.frame(x = v1,y = v2)),
  #    predict(lm(green ~ x + y),data.frame(x = v1,y = v2)),
  #    predict(lm(blue   ~ x + y),data.frame(x = v1,y = v2)))
  
}

nColor <- function(number = 2,base.col = '#AF2F03',trans = 1,distinct = FALSE){
  b.col = rgb2hsv(col2rgb(base.col))
  spread = 1/number
  cols = matrix(NA,3,number)
  cols[,1] = b.col
  if(number > 1){
    for(i in 2:number){
      cols[1,i] = b.col[1]  + (i-1)/number
      cols[2:3,i] = b.col[2:3]
    }
  }
  cols[1,cols[1,] > 1] = cols[1,cols[1,] > 1] - 1
  order = 1:number
  
  if(number > 2 & distinct == TRUE){
    if(number %% 2 == 0){
      order = c()
      one = 1:(number/2)
      two = (number/2+1):number
      for(i in 1:length(one)){
        order = c(order, one[i] , two[i])
      }
    }else{
      order = 1
      one = 2:(floor(number/2)+1)
      two = (floor(number/2)+2):number  
      for(i in 1:length(one)){
        order = c(order, two[i] , one[i])
      }
    }
  }
  
  
  hsv(cols[1,order],cols[2,order],cols[3,order],trans)
  
}

plot.data <- function(x,y,data,...){
  opar <- par()
  par(mar=c(3, 4, 3, 2), oma=c(0,0,0,0), bg="white", 
      xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), 
      las=1, col.axis="#434343", col.main="#343434", 
      tck=0, lend=1,pch=20)
  xx = data[,which(colnames(data) == x)]
  yy = data[,which(colnames(data) == y)]
  plot(xx,yy,bty="n", las=1,xlab=bquote(bold(.(x))), 
       ylab=bquote(bold(.(y))), family="Helvetica", 
       cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", 
       yaxt="n",...)
  grid(NULL, NULL, col="#DEDEDE", lty=3, lwd=0.9)
  axis(1, tick=FALSE, cex.axis=0.9)
  axis(2, tick=FALSE, cex.axis=0.9)
  par(opar)

}


five38.plot <- function(multiLine,minX,maxX,minY,maxY,Title ="",xTit ="",yTit ="",colUsed = "",embolden = 25){ # FiveThirtyEight style many line plot
  opar <- par()
  if(colUsed == "")colUsed <- minX:maxX
  par(mar=c(3, 4, 3, 2), oma=c(0,0,0,0), bg="#F0F0F0", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
  grid(NULL, NULL, col="#DEDEDE", lty="solid", lwd=0.9)
  axis(1, tick=FALSE, cex.axis=0.9)
  axis(2, tick=FALSE, cex.axis=0.9)
  ManyLinePlots(series = multiLine,NameColUsed=colUsed,boldSeries = embolden, col="dark grey", lwd=1, hlwd=3, hcol="#008ED4")
  par(opar)
}

dMvn <- function(X,mu,Sigma,log=F) { #multivariate normal density
  #X values, mu = means, sigma = cov
  k <- ncol(X)
  rooti <- backsolve(chol(Sigma),diag(k))
  quads <- colSums((crossprod(rooti,t(X-mu)))^2)
  if(log==T){
    return(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads)
  }else{
    return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads)) 
  }
}

cart2pol <- function(x, y)
{
  r <- sqrt(x^2 + y^2)
  t <- atan(y/x)
  
  c(r,t)
}

polarTOcart <- function(angle,radius,inner.c=0){#angle in degrees adds an inner radius
  radius <- radius + inner.c
  if(angle <= 90){  
    x <- radius*sin(angle*2*pi/360)
    y <- radius*cos(angle*2*pi/360)}
  if(90 < angle & angle <= 180){  
    x <- radius*cos((angle-90)*2*pi/360)
    y <- -radius*sin((angle-90)*2*pi/360)}
  if(180 < angle & angle <= 270){  
    x <- -radius*sin((angle-180)*2*pi/360)
    y <- -radius*cos((angle-180)*2*pi/360)}
  if(270 < angle & angle <= 360){  
    x <- -radius*cos((angle-270)*2*pi/360)
    y <- radius*sin((angle-270)*2*pi/360)}
  list(x=x,y=y)
}


####radar Plot
radarPlot <- function(values,val.l,val.h,scale.within = F,legend = F,color='black'){
  if(!is.matrix(values)){ 
    values = t(as.matrix(values))
    val.l = t(as.matrix(val.l))
    val.h = t(as.matrix(val.h))
  }
  if(is.null(colnames(values))) colnames(values) = 1:ncol(values)
  
  opar <- par()
  on.exit(par(opar))
  #scale.within = F
  
  #values <- allQ[,-1] 
  #val.l <- allQ.lo[,-1]
  #val.h <- allQ.hi[,-1]
  
  #need to alter colors
  
  max.v <- max(val.h,na.rm=T)
  min.v <- min(val.l,na.rm=T)
  
  dims = c(ceiling(sqrt(nrow(values))),max(1,nrow(values) -ceiling(sqrt(nrow(values)))))
  if(color == 'black'){col1 = 'black';col2 = 'white'}else{col1 = 'white';col2 = 'black'}
  
  par(mfrow=c(dims[1],dims[2]),pty='s',mar=c(1,1,2,1),bg=col1)
  for(i in 1:dim(values)[1]){  #if scaleing within
    
    if(scale.within==T){
      max.v <- max(val.h[i,],na.rm=T)
      min.v <- min(val.l[i,],na.rm=T)
    }
    
    inner.c <- (max.v + abs(min.v))/10
    outer.r <- max.v + abs(min.v) + inner.c + inner.c
    values[is.na(values)] <- 0
    val.l[is.na(val.l)] <- 0
    val.h[is.na(val.h)] <- 0
    npara <- dim(values)[2]
    angles <- 360/npara
    
    polarTOcart <- function(angle,radius,inner.c=0){#angle in degrees adds an inner radius
      radius <- radius + inner.c
      if(angle <= 90){  
        x <- radius*sin(angle*2*pi/360)
        y <- radius*cos(angle*2*pi/360)}
      if(90 < angle & angle <= 180){  
        x <- radius*cos((angle-90)*2*pi/360)
        y <- -radius*sin((angle-90)*2*pi/360)}
      if(180 < angle & angle <= 270){  
        x <- -radius*sin((angle-180)*2*pi/360)
        y <- -radius*cos((angle-180)*2*pi/360)}
      if(270 < angle & angle <= 360){  
        x <- -radius*cos((angle-270)*2*pi/360)
        y <- radius*sin((angle-270)*2*pi/360)}
      list(x=x,y=y)
    }
    
    
    #split.screen(c(4,3))
    #par(mfrow=c(1,1),pty='s',mar=c(2,1,2,1))
    #for(i in 1:dim(values)[1]){  #if scaling amoung
    plot( 1,1,xlim=c(-outer.r,outer.r),ylim=c(-outer.r,outer.r),
          type='n',xlab='',ylab='',main=rownames(values)[i],
          bty='n',xaxt='n',yaxt='n',col.main=col2)
    
    point.mat <- val.mat <- out.mat <- inn.mat <- zer.mat <- matrix(NA,dim(values)[2],2) #hold vertices for later use
    for(j in 1:dim(values)[2]){ #plots values
      if(j == dim(values)[2]){
        one <- polarTOcart(angles*(j-1),values[i,j],inner.c+abs(min.v))
        two <- polarTOcart(angles*0,values[i,1],inner.c+abs(min.v))
        # lines(c(one$x,two$x),c(one$y,two$y),lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=1,lwd=1)
      }else{
        one <- polarTOcart(angles*(j-1),values[i,j],inner.c+abs(min.v))
        two <- polarTOcart(angles*j,values[i,j+1],inner.c+abs(min.v))
        #lines(c(one$x,two$x),c(one$y,two$y),lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=1,lwd=1)
      }
      val.mat[j,] <- c(one$x,one$y)
      
      one <- polarTOcart(angles*(j-1),val.l[i,j],inner.c+abs(min.v)) #plots credible intervals
      two <- polarTOcart(angles*(j-1),val.h[i,j],inner.c+abs(min.v))
      lines(c(one$x,two$x),c(one$y,two$y),lty=1,lwd=4,col=4)
      
      if(j == dim(values)[2]){ #plot outer CI
        one <- polarTOcart(angles*(j-1),val.h[i,j],inner.c+abs(min.v))
        two <- polarTOcart(angles*0,val.h[i,j],inner.c+abs(min.v))
        #lines(c(one$x,two$x),c(one$y,two$y),lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=1,lwd=1)
      }else{
        one <- polarTOcart(angles*(j-1),val.h[i,j],inner.c+abs(min.v))
        two <- polarTOcart(angles*j,val.h[i,j],inner.c+abs(min.v))
        #lines(c(one$x,two$x),c(one$y,two$y),lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=1,lwd=1)
      }  
      out.mat[j,] <- c(one$x,one$y)
      
      if(j == dim(values)[2]){ #plot inner CI
        one <- polarTOcart(angles*(j-1),val.l[i,j],inner.c+abs(min.v))
        two <- polarTOcart(angles*0,val.l[i,j],inner.c+abs(min.v))
        #lines(c(one$x,two$x),c(one$y,two$y),lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=1,lwd=1)
      }else{
        one <- polarTOcart(angles*(j-1),val.l[i,j],inner.c+abs(min.v))
        two <- polarTOcart(angles*j,val.l[i,j],inner.c+abs(min.v))
        #lines(c(one$x,two$x),c(one$y,two$y),lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=1,lwd=1)
      }  
      inn.mat[j,] <- c(one$x,one$y)
      
      if(j == dim(values)[2]){ #plots inner circle
        one <- polarTOcart(angles*(j-1),inner.c)
        two <- polarTOcart(angles*0,inner.c)
        lines(c(one$x,two$x),c(one$y,two$y)) 
      }else{
        one <- polarTOcart(angles*(j-1),inner.c)
        two <- polarTOcart(angles*j,inner.c)
        lines(c(one$x,two$x),c(one$y,two$y))       
      }
      point.mat[j,] <- c(one$x,one$y)
      
      if(j == dim(values)[2]){ #plots zerolines
        one <- polarTOcart(angles*(j-1),inner.c+abs(min.v))
        two <- polarTOcart(angles*0,inner.c+abs(min.v))
        lines(c(one$x,two$x),c(one$y,two$y),lty=2,col='red') 
      }else{
        one <- polarTOcart(angles*(j-1),inner.c+abs(min.v))
        two <- polarTOcart(angles*j,inner.c+abs(min.v))
        lines(c(one$x,two$x),c(one$y,two$y),lty=2,col='red')       
      }
      zer.mat[j,] <- c(one$x,one$y)
      
      if(j == dim(values)[2]){ #plots outer lines
        one <- polarTOcart(angles*(j-1),outer.r)
        two <- polarTOcart(angles*0,outer.r)
        lines(c(one$x,two$x),c(one$y,two$y),lty=1,col=col2,lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=2) 
      }else{
        one <- polarTOcart(angles*(j-1),outer.r)
        two <- polarTOcart(angles*j,outer.r)
        lines(c(one$x,two$x),c(one$y,two$y),lty=1,col=col2,lwd=2) 
        #lines(c(one$x,0),c(one$y,0),lty=2)
      }
      
      if(values[i,j] == 0)next
      par(xpd=T)
      if(colnames(values)[j]=='SL'){
        text(one$x*(outer.r*1.05)/outer.r,one$y*(outer.r*1.1)/outer.r,GDsGS[i],srt=0 - angles*(j-1),col=col2)  
      }else{
        text(one$x*(outer.r*1.05)/outer.r,one$y*(outer.r*1.1)/outer.r,colnames(values)[j],srt=0 - angles*(j-1),col=col2)
      }
      par(xpd=F)
    }
    polygon(x=as.numeric(val.mat[,1]),y=as.numeric(val.mat[,2]),col=rgb(.2,.2,.5,.5),border=NA)
    #polygon(x=as.numeric(out.mat[,1]),y=as.numeric(out.mat[,2]),col=rgb(0,0,1,.25),border=NA)
    #polygon(x=as.numeric(out.mat[,1]),y=as.numeric(inn.mat[,2]),col='white',border=NA)
    #polygon(x=as.numeric(zer.mat[,1]),y=as.numeric(out.mat[,2]),col=rgb(1,0,0,.25),border=NA)
    polygon(x=as.numeric(point.mat[,1]),y=as.numeric(point.mat[,2]),col=col1,border=NA) 
  }
  
  ###legend
  if(legend==T){
    plot( 1,1,xlim=c(-outer.r,outer.r),ylim=c(-outer.r,outer.r),
          type='n',xlab='',ylab='',
          bty='n',xaxt='n',yaxt='n')
    point.mat <- zero.mat <- matrix(NA,dim(values)[2],2) #hold vertices for later use
    for(j in 1:dim(values)[2]){ 
      if(j == dim(values)[2]){ #plots inner circle
        one <- polarTOcart(angles*(j-1),inner.c)
        two <- polarTOcart(angles*0,inner.c)
        lines(c(one$x,two$x),c(one$y,two$y))
      }else{
        one <- polarTOcart(angles*(j-1),inner.c)
        two <- polarTOcart(angles*j,inner.c)
        lines(c(one$x,two$x),c(one$y,two$y))
      }
      point.mat[j,] <- c(one$x,one$y)
      
      if(j == dim(values)[2]){ #plots zerolines
        one <- polarTOcart(angles*(j-1),inner.c+abs(min.v))
        two <- polarTOcart(angles*0,inner.c+abs(min.v))
        lines(c(one$x,two$x),c(one$y,two$y),lty=2,col='red') 
      }else{
        one <- polarTOcart(angles*(j-1),inner.c+abs(min.v))
        two <- polarTOcart(angles*j,inner.c+abs(min.v))
        lines(c(one$x,two$x),c(one$y,two$y),lty=2,col='red')       
      }
      text(one$x,one$y,colnames(values)[j],srt=90-angles*(j-1))
      zero.mat[j,] <- c(one$x,one$y)
      
      if(j == dim(values)[2]){ #plots outer lines
        one <- polarTOcart(angles*(j-1),outer.r)
        two <- polarTOcart(angles*0,outer.r)
        lines(c(one$x,two$x),c(one$y,two$y),lty=2) 
        lines(c(one$x,0),c(one$y,0),lty=2) 
      }else{
        one <- polarTOcart(angles*(j-1),outer.r)
        two <- polarTOcart(angles*j,outer.r)
        lines(c(one$x,two$x),c(one$y,two$y),lty=2) 
        lines(c(one$x,0),c(one$y,0),lty=2)
      }
      text(one$x*(outer.r*1.05)/outer.r,one$y*(outer.r*1.1)/outer.r,colnames(values)[j],srt=0 - angles*(j-1))
    }
    polygon(x=as.numeric(point.mat[,1]),y=as.numeric(point.mat[,2]),col=1)
  }
}



#pictogram
pictogram <- function(data,size.gap = 3, dtree, ctree = NULL,aspect = .7,Main.Title = '',xTitle = '', yTitle='',ypos=NA,Unit.Name = ''){ #plot a pictogram
  #data: vector of counts with names
  #size.gap: fraction of column height is gap fraction = icon.height/size.gap, 
  #dtree: primary icon.png 
  #ctree: secondary icon.pnd
  #aspect:width/height of icon
  #ypos: positions the ylabel vertically [0,1] default is the highest column
  #make sure icon is croped to the edges
  #dtree<-readPNG("/home/mckwit/Downloads/decid45-64.png") #0.7 X 1.0
  #ctree<-readPNG("/home/mckwit/Downloads/conif45-64.png")
  #ctree<-readPNG("/home/mckwit/Downloads/sun.png")#1X1 
  #ctree<-readPNG("/home/mckwit/Downloads/pollution.png")#1X1 
  #ctree<-readPNG("/home/mckwit/Downloads/factory.png")#1X1 
  #ctree<-readPNG("/home/mckwit/Downloads/man28X64.png")# 28/64
  #ctree<-readPNG("/home/mckwit/Downloads/car95.png")#
  require(png)  
  ndata = length(data)
  height = 1/ndata 
  gap = height / size.gap
  height = height - gap
  width = aspect*height
  
  xs = c(0,width)
  ys = c(0,height)
  xshift = width
  yaxe = seq(height/2,1,by = height+gap)[1:ndata]
  yshift = height + gap
  
  num = floor(1/xshift)
  size = ceiling(max(data)/num)
  ns = data/size
  Title = Main.Title
  xTit = yTitle
  yTit = xTitle
  
  par(mar=c(3, 4, 3, 2), oma=c(0,0,0,0), bg="#F0F0F0", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
  plot(0,0,xlim=c(0,1),ylim=c(0,1),type='n',xaxt='n',yaxt='n',bty='n',las=1,main=Title, xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8)
  
  if(is.null(ctree))ctree = dtree
  for(j in 1:length(data)){
    for(i in 1:ceiling(ns[j])){
      if(i %% 2 != 0){icon = dtree}else{icon = ctree}
      rasterImage(image=icon,xleft=xs[1]+xshift*(i-1),xright=xs[2]+xshift*(i-1),ybottom=ys[1]+yshift*(j-1),ytop=ys[2]+yshift*(j-1)) # single pack 
    }
    rect(col = "#F0F0F0",border="#F0F0F0",
         xleft=xs[1]+xshift*(i-1)+(ns[j]%%1)*diff(xs),xright=xs[2]+xshift*(i-1),ybottom=ys[1]+yshift*(j-1),ytop=ys[2]+yshift*(j-1))
  }
  grid(NULL, NULL, col="#DEDEDE", lty="solid", lwd=0.9)
  axis(1,at = c(0,(max(ns)/2)*xshift,max(ns)*xshift),labels=c(0,floor(max(data)/2),max(data)), cex.axis=0.9)
  axis(2,at =yaxe,labels=names(data),tick=F, cex.axis=0.9,line=-.4)
  
  if(is.na(ypos))ypos = ys[1]+yshift*(j-1)
  icon = dtree
  zoom = 0.9
  a = strwidth('Each ', cex = zoom)
  b = strwidth(paste(' equals',size,Unit.Name), cex = zoom)
  xpos = 1-a-b - width
  rasterImage(image=icon,xleft=xpos,xright=xpos+width,ybottom=ypos,ytop=ypos+height)
  text(xpos,ypos+height/2,labels = 'Each ',adj=c(1,0.5))
  text(xpos+width,ypos+height/2,labels = paste(' equals',size,Unit.Name),adj=c(0,.5))
  
}

squarePie <- function(data,size.gap = 10,Main.Title = 'Main Title',Sub.Title=NA,Unit.Name = 'Units',nh = 7,nw = 34){ #plot a pictogram
  
  #data: vector of counts with names
  #size.gap: fraction of column height is gap fraction = icon.height/size.gap, 

  data = rev(sort(data))
  ndata = length(data)
  colours = nColor(number=ndata+1,distinct=T)
  
  colors = c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
  hsvC=rgb2hsv(col2rgb(colors))
  colours = c("#00A0B0","#6A4A3C","#CC333F","#EDC951","#37C040")
  colours = colorRange(percent = seq(0,1,length=ndata+1),colours = colours,trans = 1)
  width = .7
  height = width/4
  n = nh*nw
  unit = round(sum(data,na.rm=T) / n)
  allNum = round(data/unit)
  cumNum = cumsum(allNum)
  
  iSide=min(c(height/nh,width/nw))
  gap = iSide / size.gap
  
  iSide = iSide - gap
  shift = iSide + gap  
  
  xs = c(shift,shift+iSide)
  ys = c(.45,.45+iSide)
  dotted = c(1,ceiling(cumNum/nh)+1)

  par(mar=c(0, 0, 0, 0), oma=c(0,0,0,0), pty='s',bg="white", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
  plot(0,0,xlim=c(0,1),ylim=c(0,1),type='n',xaxt='n',yaxt='n',bty='n',las=1,xlab="", ylab="", family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8)
  rect(xleft = 0, xright = shift*(nw + .4*nw) + xs[1], 
       ybottom = ys[1] - shift*3.5 - 2*strheight(names(data)[1],cex=.75), ytop =  shift*(nh) + ys[2] + strheight(Main.Title) + shift,
       col = "#E8F4F8" ,border="#E8F4F8")
  k = 0;l=1
  for(j in 1:(nw+3)){
    if(j==dotted[l]){
      if(l %% 3 == 0){mul = 2.6}else{if(l %% 3 == 2){mul = 1.3}else{mul = 0}}
      lines(x = rep(xs[1]+shift*(j-1+.5),2),y = c(ys[1]-shift*.1,ys[1] - shift*0.8 - mul*strheight(names(data)[1],cex=.75)),lty=3,lwd=1,col="#000026")
      text(xs[1]+shift*(j-1),ys[1] - shift - mul*strheight(names(data)[1],cex=.75),names(data)[l],cex=.75,adj=c(0,.9),col="#000026")
      l=l+1
    }
    for(i in 1:nh){
      k = k + 1
      colNum = rep(1:ndata,allNum)
      rect(xleft=xs[1]+shift*(j-1),xright=xs[2]+shift*(j-1),ybottom=ys[1]+shift*(i-1),ytop=ys[2]+shift*(i-1),
           col = colours[colNum[k]] ,border="#E8F4F8") # single pack 
      if(k == sum(allNum))break
    }
    if(k == sum(allNum))break  
  }
  text(shift, 
       shift*(nh) + ys[2] ,
       labels = Main.Title,
       pos=4,offset=0,col="#000026")
  
  sub.width = (shift*(nw + .4*nw) + xs[1]) - (xs[2]+shift*(j-1)) - shift
  if(is.na(Sub.Title)){
    Sub.Title = str_wrap(paste(prettyNum(sum(data,na.rm=T),big.mark=",",scientific=F),'total', Unit.Name), width = floor(2*sub.width/strwidth('me', cex = 1)), indent = 0, exdent = 0)
  }else{
    Sub.Title = str_wrap(Sub.Title, width = floor(2*sub.width/strwidth('me', cex = 1)), indent = 0, exdent = 0)
  }
    text(shift*(nw + .4*nw) + xs[1] - shift, 
       shift*(nh) + ys[2] + strheight(Main.Title),
       labels = Sub.Title,
       adj=c(1,1),col="#000026")
  
  rect(xleft =  (shift*(nw + .4*nw) + xs[1]) - (2*shift + strwidth(paste(unit,Unit.Name))/2) ,
       xright = (shift*(nw + .4*nw) + xs[1]) - (shift + strwidth(paste(unit,Unit.Name))/2),
       ybottom=ys[1],ytop=ys[2],
       col = colours[ndata+1] ,border="#E8F4F8") 
  lines(x = rep((shift*(nw + .4*nw) + xs[1]) - (1.5*shift + strwidth(paste(unit,Unit.Name))/2),2),
        y = c(ys[1]-shift*.1,ys[1] - shift*.8),lty=3,lwd=1,col="#000026")
  text(x = (shift*(nw + .4*nw) + xs[1]) - (1.5*shift + strwidth(paste(unit,Unit.Name))/2),
       y = ys[1] - shift,
       paste(prettyNum(unit,big.mark=",",scientific=F),Unit.Name),cex=.75,col="#000026")
  
  
}

#short as .numeric
a.n <- function(x){as.numeric(x)}#short as.numeric

#short as.character
a.c <- function(x){as.character(x)}#short as.character

##shorten Length
len <- function(x){length(x)}

#return last value
last <- function(vec){rev(vec)[1]}

#turns ones to zeros zeros to ones
invert <- function(x){
  whr = which(x == 1)
  x[x == 0] = 1
  x[whr] = 0
  x
}

#split vector string and return matrix
unlistSplit = function(vec,sep="-"){
  matrix(unlist(strsplit(vec,split=sep)),nrow = length(vec),byrow=T)
}

getScoreNorm <- function(x,mu,xvar){ #proper predictive score
  # Gneiting and Raferty's proper scoring rule
  # outcome x, prediction mean and variance (mu, xvar)
  #higher scores are better
  - ( (x - mu)^2)/xvar - log(xvar)
}

#
# Splining a polygon.
#
#   The rows of 'xy' give coordinates of the boundary vertices, in order.
#   'vertices' is the number of spline vertices to create.
#              (Not all are used: some are clipped from the ends.)
#   'k' is the number of points to wrap around the ends to obtain
#       a smooth periodic spline.
#
#   Returns an array of points. 
# 
spline.poly <- function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}

#distance between a point and a matrix
distance = function(valM,M){
  dims = length(valM)
  resp = rep(0,dim(M)[1])
  for(i in 1:dims){
    resp = resp + (M[,i]-valM[i])^2
  }
  sqrt(resp) 
}


angleTo <- function(valM,M){ #output in degrees
  
  tx = M[,1]-valM[1]
  ty = M[,2]-valM[2]
  ang = atan(ty/tx) * 360/(2*pi)
  tx = (tx > 0) - (tx < 0) 
  ty = (ty > 0) - (ty < 0) 
  list(ang = ang,tx=tx,ty=ty)
}


#
spatialQuant = function(xyMat,nbrks = 10,quants = c(.025,.975)){
  stdX = sd(xyMat[,1],na.rm=T)
  stdY = sd(xyMat[,2],na.rm=T)
  
  
  xyMat[,1] = (xyMat[,1]) /stdX
  xyMat[,2] = (xyMat[,2]) /stdY
  
  mx = median(xyMat[,1])
  my = median(xyMat[,2])
  
  tmp = angleTo(c(mx,my),xyMat)
  angs = tmp$ang
  tx = tmp$tx
  ty = tmp$ty
  quad = rep(0,len(tx))
  quad[tx ==  1 & ty ==  1] = 1
  quad[tx ==  1 & ty == -1] = 1
  quad[tx == -1 & ty == -1] = -1
  quad[tx == -1 & ty ==  1] = -1
  
  dis = distance(c(mx,my),xyMat)
  
  
  breaks <- seq(-90,90,length = nbrks)
  bins <- cut(angs,breaks = breaks)
  
  iBin = sort(unique(bins))
  tmp = matrix(NA,len(iBin),3)
  tmp[,1] = breaks[-1] - diff(breaks)[1]/2
  for(i in 1:len(iBin)){
    whr = which(bins == iBin[i])
    tmp[i,c(2,3)] = quantile(dis[whr] * quad[whr],quants)
  }
  xs = (cos(c(tmp[,1],tmp[,1])*2*pi/360)*c(tmp[,2],tmp[,3]) + mx)*stdX 
  ys = (sin(c(tmp[,1],tmp[,1])*2*pi/360)*c(tmp[,2],tmp[,3]) + my)*stdY 
  cbind(xs,ys)
}


#stretch to fit new range
stretch <- function(dat,finMax,finMin,datMax,datMin){#stretch to fit new range
  
  res = finMax - (finMax - finMin) * (datMax- dat)/(datMax-datMin)
  res
}

designCombn <- function(mat){
  comb = numeric()
  for(i in 1:ncol(mat)){
    
    whr = t(combn(ncol(mat),ncol(mat)-i+1))
    for(j in 1:nrow(whr)){
      tmp = rep(0,ncol(mat))  
      tmp[whr[j,]] = 1
      comb = rbind(comb,tmp)
      
    }
  }
  comb
}

allLm <- function(resp,mat,Int = T){
  combs = designCombn(mat)
  if(Int == T){
    combs = cbind(rep(1,nrow(combs)),combs)
    mat = cbind(rep(1,nrow(mat)),mat)
  }
  ICs = numeric()
  for(i in 1:nrow(combs)){
    whr = which(combs[i,] == 1)
    tmp = as.matrix(mat[,whr])
    fit = lm(resp ~ tmp - 1 , y=T)
    combs[i,whr] = coefficients(fit)
    fitted = predict(fit, se.fit = TRUE)
    ps = mean(getScoreNorm(fit$y,fitted.values(fit),fitted$se.fit^2),na.rm=T)
    ICs = rbind(ICs, c(AIC(fit),BIC(fit),ps)) 
    
  }
  colnames(ICs) = c('AIC','BIC','PS')
  w = (1/ICs[,'PS'])/ sum(1/ICs[,'PS'])
  avgB = apply(combs,2,function(x) weighted.mean(x,w))
  list(coefs = combs,ICs = ICs, avgB = avgB)
}

fillNAs <- function(mat,resp = mat[,1]){
  
  for(i in 1:ncols(design)){
    
    resp = design[,i]
    allLm(resp,design)
    
  }
  
}

#progress bar
progress = function(count, total = null, interval = 1,startTime = NA){ #progress bar
  #progress bar
  #count: counter of loop
  #total: number when loop completes
  #interval: interval to update counter
  width = getOption("width")/2 - 4
  if(count%%interval == 0){
    Ltotal = nchar(total)
    Lcount = nchar(count)
    if(is.na(startTime)){
      cat("\r",paste(count, rep(" ",Ltotal-Lcount),
                   " [",
                   paste(rep("~",floor(width*count/total)),collapse = ""),
                   paste(rep(" ",ceiling(width*(1-(count/total)))),collapse = ""),
                   "] ",round(100*count/total,2),"%",sep=""))
    }else{
      print((total-count)*(Sys.time() - startTime)/count)
    }
    flush.console()    
  }
  if(count == total){cat("\n Complete \n")}
}

#print(paste("It has been", as.numeric(as.Date(Sys.time()) - as.Date('1980-11-09 04:06:00') ),"days."))

matchAllCol = function(mat,vect,TorF = FALSE){

  if(ncol(mat) != length(vect))return("dimensions don't match")
  apply(mat,1,function(x) paste(x,collapse=""))
  if(TorF){
    apply(mat,1,function(x) paste(x,collapse="")) == paste(vect,collapse="") 
  }else{
    which(apply(mat,1,function(x) paste(x,collapse="")) == paste(vect,collapse=""))
  }  
}

#percent overlap of pop a and b
overlap = function(a,b){
  n = max(length(a),length(b),10000)
  a = sample(a,n,replace = T) #make samples equal size
  b = sample(b,n,replace = T)
  a95 = quantile(a,.95)
  b95 = quantile(b,.05)
  c = c(a,b)
  overlap = (length(which(b95 < a ))+length(which(b < a95 ))) / length(c) #how much overlap between 95% and other distribution twice
  aOverB = length(which(b95 < a )) / length(a) #how much of a is greater than 5%b
  bOverA = length(which(b < a95 )) / length(b) #how much of b is less than 95%a
  overlap95 = (length(which(b95 < c & c < a95 ))) / length(c) #how much overlap between 95% and 95% of other distribution
  
  list(overlap = overlap,aOverB = aOverB,bOverA = bOverA,
        overlap95 = overlap95, probOverlap = aOverB*bOverA, sampOver = 1- (length(which(a<b))/length(a)))
}

#square root of matrix a
sqrtMat <- function(a){  #square root of matrix a
  a.eig <- eigen(a)
  a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
}

#sqaure root of matrix
denman.beavers <- function(mat,maxit=50) {
  stopifnot(nrow(mat) == ncol(mat))
  niter <- 0
  y <- mat
  z <- diag(rep(1,nrow(mat)))
  for (niter in 1:maxit) {
    y.temp <- 0.5*(y+solve(z))
    z <- 0.5*(z+solve(y))
    y <- y.temp
  }
  return(list(sqrt=y,sqrt.inv=z))
}


#random wishart from clark
rwish <- function(df,S){
  
  z  <- matrix(rnorm(df*nrow(S)),df,nrow(S))%*%chol(S)
  crossprod(z)
}

#random inverse wishart from clark
riwish <- function(v,S){
  
  solve(rwish(v,solve(S)))
}

#pastes full name of file address using paste(path, "add",sep="") 
paths = function(path = path, address){ #
  
  paste(path,address,sep="")
  
}

## Size combines dim and length
sz <- function(x){
  if(is.null(dim(x))){ 
    len(x)
  }else{
    dim(x)
  }
}

#bind rows don't worry about matching names
rbind.all <- function(x1,x2,...){
  t(cbind(x1,x2,...))
} 

#Function reads zip into temp folder returns paths 
pathZip <- function(zipName,fileName = character(),readFiles=T){#Function reads zip into temp folder returns paths 
  #Function reads zip into temp folder 
  #Returns path of files
  #zipName path of zip
  #fileName: name of files in zip to return path if blank reuturns all paths
  zipdir = tempdir()
  names = unzip(zipName, exdir=zipdir, list = T)
  print(names)
  names = names[,1]
  unzip(zipName, exdir=zipdir, list = F)
  if(len(fileName) == 0 ) fileName = names
  list(tmpPath = zipdir, fileNames = names, filePath = file.path(zipdir,fileName))
}

#Converts rgb to hex
rgb2hex <- function(x){ #converts rgb to hex alpha is [0:256]
  require(stringr)
  #x = as.matrix(x)
  hex <- as.hexmode(x)
  hex <- as.character(hex)
  whr = which(str_length(hex) == 1)
  hex[whr] = paste0("0",hex[whr] )
  if(length(dim(x)) == 2){
    hex <- matrix(hex, ncol = ncol(x),byrow=T)
  }else{
    hex <- matrix(hex, ncol = length(x),byrow=T)
  }
  hex <- apply(hex,1, function(x){paste0(x,collapse = '')})
  hex <- paste0('#',hex)
  hex
}

#list all functions in package
lsp <-function(package, all.names = FALSE, pattern) {
  package <- deparse(substitute(package))
  ls(
    pos = paste("package", package, sep = ":"),
    all.names = all.names,
    pattern = pattern
  )
}


##Custom legend for continuous data meant for a map
myLegend <- function(x,y,h,w,colours = c('#ffffb2','#bd0026'),padWhite = T,highVal = "",lowVal =""){
  #x: left x all in plot scale
  #y: bottom y
  #h: height
  #w: width 
  #colours: colors in used in colorRange() 
  #padWhite: Should scale be padded with white at bottom
  #highVal lowVal: vlaues to print on legend
  
  hs = seq(y,y+h,length=100)
  cols = colorRange((hs-y)/max((hs-y)),colours = colours,trans = 1)
  if(padWhite){cols = c('white','white','white','white','white',cols)}
  for(i in 1:len(cols)) lines(c(x,x+w),c(hs[i],hs[i]),col=cols[i],lwd=1) 
  rect(x,y,x+w,y+h,border=rgb(.8,.8,.8,.5))
  text(x+w,max(hs),paste(highVal),adj=c(-0.1,1))
  text(x+w,min(hs),paste(lowVal),adj=c(-0.1,0))
}
