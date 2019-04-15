#############################
# 1/19/2015
# R script: New Boston Air Force Station
#    KTS: 1/19/2015:  start working on script for testing with preliminary data 
#    KTS: 4/12/2015:  start analysis with full data set

#############################
# PREPARE THE WORKSPACE
#############################

rm(list=ls())                  # clear the workspace
options(scipen=6)              # for viewing clarity, set scientific notation to occur only when more than six decimal places

#############################
#  FUNCTIONS
#############################

FtoC <- function(F){
  (F - 32) * 5/9
}

IntoCm <- function(In){
  if(In!="N/A"){
    as.numeric(In)/0.39370
  } else{
    0
  }
}

MPHtoMPS <- function(MPH){
  if(MPH!="calm"){
    as.numeric(MPH) * 0.44704
  } else{
    0.01
  }
} 

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

CondToCond <- function(Cond){
  ndx <- which(Cond==CondReclass$Original)
  CondReclass$Reclass[ndx]
}

##### function for summarizing total effective captures from time difference data
computecaps <- function(timedif){
  if(length(timedif)<=1){
    count <- length(timedif)
  }else{
    count <- 0
    runsum <- 16
    for(i in 1:length(timedif)){
      runsum <- runsum + timedif[i]
      if(runsum>15){       # if time difference is greater than 15 minutes!
        count <- count + 1
        runsum <- 0
      }
    }  
  }
  return(count)
}


#############################
#   LOAD PACKAGES
#############################

########################################
## GENERIC FUNCTIONs FOR INSTALLING/LOADING PACKAGES FROM CRAN AND GITHUB
########################################

loadPackage <- function(pkg){
  if(pkg %in% rownames(installed.packages()) == FALSE) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)
}

source_github <- function(baseurl,scriptname) {
  # load package
  #suppressMessages(suppressWarnings(require(RCurl)))
  loadPackage("RCurl")
  
  # read script lines from website
  url <- sprintf("%s%s",baseurl,scriptname)
  script <- getURL(url, ssl.verifypeer = FALSE)
  
  script <- gsub("\r\n", "\n", script)     # get rid of carriage returns (not sure why this is necessary...)
  
  # parse lines and evaluate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}



#############################
#   LOAD PACKAGES
#############################

loadPackage("RJSONIO")
loadPackage("R2WinBUGS")      # interface with WinBUGS
loadPackage("emdbook")         # "emdbook" package: has a useful function for converting BUGS output to CODA output (eg. for convergence diagnostics)
loadPackage("plyr")
loadPackage("devtools")
loadPackage("maptools")
loadPackage("weatherData")
loadPackage("lubridate")       # "lubridate" package: facilitates working with dates
loadPackage("abind")           # for combining multidimensional arrays
loadPackage("Hmisc")           # for displaying error bars among other things
loadPackage("oce")
loadPackage("ROCR")
loadPackage("RColorBrewer")
#loadPackage("partykit")      # try hothorn's new package



#############################
# SET USER AND DIRECTORIES
#############################

KEVIN = FALSE
OFFICE = TRUE

         # set base directory for storing results, images, etc.
if (KEVIN) baseDir <- "C:\\Users\\Kevin\\Dropbox\\Scott Reynolds\\New Boston\\"
if (OFFICE) baseDir <- "E:\\Dropbox\\Scott Reynolds\\New Boston\\"

imageDir <- paste(baseDir,"rawSVG",sep="")         # directory for storing raw image files (in SVG format)
resultsDir <- paste(baseDir,"Results",sep="")      # '' for results files
dataDir <- paste(baseDir,"Data",sep="")            # '' for raw data
dumpDir <- paste(baseDir,"Dump",sep="")            # '' for debug messages etc
BUGSDir <- paste(baseDir,"BUGS",sep="")            # '' for BUGS code

BUGSloc <- "C:\\WinBUGS\\winbugs14\\WinBUGS14" 

###########################
# SET MCMC/GIBBS SAMPLER SYSTEM
############################

ni <- 20000   # number of draws from the posterior
nt <- 10      # thinning rate
nb <- 10000    # number to discard for burn-in
nc <- 1      # number of chains

n.cores <- 1  


############################
# OTHER GLOBAL VARS
############################

Acoustics_Versions <- c("_OLD","_NEW")
Acoustics_Version <- Acoustics_Versions[2]     # use the new acoustics data

Acoustics_Datefields <- c("Adjusted.Date","Date")
Acoustics_Datefield <- Acoustics_Datefields[2]

Acoustics_Specfields <- c("Prominent.Species","SPECIES.1","SPECIES.2")
Acoustics_Specfield <- Acoustics_Specfields[3]

###########################
## RECLASS WEATHER CONDITIONS 
###########################

setwd(dataDir)
file="condMap.csv"
sink(file)
cat("
    Original,Reclass
    Clear,Clear
    Mostly Cloudy,Cloudy
    Partly Cloudy,Clear
    Overcast,Cloudy
    Scattered Clouds,Clear
    Haze,Clear
    Light Rain,Rain
    Heavy Rain,Rain
    Rain,Rain
    Fog,Cloudy
    Heavy Thunderstorms and Rain,Rain
    Light Thunderstorms and Rain,Rain
    Thunderstorm,Rain
    Light Drizzle,Rain
    Thunderstorms and Rain,Rain
    Drizzle,Rain
    Snow,Snow
    Blowing Snow,Snow
    Light Snow,Snow
    Heavy Snow,Snow
    Light Freezing Rain,Rain
    Unknown,NA
    Freezing Rain,Rain
    Light Freezing Fog,Cloudy
    Light Ice Pellets,Rain
    Ice Pellets,Rain
    Mist,Cloudy
    Light Freezing Drizzle,Rain
    Small Hail,Rain
    ")
sink()

CondReclass <- read.csv(file,header=T)
CondReclass$Reclass <- trim(CondReclass$Reclass)
CondReclass <- CondReclass[-nrow(CondReclass),]

############################
#   MAP ALTERNATIVE SPECIES NAMES
############################

setwd(dataDir)
file="specMap.csv"
sink(file)
cat("
    newname,oldname
    LABO,Lbore
    EPFU,Efusc
    MYLU,Mlucy
    LACI,Lcine
    LANO,Lnoct
    MYLE,Mleib
    MYSE,Msept
    PESU,Psubf
    MYOTIS,Myotis
    ALLMYOTIS,AllMyotis
    ALLTREE,AllTree
    ")
sink()

NameReclass <- read.csv(file,header=T)
NameReclass$oldname <- trim(NameReclass$oldname)
NameReclass$newname <- trim(NameReclass$newname)
NameReclass <- NameReclass[-nrow(NameReclass),]


############################
#   READ IN ALL ACOUSTIC DATA
############################

setwd(dataDir)
acousticEffort.df <- read.csv("AcousticEffort.csv",header=T)
acousticEffort.df$DATE <- as.character(acousticEffort.df$First.Night)
acousticEffort.df$DATE2 <- mdy(acousticEffort.df$DATE,tz = "EST")
acousticEffort.df$MONTH <- month(acousticEffort.df$DATE2)
acousticEffort.df$DATETIME <- acousticEffort.df$DATE2 + hours(0)

acousticdatafile <- sprintf("AcousticData%s.csv",Acoustics_Version)
acousticData.df <- read.csv(acousticdatafile,header=T)
acousticData.df$DATE2 <- ymd(as.character(acousticData.df$Date),tz = "EST")
acousticData.df$YEAR <- year(acousticData.df$DATE2)
acousticData.df$MONTH <- month(acousticData.df$DATE2)
acousticData.df$TIME2 <- hms(as.character(acousticData.df$Time))
acousticData.df$DATE2[acousticData.df$TIME2<hm("10:30")] <- acousticData.df$DATE2[acousticData.df$TIME2<hm("10:30")] + days(1)     # turn into REAL date....
acousticData.df$DATETIME <- acousticData.df$DATE2 + acousticData.df$TIME2
acousticData.df$DATE <- as.character(format(acousticData.df$DATE2, "%Y-%m-%d"))
acousticData.df$Prominent.Species <- gsub(". ","",acousticData.df$SPECIES.2)
acousticData.df$SPECIES <- as.character(acousticData.df$Prominent.Species)
acousticData.df$SPECIES[acousticData.df$Prominent.Species=="Unknown"] <- NA

head(acousticEffort.df)
head(acousticData.df)

table(acousticData.df$DATE)



############################
#   READ IN MISTNET DATA
############################

setwd(dataDir)
netEffort.df <- read.csv("MistnetEffort.csv",header=T)
    #acousticEffort.df <- read.csv("AcousticEffort.csv",header=T)
netData.df <- read.csv("MistnetData.csv",header=T)
    #acousticData.df <- read.csv("AcousticData.csv",header=T)

head(netEffort.df)
     #head(acousticEffort.df)
head(netData.df)
    #head(acousticData.df)

### total recaptures...
sum(netData.df$RECAPT)      # 7 total recaptures!   can't do CMR analysis!

##########
### process and summarize capture data


netData.df$DATE2 <- mdy(netData.df$DATE)
netData.df$YEAR <- year(netData.df$DATE2)
netData.df$DATE <- as.character(netData.df$DATE)                                   
netData.df$NET.CODE <- as.numeric(netData.df$NET.CODE)              
netData.df$SPECIES <- as.character(netData.df$SPECIES)
netData.df$TIME <- as.character(netData.df$TIME)
netData.df$TIME2 <- hm(netData.df$TIME)

# remove spaces from the species names
netData.df$SPECIES <- gsub(". ","",netData.df$SPECIES)

netEffort.df$DATE <- as.character(netEffort.df$First.Night) 
netEffort.df$DATE[which(netEffort.df$DATE=="6/31/2002")] <- "7/1/2002"

netEffort.df$NET.HEIGHT <- as.factor(netEffort.df$NET.HEIGHT) 
netEffort.df <- subset(netEffort.df,select=-c(First.Night,X,X.1)) 

netEffort.df$NET.HEIGHT <- as.character(netEffort.df$NET.HEIGHT)
netEffort.df$NET.HEIGHT[which(netEffort.df$NET.HEIGHT=="C")] <- "4"  
netEffort.df$NET.HEIGHT <- as.numeric(netEffort.df$NET.HEIGHT)    


table(acousticData.df$Prominent.Species,year(ymd(acousticData.df$DATE,tz="EST")))
table(netData.df$SPECIES,year(mdy(netData.df$DATE)))

### summarize 2012- only year with both acoustic and mistnet data

ssa <- subset(acousticData.df,YEAR==2012)
ssn <- subset(netData.df,YEAR==2012)

table(ssa$MONTH)

subset(ssa,MONTH==6)

################
#  daily weather data
################

setwd(dataDir)
filename <- "ManchesterAirportWeatherStationDaily.csv"
weatherData.df <- read.csv(filename,header=T)

head(weatherData.df)

############################
#   SET GLOBAL PARAMS
############################
 
maxAbund <- 500      # only applies to acoustics data

activitySeason <- c(5,6,7,8,9)    #c("May","June","July","August","September")

beginRoosting <- "5/15"
endRoosting <- "8/15"


#####

head(netEffort.df)
head(acousticEffort.df)
head(netData.df)
head(acousticData.df)

######################
#  Process weather data

head(weatherData.df)
weatherData.df$DATE = ymd(weatherData.df$DATE,tz="EST")
weatherData.df$PRCP = ifelse(weatherData.df$PRCP==-9999,NA,weatherData.df$PRCP/10)
weatherData.df$PRCP = ifelse(weatherData.df$PRCP==-9999,NA,weatherData.df$PRCP/10)
weatherData.df$TMAX = ifelse(weatherData.df$TMAX==-9999,NA,weatherData.df$TMAX/10) 
weatherData.df$TMIN = ifelse(weatherData.df$TMIN==-9999,NA,weatherData.df$TMIN/10)
weatherData.df$AWND = ifelse(weatherData.df$AWND==-9999,NA,weatherData.df$AWND/10)
weatherData.df$WSF2 = ifelse(weatherData.df$WSF2==-9999,NA,weatherData.df$WSF2/10)

head(weatherData.df)

#################################
#     PROCESS AND COMBINE DATA FOR ANALYSIS
#################################

table(netData.df$SPECIES)
table(acousticData.df$Prominent.Species)

#########
  ##   find unique species that are in both datasets

 ### unique species (mistnet):
allSpecies_mist <- unique(netData.df$SPECIES)

 ### unique species (acoustic):
allSpecies_acoust <- unique(acousticData.df$Prominent.Species)

allSpecies_acoust <- NameReclass$oldname[match(allSpecies_acoust,NameReclass$newname)]

allSpecies_mist
allSpecies_acoust

allSpecies <- intersect(allSpecies_mist,allSpecies_acoust)   ## names of all species in common in the two data sets

allSpecies

#########
# Make the "Mistnet" dataset

Mistnet.df <- netEffort.df
 
  # add columns for each species of interest
for(i in 1:length(allSpecies)){
  colname <- paste(allSpecies[i],"_caps",sep="")
  expr <- paste("Mistnet.df$",colname," <- 0",sep="")
  eval(parse(text=expr))
}

Mistnet.df$DATE2 <- mdy(Mistnet.df$DATE)

head(Mistnet.df)

###############
### Summarize captures of each species by net

i=2
i=5
for(i in 1:nrow(Mistnet.df)){
  ss <- subset(netData.df,NET.CODE==Mistnet.df$NET[i])
  s=1
  for(s in 1:length(allSpecies)){
    colname <- paste(allSpecies[s],"_caps",sep="")
    colnum <- which(names(ss)=="SPECIES")
    value <- length(which(ss[,colnum]==allSpecies[s]))
    expr <- paste("Mistnet.df$",colname,"[i] <- ",value,sep="")
    eval(parse(text=expr))   
  }
}

head(Mistnet.df)


   ## add column representing all captures...
Mistnet.df$All_caps <- apply(Mistnet.df[,c(13:18)],1,sum)


### Explore final mistnet data

sum(Mistnet.df[,c(13:18)])    # total of 4012 captures...
Mistnet.df[which(Mistnet.df$YEAR==2012),] 
with(Mistnet.df,tapply(All_caps,YEAR,sum) )     # most mistnet captures were in 2006 and 2007



###############
### Make the final acoustic data frame

Acoustic.df <- acousticEffort.df
 
  # add columns for each species of interest
for(i in 1:length(allSpecies)){
  colname <- paste(allSpecies[i],"_caps",sep="")
  expr <- paste("Acoustic.df$",colname," <- 0",sep="")
  eval(parse(text=expr))
}

head(Acoustic.df)

##############
##  Add rows for the nonconsecutive surveys in 2011

ss <- subset(Acoustic.df,(YEAR==2011)&(TOTAL.NIGHTS==8))

ndx <- which(acousticData.df$CODE%in%ss$CODE)

ncdates <- sort(ymd(unique(acousticData.df[ndx,]$DATE),tz="EST"))

i=1
for(i in 1:nrow(ss)){
  df <- ss[i,]
  df$TOTAL.NIGHTS[1] <- 1
  oldcode <- df$CODE[1]
  newcode <- as.numeric(paste(oldcode,"1",sep=""))
  df$CODE[1] <- newcode
  ndx <- which((acousticData.df$CODE==oldcode)&(acousticData.df$DATE2==ncdates[1]))
  if(length(ndx)>0) acousticData.df$CODE[ndx] <- newcode
  j=2
  for(j in 2:ss$TOTAL.NIGHTS[i]){
    df <- rbind(df,ss[i,])
    df$DATE2[j] <- ncdates[j]
    df$TOTAL.NIGHTS[j] <- 1
    oldcode <- df$CODE[j]
    newcode <- as.numeric(paste(oldcode,j,sep=""))
    df$CODE[j] <- newcode
    ndx <- which((acousticData.df$CODE==oldcode)&(acousticData.df$DATE2==ncdates[j]))
    if(length(ndx)>0) acousticData.df$CODE[ndx] <- newcode
  }
  delete <- which(Acoustic.df$CODE==ss$CODE[i])
  Acoustic.df <- Acoustic.df[-delete,]
  Acoustic.df <- rbind(Acoustic.df,df)
}

 # Acoustic.df
head(Acoustic.df,30)


###############
### Remove acoustic captures outside of the main survey window...

head(acousticData.df)

acoustYears <- sort(unique(acousticData.df$YEAR))
y=acoustYears[1]
ndx <- numeric(0)
ndx2 <- numeric(0)  # survey effort
for(y in acoustYears){
   start <- mdy(paste(beginRoosting,"/",y,sep=""),tz="EST")
   end <- mdy(paste(endRoosting,"/",y,sep=""),tz="EST")
   ndx <- c(ndx,which((acousticData.df$DATE2>=start)&(acousticData.df$DATE2<=end)))
   #ndx2 <- c(ndx2,which((Acoustic.df$DATE2>=(start+days(Acoustic.df$TOTAL.NIGHTS)))&(Acoustic.df$DATE2<=end)))
}

acousticData.df <- acousticData.df[ndx,]
head(acousticData.df)

#Acoustic.df <- Acoustic.df[ndx2,]

head(Acoustic.df,20)


###############
### summarize raw acoustic captures

head(acousticData.df)

summaryAA <- as.data.frame(array(0,dim=c(length(acoustYears),length(allSpecies_acoust))))
names(summaryAA) <- allSpecies_acoust
rownames(summaryAA) <- acoustYears
s="Mlucy"
for(s in allSpecies_acoust){
  newname <- NameReclass$newname[NameReclass$oldname==s]
  ss <- subset(acousticData.df,Prominent.Species==newname)
  col <- which(allSpecies_acoust==s)
  y=2012
  for(y in acoustYears){
    row = which(acoustYears==y)
    sss = subset(ss,YEAR==y)
    summaryAA[row,col] <- nrow(sss)
  } 
}
summaryAA

apply(summaryAA,1,sum)
apply(summaryAA,2,sum)
sum(apply(summaryAA,2,sum))

################  
###  Deal with acoustic surveys that start before the summer roosting season

i=29
for(i in 1:nrow(Acoustic.df)){
   start <- mdy(paste(beginRoosting,"/",Acoustic.df$YEAR[i],sep=""),tz="EST")
   end <- mdy(paste(endRoosting,"/",Acoustic.df$YEAR[i],sep=""),tz="EST")
   if(Acoustic.df$DATE2[i]<start){
     dif <- days(start-Acoustic.df$DATE2[i])@day   # difference in days
     Acoustic.df$DATE2[i]<-start
     Acoustic.df$TOTAL.NIGHTS[i] <- Acoustic.df$TOTAL.NIGHTS[i] - dif 
   }
   if((Acoustic.df$DATE2[i]+days(Acoustic.df$TOTAL.NIGHTS[i]))>end){
     dif <- days((Acoustic.df$DATE2[i]+days(Acoustic.df$TOTAL.NIGHTS[i]))-end)@day
     Acoustic.df$TOTAL.NIGHTS[i] <- Acoustic.df$TOTAL.NIGHTS[i] - dif
   }
}

ndx <- which(Acoustic.df$TOTAL.NIGHTS<1)   # remove mics that never got into roosting season
Acoustic.df[ndx,]

Acoustic.df<-Acoustic.df[-ndx,]
head(Acoustic.df,20)

  ### redo the DATETIME variable
Acoustic.df$DATETIME <- Acoustic.df$DATE2 + hours(0)
Acoustic.df$DATE <- as.character(format(Acoustic.df$DATE2, "%Y-%m-%d"))


###############
### Summarize acoustic captures of each species by (station?)

acousticEffort.df <- Acoustic.df
head(acousticEffort.df,20)

i=7
i=1
i=11
i=26
for(i in 1:nrow(acousticEffort.df)){
  ss <- subset(acousticData.df,CODE==acousticEffort.df$CODE[i])
  ss <- ss[order(ss$DATETIME),]    # sort by date-time
  s=1
  if(acousticEffort.df$TOTAL.NIGHTS[i]<=3){
    for(s in 1:length(allSpecies)){
      colname <- paste(allSpecies[s],"_caps",sep="")
      colnum <- which(names(ss)=="Prominent.Species")
      newname <- NameReclass$newname[NameReclass$oldname==allSpecies[s]]
      sss <- subset(ss,Prominent.Species==newname)     
      if(nrow(sss)>1){
         #mindif <- round(as.numeric(sss$DATETIME[2:nrow(sss)]-sss$DATETIME[1:(nrow(sss)-1)])/60)   # should be in minutes
         mindif <- round(as.numeric(difftime(sss$DATETIME[2:nrow(sss)],sss$DATETIME[1:(nrow(sss)-1)],units="mins")),0)
      }else{
         mindif <- numeric(nrow(sss))
      }
      value <- computecaps(mindif)
      if(value>maxAbund) value <- maxAbund
          #value <- length(which(ss[,colnum]==allSpecies[s]))
      row <- which(Acoustic.df$CODE == acousticEffort.df$CODE[i])
      expr <- paste("Acoustic.df$",colname,"[row] <- ",value,sep="")
      eval(parse(text=expr))
    }
             
  }else{    ## otherwise break into 3-day segments (similar to mistnet)
    startDate <- acousticEffort.df$DATE2[i]
    totNights <-  acousticEffort.df$TOTAL.NIGHTS[i]
    endDate <- startDate + days(totNights)
    interval <- interval(startDate,endDate)
    numSegments <- interval %/% days(3)   # note that the last period may not be 3 days...
    z=1

    alldates <- ymd(ss$DATE2,tz="EST")
    flag = TRUE
    for(z in 1:(numSegments+1)){
      newrow <- Acoustic.df[i,]    #new row to be added
      newstart <- startDate+days(3*(z-1))
      newdur <- days(endDate-newstart)$day
      newrow$DATE2[1] <- newstart  #format(newstart, format="%m/%d/%Y") 
      newrow$YEAR[1] <- year(newstart)
	    flag=TRUE
      if(newdur<3){ 
	      newrow$TOTAL.NIGHTS[1] <- newdur
	      #interval(newstart,(newstart+newdur))
	      newend <- newstart + newdur
	      if(newdur==0) flag=FALSE
	    }else{ 
	      newrow$TOTAL.NIGHTS[1] <- 3
        newinterval <- interval(newstart,(newstart+days(3)))
	      newend <- newstart+days(3)
      }		    
      sss <- ss[alldates >= newstart & alldates <= newend , ]       
      s=2
      for(s in 1:length(allSpecies)){
        colname <- paste(allSpecies[s],"_caps",sep="")
        colnum <- which(names(ss)=="Prominent.Species")
        newname <- NameReclass$newname[NameReclass$oldname==allSpecies[s]]
        ssss <- subset(sss,Prominent.Species==newname)
        if(nrow(ssss)>1){
          #mindif <- round(as.numeric(ssss$DATETIME[2:nrow(ssss)]-ssss$DATETIME[1:(nrow(ssss)-1)])/60)
          mindif <- round(as.numeric(difftime(ssss$DATETIME[2:nrow(ssss)],ssss$DATETIME[1:(nrow(ssss)-1)],units="mins")),0)
        }else{
          mindif <- numeric(nrow(ssss))
        }
        value <- computecaps(mindif)
          #value <- length(which(sss[,colnum]==allSpecies[s]))
        if(value>maxAbund) value <- maxAbund
        expr <- paste("newrow$",colname,"[1] <- ",value,sep="")
        eval(parse(text=expr))        
      }
      if(flag) Acoustic.df <- rbind(Acoustic.df,newrow)
    } 
  }
}

  ### 3 warnings: failed to parse... 

Acoustic.df <- subset(Acoustic.df,TOTAL.NIGHTS<4)    # remove surveys greater than 3 days in length

   ## add column representing all captures...
ndx <- grep("_caps",names(Acoustic.df))
Acoustic.df$All_caps <- apply(Acoustic.df[,ndx],1,sum)

head(Acoustic.df,20)

### redo the DATETIME variable
Acoustic.df$DATETIME <- Acoustic.df$DATE2 + hours(0)
Acoustic.df$DATE <- as.character(format(Acoustic.df$DATE2, "%Y-%m-%d"))

###################
##### Explore final acoustic dataset...

sum(Acoustic.df[,ndx])    # total of 4759 captures...
Acoustic.df[which(Acoustic.df$YEAR==2012),]    # only 9 bouts with 94 captures in 2012- not enough to serve as a good anchor?
with(Acoustic.df,tapply(All_caps,YEAR,sum) )     # most acoustic captures were in 2011. Acoustic data not all that useful





##########################
   ##  check to see how total captures compares to final

###########################
#  ADD WEATHER DATA
###########################

#######
### Add weather data to Mistnet.df

# i=1
# Mistnet.df$precip <- NA
# Mistnet.df$tmin <- NA
# Mistnet.df$awnd <- NA
# Mistnet.df$goodcond_days <- NA
# 
# for(i in 1:nrow(Mistnet.df)){
#   startdate <- mdy(Mistnet.df$DATE[i])
#   enddate <- startdate + days(Mistnet.df$NIGHTS[i])
#   which((weatherData.df$DATE>=startdate)&(weatherData.df$DATE<=enddate))
#   ndx <- which((weatherData.df$DATE>=startdate)&(weatherData.df$DATE<=enddate))
#   prcp <- ifelse(weatherData.df$PRCP[ndx]<5,1,0)
#   tmin <- ifelse(weatherData.df$TMIN[ndx]<=10,0,weatherData.df$TMIN[ndx]-10)
#   awnd <- ifelse(weatherData.df$AWND[ndx]<3,abs(3-weatherData.df$AWND[ndx]),0)
#   Mistnet.df$goodcond_days[i] <- sum(prcp * tmin * awnd)
#   Mistnet.df$precip[i] <- mean(weatherData.df$PRCP[ndx])
#   Mistnet.df$tmin[i] <- mean(weatherData.df$TMIN[ndx])
#   Mistnet.df$awnd[i] <- mean(weatherData.df$AWND[ndx])
# }
# 
# cbind(mdy(Mistnet.df$DATE),Mistnet.df$goodcond_days,Mistnet.df$Mlucy_caps, Mistnet.df$awnd, Mistnet.df$tmin)

##########
#  READ IN WEATHER DATA AND APPEND TO MISTNET DATA FRAME
##########

airportCode <- "MHT"

# list of all days in the master dataset...


allIntervals <- interval(Mistnet.df$DATE2,Mistnet.df$DATE2+days(Mistnet.df$NIGHTS),tz="EST") 
uniqueIntervals <- unique(allIntervals)    # all periods of interest

startDays <- format(uniqueIntervals@start, "%Y-%m-%d")
endDays <- format(uniqueIntervals@start+seconds(uniqueIntervals@.Data), "%Y-%m-%d")

### download daily weather data for all relevant dates

i=1
for(i in 1:length(uniqueIntervals)){
  temp <- getSummarizedWeather(airportCode, startDays[i],endDays[i], station_type = "airportCode",
                               opt_temperature_columns = TRUE, opt_all_columns = TRUE)
  if(i==1){   # add to master data set
    dailyweather <- temp[,-2] 
  } else {
    dailyweather <- rbind(dailyweather,temp[,-2])
  }   
}

#remove non-unique observations
a <- unique(dailyweather$Date)
ndx <- match(a,dailyweather$Date)
dailyweather <- dailyweather[ndx,]

nrow(dailyweather)

# save to CSV
setwd(dataDir)
write.csv(dailyweather, file="dailyweather_mistnet.csv", row.names=FALSE)
setwd(dataDir)
dailyweather <- read.csv("dailyweather_mistnet.csv",header=T)

head(dailyweather)
dailyweather$YEAR <- year(dailyweather$Date) 
dailyweather$MONTH <- month(dailyweather$Date)
dailyweather$DATE <- format(ymd(dailyweather$Date), "%Y-%m-%d")
dailyweather$DATE2 <- ymd(dailyweather$Date,tz = "EST")


head(dailyweather)



### add weather data to the acoustic master data set (takes a while- very inefficient!)

Mistnet.df$DATE <- format(Mistnet.df$DATE2, "%Y-%m-%d")
Mistnet.df$DAY <- day(Mistnet.df$DATE2)
Mistnet.df$TempC <- 0
Mistnet.df$WindMPS <- 0
Mistnet.df$PrecipCM <- 0
#Mistnet.df$Condition <- ""

i=1
for(i in 1:nrow(Mistnet.df)){
  all <- Mistnet.df$DATE2[i] + days(1:Mistnet.df$NIGHTS[i])
  all <- format(all, "%Y-%m-%d")
  ss <- subset(dailyweather,DATE%in%all)
  Mistnet.df$TempC[i] <- FtoC(mean(ss$Mean_TemperatureF))
  Mistnet.df$WindMPS[i] <- MPHtoMPS(mean(ss$Mean_Wind_SpeedMPH))
  Mistnet.df$PrecipCM[i] <- IntoCm(mean(ss$PrecipitationIn))
  #Acoustic.df$Condition[i] <- CondToCond(as.character(ss[ndx,]$Conditions))
}


############
### Add weather data to Acoustic.df
############

# i=1
# Acoustic.df$precip <- NA
# Acoustic.df$tmin <- NA
# Acoustic.df$awnd <- NA
# Acoustic.df$goodcond_days <- NA
# 
# i=1
# for(i in 1:nrow(Acoustic.df)){
#   startdate <- Acoustic.df$DATE2[i] #ymd(Acoustic.df$DATE[i])
#   enddate <- startdate + days(Acoustic.df$TOTAL.NIGHTS[i])
#   which((weatherData.df$DATE>=startdate)&(weatherData.df$DATE<=enddate))
#   ndx <- which((weatherData.df$DATE>=startdate)&(weatherData.df$DATE<=enddate))
#   prcp <- ifelse(weatherData.df$PRCP[ndx]<5,1,0)
#   tmin <- ifelse(weatherData.df$TMIN[ndx]<=10,0,weatherData.df$TMIN[ndx]-10)
#   awnd <- ifelse(weatherData.df$AWND[ndx]<3,abs(3-weatherData.df$AWND[ndx]),0)
#   Acoustic.df$goodcond_days[i] <- sum(prcp * tmin * awnd,na.rm=T)
#   Acoustic.df$precip[i] <- mean(weatherData.df$PRCP[ndx],na.rm=T)
#   Acoustic.df$tmin[i] <- mean(weatherData.df$TMIN[ndx],na.rm=T)
#   Acoustic.df$awnd[i] <- mean(weatherData.df$AWND[ndx],na.rm=T)
# }
# 
# cbind(ymd(Acoustic.df$DATE),Acoustic.df$goodcond_days,Acoustic.df$Mlucy_caps, Acoustic.df$awnd)
# 



################

##########
#  READ IN WEATHER DATA AND APPEND TO ACOUSTICS DATA FRAME
##########

airportCode <- "MHT"

# list of all days in the master dataset...


allIntervals <- interval(Acoustic.df$DATE2,Acoustic.df$DATE2+days(Acoustic.df$TOTAL.NIGHTS),tz="EST") 
uniqueIntervals <- unique(allIntervals)    # all periods of interest

startDays <- format(uniqueIntervals@start, "%Y-%m-%d")
endDays <- format(uniqueIntervals@start+seconds(uniqueIntervals@.Data), "%Y-%m-%d")

### download daily weather data for all relevant dates

i=1
for(i in 1:length(uniqueIntervals)){
  temp <- getSummarizedWeather(airportCode, startDays[i],endDays[i], station_type = "airportCode",
                             opt_temperature_columns = TRUE, opt_all_columns = TRUE)
  if(i==1){   # add to master data set
    dailyweather <- temp[,-2] 
  } else {
    dailyweather <- rbind(dailyweather,temp[,-2])
  }   
}

   #remove non-unique observations
a <- unique(dailyweather$Date)
ndx <- match(a,dailyweather$Date)
dailyweather <- dailyweather[ndx,]

nrow(dailyweather)

# save to CSV
setwd(dataDir)
write.csv(dailyweather, file="dailyweather_acoust.csv", row.names=FALSE)
setwd(dataDir)
dailyweather <- read.csv("dailyweather_acoust.csv",header=T)

head(dailyweather)
dailyweather$YEAR <- year(dailyweather$Date) 
dailyweather$MONTH <- month(dailyweather$Date)
dailyweather$DATE <- format(ymd(dailyweather$Date), "%Y-%m-%d")
dailyweather$DATE2 <- ymd(dailyweather$Date,tz = "EST")


head(dailyweather)



### add weather data to the acoustic master data set (takes a while- very inefficient!)

Acoustic.df$DATE <- format(Acoustic.df$DATE2, "%Y-%m-%d")
Acoustic.df$DAY <- day(Acoustic.df$DATE2)
Acoustic.df$TempC <- 0
Acoustic.df$WindMPS <- 0
Acoustic.df$PrecipCM <- 0
 #Acoustic.df$Condition <- ""

i=1
for(i in 1:nrow(Acoustic.df)){
  all <- Acoustic.df$DATE2[i] + days(1:Acoustic.df$TOTAL.NIGHTS[i])
  all <- format(all, "%Y-%m-%d")
  ss <- subset(dailyweather,DATE%in%all)
  Acoustic.df$TempC[i] <- FtoC(mean(ss$Mean_TemperatureF))
  Acoustic.df$WindMPS[i] <- MPHtoMPS(mean(ss$Mean_Wind_SpeedMPH))
  Acoustic.df$PrecipCM[i] <- IntoCm(mean(ss$PrecipitationIn))
  #Acoustic.df$Condition[i] <- CondToCond(as.character(ss[ndx,]$Conditions))
}


################
           #  Make "effort" covariate

countCols <- grep("_caps",names(Mistnet.df))

Mistnet.df$effort <- Mistnet.df$NIGHTS*(Mistnet.df$NET.SIZE/18)*(Mistnet.df$NET.HEIGHT/2)   ## net-nights represented by each survey [note: net.height indicates the size of each net]

  ## NOTE: for acoustic data set, effort is just represented by number of nights.

##################################
#    VISUALIZE RELATIONSHIPS
###################################

#### 
# Captures by species

mn_allspec <- c(Mistnet.df$Mlucy_caps,Mistnet.df$Efusc_caps,
                Mistnet.df$Mleib_caps, Mistnet.df$Lcine_caps,
                Mistnet.df$Msept_caps, Mistnet.df$Lbore_caps)
mn_year <- rep(c(Mistnet.df$YEAR),times=nSpecies)
mn_specnames <- rep(allSpecies,each=nrow(Mistnet.df))

tempdf <- data.frame(counts=mn_allspec,names=mn_specnames,year=mn_year)
mistnet_byspec <- tapply(tempdf$counts,tempdf$names,sum)
mistnet_byspec_year <- tapply(tempdf$counts,list(spec=tempdf$names,year=tempdf$year),sum)


a_allspec <- c(Acoustic.df$Mlucy_caps,Acoustic.df$Efusc_caps,
                Acoustic.df$Mleib_caps, Acoustic.df$Lcine_caps,
                Acoustic.df$Msept_caps, Acoustic.df$Lbore_caps)
a_year <- rep(c(Acoustic.df$YEAR),times=nSpecies)
a_specnames <- rep(allSpecies,each=nrow(Acoustic.df))
tempdf <- data.frame(counts=a_allspec,names=a_specnames,year=a_year)
acoust_byspec <- tapply(tempdf$counts,tempdf$names,sum)
acoust_byspec_year <- tapply(tempdf$counts,list(spec=tempdf$names,year=tempdf$year),sum)

visualize = FALSE

if(visualize){

      head(Acoustic.df)
 
      andx <- which((Acoustic.df$YEAR==2011))

	plot(Acoustic.df$Mlucy_caps ~ Acoustic.df$DATE2)
	plot(Mistnet.df$Mlucy_caps ~ Mistnet.df$DATE2)

      plot(Acoustic.df$Efusc_caps ~ Acoustic.df$precip)     # lower capture rate with higher precip
	plot(Mistnet.df$Efusc_caps ~ Mistnet.df$precip)

      plot(Acoustic.df$Efusc_caps ~ Acoustic.df$tmin)     # optimal 15-17 C?  
	plot(Mistnet.df$Efusc_caps ~ Mistnet.df$tmin)

      plot(Acoustic.df$Efusc_caps ~ Acoustic.df$awnd)     # optimal 15-17 C?  
	plot(Mistnet.df$Efusc_caps ~ Mistnet.df$awnd)


	plot(Acoustic.df$Efusc_caps ~ Acoustic.df$DATE2)
	plot(Mistnet.df$Efusc_caps ~ Mistnet.df$DATE2)

	plot(Acoustic.df$Mlucy_caps ~ Acoustic.df$DATE2)
	plot(Mistnet.df$Mlucy_caps ~ Mistnet.df$DATE2)

	plot(Acoustic.df$Mleib_caps ~ Acoustic.df$DATE2)
	plot(Mistnet.df$Mleib_caps ~ Mistnet.df$DATE2)

      plot(Acoustic.df$Efusc_caps[andx] ~ Acoustic.df$DATE2[andx])
}



##################################
#      PREPARE DATA FOR BUGS    [NOW: REMOVE YEAR 2002]
##################################


#####
#  first determine the global control parameters

nSpecies <- length(allSpecies)
specRow_mist <- numeric(nSpecies)
specRow_acoust <- numeric(nSpecies)

for(s in 1:nSpecies){
  specRow_mist[s] <- grep(allSpecies[s],names(Mistnet.df))
  specRow_acoust[s]  <- grep(allSpecies[s],names(Acoustic.df))
}

allHabitats <- sort(union(as.character(unique(Mistnet.df$HABITAT)),as.character(unique(Acoustic.df$HABITAT))))
nHabitats <- length(allHabitats)
habNums <- c(1:nHabitats)

Mistnet.df$HABNUM <- match(as.character(Mistnet.df$HABITAT),allHabitats)
Acoustic.df$HABNUM <- match(as.character(Acoustic.df$HABITAT),allHabitats)


#########
##  Reclassify habitats into new categories...  

allHabitats
habNums
habNums2 <- c(1,1,2,3,2,3,2) 

Mistnet.df$HABNUM2 <- habNums2[Mistnet.df$HABNUM]
Acoustic.df$HABNUM2 <- habNums2[Acoustic.df$HABNUM]

nHabitats <- max(habNums2)

mistnetYears <- unique(Mistnet.df$YEAR)    
acoustYears <- unique(Acoustic.df$YEAR)
survYears <- sort(union(mistnetYears,acoustYears))[-1]
nSurvYears <- length(survYears)

firstYear <- min(survYears)
lastYear <- max(survYears)

mistnetVisits <- numeric(nSurvYears)
acoustVisits <- numeric(nSurvYears)
for(i in 1:nSurvYears){
  mistnetVisits[i] <- nrow(subset(Mistnet.df,YEAR==survYears[i]))
  acoustVisits[i] <- nrow(subset(Acoustic.df,YEAR==survYears[i]))
}


#################
#  Link env stochasticity among the two datasets?   [in this case not really possible, no time period where both survey methods were in effect]

### break into 6 day periods?
per_len <- 6    # days

y=6
Acoustic.df$PERIOD <- NA
Mistnet.df$PERIOD <- NA
for(y in 1:nSurvYears){
   m_temp <- subset(Mistnet.df,YEAR==survYears[y])
   a_temp <- subset(Acoustic.df,YEAR==survYears[y])
   alldates <- sort(c(m_temp$DATE2,a_temp$DATE2))
   mindate <- alldates[1]
   maxdate <- alldates[length(alldates)]
   interval <- interval(mindate,maxdate)
   period <- 1
   np <- interval%/%days(per_len)+1
   p=1
   for(p in 1:np){
     st <- mindate + days((p-1)*per_len)
     en <- mindate + days((p)*per_len)
     # newint <- interval(st,en)
     if(any((alldates>=st)&(alldates<en))){
       mdates <- which((Mistnet.df$YEAR==survYears[y])&(Mistnet.df$DATE2>=st)&(Mistnet.df$DATE2<en))
       adates <- which((Acoustic.df$YEAR==survYears[y])&(Acoustic.df$DATE2>=st)&(Acoustic.df$DATE2<en))
       if(length(mdates)>0) Mistnet.df$PERIOD[mdates] <- p
       if(length(adates)>0) Acoustic.df$PERIOD[adates] <- p
       p=p+1
     }
   }
}

cbind(Mistnet.df$DATE,Mistnet.df$PERIOD)

nperiods <- numeric(nSurvYears)
for(i in 1:nSurvYears){
  nperiods[i] <- max(c(subset(Mistnet.df,YEAR==survYears[i])$PERIOD,subset(Acoustic.df,YEAR==survYears[i])$PERIOD))
}


#####
#  set up data structures

counts.mistnet <- array(NA,dim=c(nSpecies,nSurvYears,max(mistnetVisits)))
counts.acoust <- array(NA,dim=c(nSpecies,nSurvYears,max(acoustVisits)))

netdistToPond <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
micdistToPond <- array(0,dim=c(nSurvYears,max(acoustVisits)))
netduration <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
micduration <- array(0,dim=c(nSurvYears,max(acoustVisits)))
netheight <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
netlength <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
micheight <- array(0,dim=c(nSurvYears,max(acoustVisits)))
nethabNum <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
michabNum <- array(0,dim=c(nSurvYears,max(acoustVisits)))
nettavg <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
mictavg <- array(0,dim=c(nSurvYears,max(acoustVisits)))
netawnd <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
micawnd <- array(0,dim=c(nSurvYears,max(acoustVisits)))
netprcp <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
micprcp <- array(0,dim=c(nSurvYears,max(acoustVisits)))
netperiod <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
micperiod <- array(0,dim=c(nSurvYears,max(acoustVisits)))
netmonth <- array(0,dim=c(nSurvYears,max(mistnetVisits)))
micmonth <- array(0,dim=c(nSurvYears,max(acoustVisits)))

neteffort <- array(0,dim=c(nSurvYears,max(mistnetVisits)))

  ## prepare acoustic effort
Acoustic.df$HOURSperNIGHT <- ifelse(Acoustic.df$HOURSperNIGHT>1,5,Acoustic.df$HOURSperNIGHT)

s=1;y=2;v=1
for(y in 1:nSurvYears){
  m_temp <- subset(Mistnet.df,YEAR==survYears[y])
  m_visits <- nrow(m_temp)
  a_temp <- subset(Acoustic.df,YEAR==survYears[y])
  a_visits <- nrow(a_temp)
  for(s in 1:nSpecies){
    colndx <- specRow_mist[s] 
    for(v in 1:mistnetVisits[y]){
      counts.mistnet[s,y,v] <- m_temp[v,colndx]
    }
    colndx <- specRow_acoust[s]
    for(v in 1:acoustVisits[y]){
      counts.acoust[s,y,v] <- a_temp[v,colndx]
    }
  }
}

y=1;v=1
for(y in 1:nSurvYears){
  m_temp <- subset(Mistnet.df,YEAR==survYears[y])
  m_visits <- nrow(m_temp)
  a_temp <- subset(Acoustic.df,YEAR==survYears[y])
  a_visits <- nrow(a_temp)
  if(m_visits>0){
    for(v in 1:mistnetVisits[y]){
      netdistToPond[y,v] <- m_temp$PONDm[v]
      netduration[y,v] <- m_temp$NIGHTS[v]
      netheight[y,v] <-  m_temp$NET.HEIGHT[v]
      netlength[y,v] <- m_temp$NET.SIZE[v]
      nethabNum[y,v] <- m_temp$HABNUM2[v]
      neteffort[y,v] <- m_temp$effort[v]
      nettavg[y,v] <- m_temp$TempC[v]
      netawnd[y,v] <- m_temp$WindMPS[v]
      netprcp[y,v] <- m_temp$PrecipCM[v]
      netperiod[y,v] <- m_temp$PERIOD[v]
      netmonth[y,v] <- month(m_temp$DATE2[v])
    }
  }
  if(a_visits>0){
    for(v in 1:acoustVisits[y]){
      micdistToPond[y,v] <- a_temp$PONDm[v] 
      micduration[y,v] <- (a_temp$HOURSperNIGHT[v]/5)*a_temp$TOTAL.NIGHTS[v]   ### effectively, total # hours surveyed
          #micheight[y,v] <-  a_temp$
      michabNum[y,v] <- a_temp$HABNUM2[v]
      mictavg[y,v] <- a_temp$TempC[v]
      micawnd[y,v] <- a_temp$WindMPS[v]
      micprcp[y,v] <- a_temp$PrecipCM[v]
      micperiod[y,v] <- a_temp$PERIOD[v]
      micmonth[y,v] <- month(a_temp$DATE2[v])
    }
  }
}

allmonths <- sort(union(unique(month(Mistnet.df$DATE2)),unique(month(Acoustic.df$DATE2)) ))
nmonths <- length(allmonths)

#################################
# NORMALIZE COVARIATES (preprocessing for BUGS)
#################################

normalizeCovariate <- function(t){
  t_na <- apply(t,c(1,2),function(z) ifelse(z==0,NA,z))
  t2 <- apply(t_na,c(1,2),function(z) (z-mean(t_na,na.rm=T))/sd(as.vector(t_na),na.rm=T))
  r <- apply(t2,c(1,2),function(z) ifelse(is.na(z),0,z))
  return(r)
}

netdistToPond_std <- normalizeCovariate(netdistToPond)
micdistToPond_std <- normalizeCovariate(micdistToPond)
#netduration_std <- normalizeCovariate(netduration)
#micduration_std <- normalizeCovariate(micduration)
#netheight_std <- normalizeCovariate(netheight)
#netlength_std <- normalizeCovariate(netlength)
#micheight_std <- normalizeCovariate(micheight)

mictavg_std <- normalizeCovariate(mictavg)
micawnd_std <- normalizeCovariate(micawnd)
micprcp_std <- normalizeCovariate(micprcp)
nettavg_std <- normalizeCovariate(nettavg)
netawnd_std <- normalizeCovariate(netawnd)
netprcp_std <- normalizeCovariate(netprcp)


####  change the scale of the acoustic data

#counts.acoust2 <- apply(counts.acoust,c(1,2,3),function(t) ifelse(is.na(t),NA,
#        ifelse(t==0,0,ceil(log(t+0.1)) ) ))

#counts.acoust <- counts.acoust2


#################################
#     SET BUGS OR JAGS DIRECTORY
##################################

if(KEVIN){
  bugs.directory="C:\\Users\\Kevin\\Documents\\Employment\\ESF\\Bog Turtle\\DATA\\software\\BUGS\\WinBUGS14"
  working.directory= BUGSDir 
  setwd(working.directory)
}

if(OFFICE){
  bugs.directory="C:\\WinBUGS\\winbugs14\\WinBUGS14"
  working.directory= BUGSDir 
  setwd(working.directory)
}

############################################################
#      # WRITE BUGS OR JAGS MODEL TO FILE: ALL MODELS IN ONE
############################################################

filename_ALL <- "NewBoston9_ALL_2012and2013"    
model.file <- paste(filename_ALL,".bug",sep="")

sink(model.file)

cat("

###################################
#    START BUGS MODEL
###################################

    model {

	
##################################
      # DEFINE DATA  
##################################	
	
	for(s in 1:ns){                                                                              # loop through species
		for (y in 1:ny){                                                                         # loop through years
			for (v in 1:nnv[y]){                                                                 # loop through mist nets
				for(k in 1:3){                                                                   # loop through models
				   counts.mistnet2[s,y,v,k] <- counts.mistnet[s,y,v]
				}
				  #counts.mistnet.BEST[s,y,v] <- counts.mistnet[s,y,v]
			}
			for (v in 1:anv[y]){     # acoustic
				for(k in 1:3){
				   counts.acoust2[s,y,v,k] <- counts.acoust[s,y,v]
				}
				  #counts.acoust.BEST[s,y,v] <- counts.acoust[s,y,v]
			}			
		}
	}
			
##################################
      # SET DATA LIKELIHOOD  
##################################
 
	## OBSERVATION PROCESS #1: mistnet     

	for(s in 1:ns){                                                                              ## loop through species
		for (y in 1:ny){                                                                         ## loop through years
			for (v in 1:nnv[y]){                                                                 ## loop through nets (for each year)  

					logit.pcap.mistnet[s,y,v] <- logit.p0.mistnet  
					 + tempeff * nettavg[y,v]            ## probability of being captured by mistnet, standard effort
					 + precipeff * netprcp[y,v]
					 + windeff * netawnd[y,v]

					pcap.mistnet[s,y,v] <- (1/(1+exp(-1*(logit.pcap.mistnet[s,y,v]))))                 
																		 ## precision of the observation: if not occupied, this precision should be quite high!                            
				for(k in 1:3){                                                                   ## loop through models																		 ## OBSERVATION MODEL: 
					counts.mistnet.exp[s,y,v,k] <- (abund.total[s,y,k] * pcap.mistnet[s,y,v]) * neteffort[y,v] * (1-netzero[y,v])  
					counts.mistnet2[s,y,v,k] ~ dpois(counts.mistnet.exp[s,y,v,k])                ## DATA NODE
				}
			      
        #counts.mistnet.BEST.exp[s,y,v] <- (abund.total.BEST[s,y] * pcap.mistnet[s,y,v]) * neteffort[y,v] * (1-netzero[y,v])
        #counts.mistnet.BEST[s,y,v] ~ dpois(counts.mistnet.BEST.exp[s,y,v])                ## MODEL SELECTION  - BEST MODEL
			
			}
		}
	}	

## OBSERVATION PROCESS #2: acoustics     
  
	for(s in 1:ns){                                                                              ## loop through species
		for (y in 1:ny){                                                                         ## loop through years
			for (v in 1:anv[y]){                                                                 ## loop through nets (for each year)   
				
				logit.pcap.acoust[s,y,v] <- logit.p0.acoust 
					+ tempeff * mictavg[y,v]                ## probability of being captured by microphone
					+ treebateff * acoust.det[s]            ## treebats are less detectable
					+ precipeff * micprcp[y,v]
					+ windeff * micawnd[y,v]
				
				pcap.acoust[s,y,v] <- (1/(1+exp(-1*(logit.pcap.acoust[s,y,v])))) 
				for(k in 1:3){
																	 ## OBSERVATION MODEL:  
					counts.acoust.exp[s,y,v,k] <- (abund.total[s,y,k] * pcap.acoust[s,y,v]) * miceffort[y,v] * (1-miczero[y,v])
					counts.acoust2[s,y,v,k] ~ dpois(counts.acoust.exp[s,y,v,k])     # DATA NODE
				}

				#counts.acoust.BEST.exp[s,y,v] <- (abund.total.BEST[s,y] * pcap.acoust[s,y,v]) * miceffort[y,v] * (1-miczero[y,v]) 
				#counts.acoust.BEST[s,y,v] ~ dpois(counts.acoust.BEST.exp[s,y,v])                ## MODEL SELECTION  - BEST MODEL				
			}
		}
	}
   
###########################   
      # SET PRIORS
###########################

	#########################################################################
	########################## BAYESIAN MODEL SELECTION
	##########################################################################								
									# goal: select the process that is most likely to have generated the observed data
									
											       # assign prior model weights
# 	pM[1] <- 10/30     # fully time-dependent        # assign 50% of prior model weight to the null model to ensure that this model is selected for non-declining time series
# 	pM[2] <- 10/30  # trend                       # assign equal model weight to the habitat loss and harvest scenarios.
# 	pM[3] <- 10/30  # WNS
# 
# 	for(s in 1:ns){
# 		mod[s] ~ dcat(pM[1:3])     #  MODEL SELECTION VARIABLE: indicates which candidate model is most likely to have generated the observed data
# 	}	

	#########################################
	#  RELATIVE ABUNDANCE
	#########################################

  ## Set categorical abundance variable...
  
	incr <- (upperBoundLog-lowerBoundLog)/npos
	for(i in 1:npos){
		priorNp[i] <- 1/npos                         # probability of each possible abundance value is equivalent- uniformly distributed on log scale
		NvalLog[i] <- lowerBoundLog + incr*(i-1)
		Nval[i] <- trunc(exp(NvalLog[i]))            # all possible N values      
	}

	 ## to set 2011/2012 eFuscus baseline (set at 1000 arbitrarily)
	for(i in 1:npos){
		Np2[i] <- 1/npos
		Nval2[i] <- abund.baseline
	}

		############
		# MODEL 1: fully time-dependent (annual abundance)
		#############
		for(s in 1:1){
			for(y in 1:ny){
				abund.intX[s,y,1] ~ dcat(priorNp[])
				abund.total[s,y,1] <- Nval[abund.intX[s,y,1]]   # discrete N value 
			}
		}
		for(s in 2:2){     # single out E. fuscus
			for(y in 1:3){
				abund.intX[s,y,1] ~ dcat(priorNp[])
				abund.total[s,y,1] <- Nval[abund.intX[s,y,1]]   # discrete N value
			}
			abund.intX[s,4,1] ~ dcat(priorNp[]) #dcat(Np2[])
			abund.total[s,4,1] <- Nval[abund.intX[s,4,1]]
			abund.intX[s,5,1] ~ dcat(Np2[]) #dcat(Np2[])
			abund.total[s,5,1] <- Nval2[abund.intX[s,5,1]]    # use years 2012 and 2013, big brown bat, as the reference year/species...
			abund.intX[s,6,1] ~ dcat(Np2[]) #dcat(priorNp[])
			abund.total[s,6,1] <- Nval2[abund.intX[s,6,1]]
			abund.intX[s,ny,1] ~ dcat(priorNp[])
			abund.total[s,ny,1] <- Nval[abund.intX[s,ny,1]]
		}
		for(s in 3:ns){
			for(y in 1:ny){
				abund.intX[s,y,1] ~ dcat(priorNp[])
				abund.total[s,y,1] <- Nval[abund.intX[s,y,1]]   # discrete N value
			}
		}
		
		############
		# MODEL 2: Linear trend   (intercept of abundance/time relationship)
		#############
		
		for(s in 1:1){
			abund.intX[s,1,2] ~ dcat(priorNp[])
			abund.int[s,1,2] <- Nval[abund.intX[s,1,2]]   # discrete N value 
		}  
		abund.intX[2,1,2] ~ dcat(Np2[])                   # set e fuscus intercept to 1000
		abund.int[2,1,2] <- Nval2[abund.intX[2,1,2]]
		for(s in 3:ns){ 
			abund.intX[s,1,2] ~ dcat(priorNp[])
			abund.int[s,1,2] <- Nval[abund.intX[s,1,2]]   # discrete N value 
		}

		for(s in 1:ns){                               ## loop through species    [try log scale?]
			log.abund.int[s] <- log(abund.int[s,1,2])
			for (y in 1:ny){            
				log.abund[s,y] <- log.abund.int[s] + trend[s] * (survYears[y]-2012) 
				abund.total[s,y,2] <- exp(log.abund[s,y])
			}
		}
		
		############
		# MODEL 3: WNS (abundance estimated separately pre- and post- WNS)
		#############		
		for(s in 1:1){
			for(j in 1:2){
				abund.intX[s,j,3] ~ dcat(priorNp[])
				abund.int[s,j,3] <- Nval[abund.intX[s,j,3]]   # discrete N value
			} 
		}  
		abund.intX[2,1,3] ~ dcat(priorNp[])                     # fix post-WNS abundance of e.fuscus
		abund.int[2,1,3] <- Nval[abund.intX[2,1,3]]
		abund.intX[2,2,3] ~ dcat(Np2[])
		abund.int[2,2,3] <- Nval2[abund.intX[2,2,3]]
		for(s in 3:ns){ 
			for(j in 1:2){
				abund.intX[s,j,3] ~ dcat(priorNp[])
				abund.int[s,j,3] <- Nval[abund.intX[s,j,3]]   # discrete N value
			}
		}

		for(y in 1:ny){
			WNSndx[y] <- step(y-3)+1    # index- is WNS in effect (2) or not in effect (1)
		}

		for(s in 1:ns){                                                                              ## loop through species
			for (y in 1:ny){
				abund.total[s,y,3] <- abund.int[s,WNSndx[y],3] 
			}
		}

        #############  MODEL AVERAGED RELATIVE ABUNDANCE

# 	for(s in 1:ns){                                                                              ## loop through species
# 		for (y in 1:ny){
# 			abund.total.BEST[s,y] <- abund.total[s,y,mod[s]]           # MODEL-AVERAGED ESTIMATE OF TOTAL ABUNDANCE...  
# 		}
# 	}

	#################################################
    ## ACCOMMODATE ZERO INFLATION: # don't let zeros contribute too much to the likelihood
	#################################################
           
	pZero.int[1] ~ dbeta(1,1)  # mistnet
	pZero.int[2] ~ dbeta(1,1)  # acoust
	pZero.int.logit[1] <- log(pZero.int[1]/(1-pZero.int[1]))
	pZero.int.logit[2] <- log(pZero.int[2]/(1-pZero.int[2]))
	for(y in 1:ny){
		for(v in 1:nnv[y]){      # mistnet
			pZero.exp.logit[y,v,1] <- pZero.int.logit[1] + tempeff2*nettavg[y,v] + montheff2[netmonth[y,v]-4] 
                  pZero.exp[y,v,1] <- (1/(1+exp(-1*(pZero.exp.logit[y,v,1]))))
			netzero[y,v] ~ dbern(pZero.exp[y,v,1])   # probability of zero captures
		}
		for(v in 1:anv[y]){     # acoust
			pZero.exp.logit[y,v,2] <- pZero.int.logit[2] + tempeff2*mictavg[y,v] + montheff2[micmonth[y,v]-4] 
                  pZero.exp[y,v,2] <- (1/(1+exp(-1*(pZero.exp.logit[y,v,2]))))
			miczero[y,v] ~ dbern(pZero.exp[y,v,2])   # probability of zero captures
		}   
	}
	
	tempeff2 ~ dnorm(0,1)   # temperature effect on zero captures?
	
	for(m in 1:nmonths){
	  montheff2[m] ~ dnorm(0,1)
	}

	#################################################
    ## BASE PROBABILITY OF CAPTURE
	#################################################	
	  
	## PROBCAP TERMS
	p0.mistnet ~ dbeta(1,1)  # start with probability weight concentrated near 0
	logit.p0.mistnet <- log(p0.mistnet/(1-p0.mistnet))
	p0.acoust ~ dbeta(1,1)
	logit.p0.acoust <- log(p0.acoust/(1-p0.acoust))

	
	#################################################
    ## HABITAT AND WEATHER EFFECTS    ** logit scale **
	#################################################
	

    ##  BETA TERMS (weather)  
	tempeff ~ dnorm(0,1)
	treebateff ~ dunif(1,10)
	precipeff ~ dnorm(0,1)
	windeff ~ dnorm(0,1)
	#watereff ~ dnorm(0,1)

    #########
    ## TREND  (real scale)
    #########

     for(s in 1:ns){
       trend[s] ~ dnorm(0,10)   # set temporal trend  [RECORD THIS FOR PLOTTING!]
     }
	
     ################## 
     ## WNS EFFECT
     ##################
     for(s in 1:ns){
       WNSEff[s] <- abund.int[s,2,3] - abund.int[s,1,3]    # [RECORD WNS Effect (change in abundance) FOR PLOTTING!]
     }

	############################
	#  CONSTANTS
	############################
	
	abund.baseline <- 1000                  ## initial abundance: only relative abundance is important, so set an an arbitrary (but high) value
		#log.abund.baseline <- log(1000)

	lowerBoundLog <- 4.6   #4    #4.5        # bounds on abundance/intercept variables
	upperBoundLog <- 9.5  #10   #13.8
	
			 ## acoustic detectability:  in order: Mlucy Efusc Mleib Lcine Msept Lbore
	acoust.det[1] <- 0
	acoust.det[2] <- 0
	acoust.det[3] <- 0
	acoust.det[4] <- 1   # tree bats are very differently detectable
	acoust.det[5] <- 0
	acoust.det[6] <- 1
	 
###################################
#    END BUGS MODEL
###################################
    }  
###################################
    
    ",fill = TRUE)

sink()                   ## Writes BUGS code to the appropriate file



############################################################
#      # WRITE BUGS OR JAGS MODEL TO FILE: NO ACOUSTIC DATA
############################################################

filename_noAcoust <- "NewBoston8_noAcoust"    
model.file <- paste(filename_noAcoust,".bug",sep="")

sink(model.file)

cat("

###################################
#    START BUGS MODEL
###################################

    model {

	
##################################
      # DEFINE DATA  
##################################	
	
	for(s in 1:ns){                                                                              # loop through species
		for (y in 1:ny){                                                                         # loop through years
			for (v in 1:nnv[y]){                                                                 # loop through mist nets
				for(k in 1:3){                                                                   # loop through models
				   counts.mistnet2[s,y,v,k] <- counts.mistnet[s,y,v]
				}
				#counts.mistnet.BEST[s,y,v] <- counts.mistnet[s,y,v]
			}		
		}
	}
			
##################################
      # SET DATA LIKELIHOOD  
##################################
 
	## OBSERVATION PROCESS #1: mistnet     

	for(s in 1:ns){                                                                              ## loop through species
		for (y in 1:ny){                                                                         ## loop through years
			for (v in 1:nnv[y]){                                                                 ## loop through nets (for each year)  

					logit.pcap.mistnet[s,y,v] <- logit.p0.mistnet  
					 + tempeff * nettavg[y,v]            ## probability of being captured by mistnet, standard effort
					 + precipeff * netprcp[y,v]
					 + windeff * netawnd[y,v]
					 #+ watereff * netdistToPond[y,v]                                          
					 #+ habeff[s,nethabitat[y,v]]
					 #+ condeff[s,y,netperiod[y,v],k] 

					pcap.mistnet[s,y,v] <- (1/(1+exp(-1*(logit.pcap.mistnet[s,y,v]))))                 
																		 ## precision of the observation: if not occupied, this precision should be quite high!                            
				for(k in 1:3){                                                                   ## loop through models																		 ## OBSERVATION MODEL: 
					counts.mistnet.exp[s,y,v,k] <- (abund.total[s,y,k] * pcap.mistnet[s,y,v]) * neteffort[y,v] * (1-netzero[y,v])  
					counts.mistnet2[s,y,v,k] ~ dpois(counts.mistnet.exp[s,y,v,k])                ## DATA NODE
				}
			      
                counts.mistnet.BEST.exp[s,y,v] <- (abund.total.BEST[s,y] * pcap.mistnet[s,y,v]) * neteffort[y,v] * (1-netzero[y,v])
				counts.mistnet.BEST[s,y,v] ~ dpois(counts.mistnet.BEST.exp[s,y,v])                ## MODEL SELECTION  - BEST MODEL
			
			}
		}
	}	


###########################   
      # SET PRIORS
###########################

	#########################################################################
	########################## BAYESIAN MODEL SELECTION
	##########################################################################								
									# goal: select the process that is most likely to have generated the observed data
									
											       # assign prior model weights
	pM[1] <- 10/30     # fully time-dependent        # assign 50% of prior model weight to the null model to ensure that this model is selected for non-declining time series
	pM[2] <- 10/30  # trend                       # assign equal model weight to the habitat loss and harvest scenarios.
	pM[3] <- 10/30  # WNS

	for(s in 1:ns){
		mod[s] ~ dcat(pM[1:3])     #  MODEL SELECTION VARIABLE: indicates which candidate model is most likely to have generated the observed data
	}	

	#########################################
	#  RELATIVE ABUNDANCE
	#########################################

  ## Set categorical abundance variable...
  
	incr <- (upperBoundLog-lowerBoundLog)/npos
	for(i in 1:npos){
		priorNp[i] <- 1/npos                         # probability of each possible abundance value is equivalent- uniformly distributed on log scale
		NvalLog[i] <- lowerBoundLog + incr*(i-1)
		Nval[i] <- trunc(exp(NvalLog[i]))            # all possible N values      
	}

	 ## to set 2011/2012 eFuscus baseline (set at 1000 arbitrarily)
	for(i in 1:npos){
		Np2[i] <- 1/npos
		Nval2[i] <- abund.baseline
	}

		############
		# MODEL 1: fully time-dependent (annual abundance)
		#############
		for(s in 1:1){
			for(y in 1:ny){
				abund.intX[s,y,1] ~ dcat(priorNp[])
				abund.total[s,y,1] <- Nval[abund.intX[s,y,1]]   # discrete N value 
			}
		}
		for(s in 2:2){     # single out E. fuscus
			for(y in 1:2){
				abund.intX[s,y,1] ~ dcat(priorNp[])
				abund.total[s,y,1] <- Nval[abund.intX[s,y,1]]   # discrete N value
			}
			abund.intX[s,3,1] ~ dcat(Np2[])
			abund.total[s,3,1] <- Nval2[abund.intX[s,4,1]]    # use year 2011 and 2012, big brown bat, as the reference year/species...
			abund.intX[s,ny,1] ~ dcat(priorNp[])
			abund.total[s,ny,1] <- Nval[abund.intX[s,ny,1]]
		}
		for(s in 3:ns){
			for(y in 1:ny){
				abund.intX[s,y,1] ~ dcat(priorNp[])
				abund.total[s,y,1] <- Nval[abund.intX[s,y,1]]   # discrete N value
			}
		}
		
		############
		# MODEL 2: Linear trend   (intercept of abundance/time relationship)
		#############
		
		for(s in 1:1){
			abund.intX[s,1,2] ~ dcat(priorNp[])
			abund.int[s,1,2] <- Nval[abund.intX[s,1,2]]   # discrete N value 
		}  
		abund.intX[2,1,2] ~ dcat(Np2[])                   # set e fuscus intercept to 1000
		abund.int[2,1,2] <- Nval2[abund.intX[2,1,2]]
		for(s in 3:ns){ 
			abund.intX[s,1,2] ~ dcat(priorNp[])
			abund.int[s,1,2] <- Nval[abund.intX[s,1,2]]   # discrete N value 
		}

		for(s in 1:ns){                               ## loop through species    [try log scale?]
			log.abund.int[s] <- log(abund.int[s,1,2])
			for (y in 1:ny){            
				log.abund[s,y] <- log.abund.int[s] + trend[s] * (survYears[y]-2012) 
				abund.total[s,y,2] <- exp(log.abund[s,y])
			}
		}
		
		############
		# MODEL 3: WNS (abundance estimated separately pre- and post- WNS)
		#############		
		for(s in 1:1){
			for(j in 1:2){
				abund.intX[s,j,3] ~ dcat(priorNp[])
				abund.int[s,j,3] <- Nval[abund.intX[s,j,3]]   # discrete N value
			} 
		}  
		abund.intX[2,1,3] ~ dcat(priorNp[])                     # fix post-WNS abundance of e.fuscus
		abund.int[2,1,3] <- Nval[abund.intX[2,1,3]]
		abund.intX[2,2,3] ~ dcat(Np2[])
		abund.int[2,2,3] <- Nval2[abund.intX[2,2,3]]
		for(s in 3:ns){ 
			for(j in 1:2){
				abund.intX[s,j,3] ~ dcat(priorNp[])
				abund.int[s,j,3] <- Nval[abund.intX[s,j,3]]   # discrete N value
			}
		}

		for(y in 1:ny){
			WNSndx[y] <- step(y-3)+1    # index- is WNS in effect (2) or not in effect (1)
		}

		for(s in 1:ns){                                                                              ## loop through species
			for (y in 1:ny){
				abund.total[s,y,3] <- abund.int[s,WNSndx[y],3] 
			}
		}

        #############  MODEL AVERAGED RELATIVE ABUNDANCE

	for(s in 1:ns){                                                                              ## loop through species
		for (y in 1:ny){
			abund.total.BEST[s,y] <- abund.total[s,y,mod[s]]           # MODEL-AVERAGED ESTIMATE OF TOTAL ABUNDANCE...  
		}
	}
		
	
	#################################################
    ## GOODCOND VS BADCOND  [not explained by other factors]  **logit scale**
	#################################################
	
	#for(s in 1:ns){
	#	for(k in 1:3){     # loop through models
	#		for(y in 1:ny){
	#			for(p in 1:np[y]){
	#				condeff[s,y,p,k] ~ dnorm(0,hyperprec[k])           #dbern(pbadcond)  
	#			}
	#		}
	#	}
	#}
	
	#for(k in 1:3){
	#  hyperprec[k] ~ dgamma(0.1,0.1)     # set precision of random effect on capture probability...
	#  hypersd[k] <- pow(1/hyperprec[k],0.5)
	#}

	#################################################
    ## ACCOMMODATE ZERO INFLATION: # don't let zeros contribute too much to the likelihood
	#################################################
           
	pZero.int[1] ~ dbeta(1,1)  # mistnet
	pZero.int[2] ~ dbeta(1,1)  # acoust
	pZero.int.logit[1] <- log(pZero.int[1]/(1-pZero.int[1]))
	pZero.int.logit[2] <- log(pZero.int[2]/(1-pZero.int[2]))
	for(y in 1:ny){
		for(v in 1:nnv[y]){      # mistnet
			pZero.exp.logit[y,v,1] <- pZero.int.logit[1] + tempeff2*nettavg[y,v] + montheff2[netmonth[y,v]-4] 
                  pZero.exp[y,v,1] <- (1/(1+exp(-1*(pZero.exp.logit[y,v,1]))))
			netzero[y,v] ~ dbern(pZero.exp[y,v,1])   # probability of zero captures
		}
	}
	
	tempeff2 ~ dnorm(0,1)   # temperature effect on zero captures?
	
	for(m in 1:nmonths){
	  montheff2[m] ~ dnorm(0,1)
	}

	#################################################
    ## BASE PROBABILITY OF CAPTURE
	#################################################	
	  
	## PROBCAP TERMS
	p0.mistnet ~ dbeta(1,1)  # start with probability weight concentrated near 0
	logit.p0.mistnet <- log(p0.mistnet/(1-p0.mistnet))

	
	#################################################
    ## HABITAT AND WEATHER EFFECTS    ** logit scale **
	#################################################
	
	#for(s in 1:ns){
	#    habeff[s,1] <- 0
	#	for(h in 2:nhabtypes){
	#		habeff[s,h] ~ dnorm(0,1) 
	#	}
	#}

    ##  BETA TERMS (weather)  
	tempeff ~ dnorm(0,1)
	treebateff ~ dunif(1,10)
	precipeff ~ dnorm(0,1)
	windeff ~ dnorm(0,1)
	#watereff ~ dnorm(0,1)

    #########
    ## TREND  (real scale)
    #########

     for(s in 1:ns){
       trend[s] ~ dnorm(0,10)   # set temporal trend  [RECORD THIS FOR PLOTTING!]
     }
	
     ################## 
     ## WNS EFFECT
     ##################
     for(s in 1:ns){
       WNSEff[s] <- abund.int[s,2,3] - abund.int[s,1,3]    # [RECORD WNS Effect (change in abundance) FOR PLOTTING!]
     }

	############################
	#  CONSTANTS
	############################
	
	abund.baseline <- 1000                  ## initial abundance: only relative abundance is important, so set an an arbitrary (but high) value
		#log.abund.baseline <- log(1000)

	lowerBoundLog <- 4.6   #4    #4.5        # bounds on abundance/intercept variables
	upperBoundLog <- 9.5  #10   #13.8
	
	 
###################################
#    END BUGS MODEL
###################################
    }  
###################################
    
    ",fill = TRUE)

sink()                   ## Writes BUGS code to the appropriate file


#################################
#  BUNDLE DATA FOR WINBUGS OR JAGS
#################################

Data <- list(
  ns=nSpecies,
  ny=nSurvYears,
  survYears=survYears,
  nmonths=nmonths,
  nnv=mistnetVisits,
  anv=acoustVisits,
  npos = 50,
  counts.mistnet=counts.mistnet,
  counts.acoust=counts.acoust, 
  #netdistToPond=netdistToPond_std,
  #micdistToPond=micdistToPond_std,
  nettavg = nettavg_std,
  mictavg = mictavg_std,
  netprcp = netprcp_std,
  micprcp = micprcp_std,
  netawnd = netawnd_std,
  micawnd = micawnd_std,
  miceffort = micduration/3,
  neteffort = neteffort,
  #micperiod = micperiod,
  #netperiod = netperiod,
  micmonth=micmonth,
  netmonth=netmonth 
  #nethabitat=nethabNum,
  #michabitat=michabNum
)

 # Data

#########################
# EXPERIMENT WITH INITIAL VALUES
#########################

arb_initabund <- 1000

p0_init_mistnet <- apply(counts.mistnet,c(1),function(t) mean(t,na.rm=T)/arb_initabund)
p0_init_acoustic <- apply(counts.acoust,c(1),function(t) mean(t,na.rm=T)/arb_initabund)

#p0_init_mistnet <- mean(counts.mistnet,na.rm=T)/arb_initabund*100
#p0_init_acoustic <- mean(counts.acoust,na.rm=T)/arb_initabund

    # pZero parameter
t1 <- apply(counts.mistnet,3,function(t) sum(ifelse(as.vector(t)>0,1,0),na.rm=T))
t2 <- apply(counts.mistnet,3,function(t) sum(ifelse(!is.na(t),1,0)))
pZero_init_mistnet <- mean(t1[1:10]/t2[1:10])

t1 <- apply(counts.acoust,3,function(t) sum(ifelse(as.vector(t)>0,1,0),na.rm=T))
t2 <- apply(counts.acoust,3,function(t) sum(ifelse(!is.na(t),1,0)))
pZero_init_acoustic <- mean(t1[1:20]/t2[1:20])

     


########################
# SET INITIAL VALUES
#########################

Inits <- function(){
  list(
    p0.mistnet=0.005,  
    p0.acoust=0.005, 
    tempeff = -0.4,   
    tempeff2 = 1,
    precipeff = -0.5,
    windeff = 0,
    montheff2 = rep(0,times=nmonths),
    #habeff = array(0,dim=c(nSpecies,nHabitats)),
    #condeff = array(0,dim=c(nSpecies,nSurvYears,nperiods,3)),
    #mod=rep(1,times=nSpecies),
    trend= c(-0.2,0.16,-0.05,0.01,0,0),   #rep(0,times=nSpecies),
    pZero.int = c(0.2,0.2),
    netzero = array(0,dim=c(nSurvYears,max(mistnetVisits))),
    miczero = array(0,dim=c(nSurvYears,max(acoustVisits))),
    treebateff = 2,
    #watereff = 0,
    #abund.totalX = array(30,dim=c(nSpecies,nSurvYears,1)),
    abund.intX = array(20,dim=c(nSpecies,nSurvYears,3))
  )
}
  
Inits()

############################
#  SET MONITORED PARAMS
############################

Pars <- c(
  "abund.total",
  #"abund.total.BEST",
  #"mod",
  "p0.mistnet",
  "p0.acoust",
  #"habeff",
  "tempeff",
  "trend",
  "WNSEff",
  "tempeff2",
  "montheff2",
  "treebateff",
  "windeff",   
  "pZero.int",
  "precipeff"
  #"watereff"
  #"hypersd"
  #"condeff"
)


###########################
##  NON-PARALLEL BUGS CALL, FOR DEBUGGING   [for now, this is the main way I'm running winBUGS- can't seem to get the parallel processing running]
###########################

nt=10; ni = 15000; nb = 5000; nc = 1                                            # perfect convergence! 3/2/2016
model.file <- paste(filename_ALL,".bug",sep="")

Mod <- bugs(   inits=Inits, n.chains=nc, model.file=model.file, bugs.directory=bugs.directory, 
               data=Data, parameters.to.save=Pars, 
               n.thin=nt, n.iter=ni, n.burnin=nb, DIC=FALSE, debug=T,over.relax = FALSE)
	
m <- as.mcmc.bugs(Mod)         # convert BUGS output to CODA object (from "emdbook" package)

setwd(BUGSDir)
save(list = ls(all.names = TRUE),file = "Results_2012and2013.RData", envir = .GlobalEnv)
?save

#################################
#  BUNDLE DATA FOR WINBUGS OR JAGS: no acoustic 
#################################

Data <- list(
  ns=nSpecies,
  ny=length(mistnetYears[-1]),
  survYears=mistnetYears[-1],
  nmonths=nmonths,
  nnv=mistnetVisits[-c(3,4)],
  #anv=acoustVisits,
  npos = 50,
  #np = nperiods,
  #nhabtypes=nHabitats,
  counts.mistnet=counts.mistnet[,-c(3,4),],
  #counts.acoust=counts.acoust, 
  #netdistToPond=netdistToPond_std,
  #micdistToPond=micdistToPond_std,
  nettavg = nettavg_std[-c(3,4),],
  #mictavg = mictavg_std,
  netprcp = netprcp_std[-c(3,4),],
  #micprcp = micprcp_std,
  netawnd = netawnd_std[-c(3,4),],
  #micawnd = micawnd_std,
  #miceffort = micduration/3,
  neteffort = neteffort[-c(3,4),],
  #micperiod = micperiod,
  #netperiod = netperiod,
  #micmonth=micmonth,
  netmonth=netmonth[-c(3,4),] 
  #nethabitat=nethabNum,
  #michabitat=michabNum
)

 # Data

#########################
# EXPERIMENT WITH INITIAL VALUES
#########################

arb_initabund <- 1000

p0_init_mistnet <- apply(counts.mistnet,c(1),function(t) mean(t,na.rm=T)/arb_initabund)
p0_init_acoustic <- apply(counts.acoust,c(1),function(t) mean(t,na.rm=T)/arb_initabund)

#p0_init_mistnet <- mean(counts.mistnet,na.rm=T)/arb_initabund*100
#p0_init_acoustic <- mean(counts.acoust,na.rm=T)/arb_initabund

    # pZero parameter
t1 <- apply(counts.mistnet,3,function(t) sum(ifelse(as.vector(t)>0,1,0),na.rm=T))
t2 <- apply(counts.mistnet,3,function(t) sum(ifelse(!is.na(t),1,0)))
pZero_init_mistnet <- mean(t1[1:10]/t2[1:10])

t1 <- apply(counts.acoust,3,function(t) sum(ifelse(as.vector(t)>0,1,0),na.rm=T))
t2 <- apply(counts.acoust,3,function(t) sum(ifelse(!is.na(t),1,0)))
pZero_init_acoustic <- mean(t1[1:20]/t2[1:20])

     


########################
# SET INITIAL VALUES
#########################

Inits <- function(){
  list(
    p0.mistnet=0.005,  
    #p0.acoust=0.005, 
    tempeff = -0.4,   
    tempeff2 = 1,
    precipeff = -0.5,
    windeff = 0,
    montheff2 = rep(0,times=nmonths),
    #habeff = array(0,dim=c(nSpecies,nHabitats)),
    #condeff = array(0,dim=c(nSpecies,nSurvYears,nperiods,3)),
    mod=rep(1,times=nSpecies),
    trend= c(-0.2,0.16,-0.05,0.01,0,0),   #rep(0,times=nSpecies),
    pZero.int = c(0.2,0.2),
    netzero = array(0,dim=c(nSurvYears,max(mistnetVisits))),
    #miczero = array(0,dim=c(nSurvYears,max(acoustVisits))),
    treebateff = 2,
    #watereff = 0,
    #abund.totalX = array(30,dim=c(nSpecies,nSurvYears,1)),
    abund.intX = array(20,dim=c(nSpecies,nSurvYears,3))
  )
}
  
Inits()

############################
#  SET MONITORED PARAMS
############################

Pars <- c(
  "abund.total",
  #"abund.total.BEST",
  "mod",
  "p0.mistnet",
  #"p0.acoust",
  #"habeff",
  "tempeff",
  "trend",
  "WNSEff",
  "tempeff2",
  "montheff2",
  "treebateff",
  "windeff",   
  "pZero.int",
  "precipeff"
  #"watereff"
  #"hypersd"
  #"condeff"
)


###########################
##  NON-PARALLEL BUGS CALL, FOR DEBUGGING   [for now, this is the main way I'm running winBUGS- can't seem to get the parallel processing running]
###########################

nt=10; ni = 15000; nb = 5000; nc = 2    
model.file <- paste(filename_noAcoust,".bug",sep="")

Mod <- bugs(   inits=Inits, n.chains=nc, model.file=model.file, bugs.directory=bugs.directory, 
               data=Data, parameters.to.save=Pars, 
               n.thin=nt, n.iter=ni, n.burnin=nb, DIC=FALSE, debug=T,over.relax = FALSE)
	
m <- as.mcmc.bugs(Mod)         # convert BUGS output to CODA object (from "emdbook" package)


setwd(BUGSDir)
filename <- sprintf("WithAcoustic_%s.Rdata",Sys.Date())
save(list = ls(all.names = TRUE), file = filename, envir = .GlobalEnv)


#####################################################################

  #### plot out the abundance trajectory
mid <- apply(Mod$sims.list$abund.total[,2,,1],2,mean)
up <- apply(Mod$sims.list$abund.total[,2,,1],2,function(t) quantile(t,0.975))
low <- apply(Mod$sims.list$abund.total[,2,,1],2,function(t) quantile(t,0.025))
errbar(x=survYears,y=mid,yplus=up,yminus=low)
#errbar(x=mistnetYears[-1],y=mid,yplus=up,yminus=low)

###############################
#  VISUALIZE RESULTS
###############################

names(Mod$sims.list)

    # variation in abundance (log scale)
hist(Mod$sims.list$tempeff)       # shouldn't be negative, right!  Maybe need a quadratic term in there?
hist(Mod$sims.list$windeff)    # straddles zero. not converging
hist(Mod$sims.list$precipeff)   # positive!

hist(Mod$sims.list$p0.mistnet)
hist(Mod$sims.list$p0.acoust)  

hist(Mod$sims.list$pZero.int[,1])  # 30% prob of zero captures for mistnet surveys
hist(Mod$sims.list$pZero.int[,2])  # 50% prob of zero captures for acoustic surveys

hist(Mod$sims.list$treebateff)  # peak is at lower bound of 1...
  

###  function for plotting out species time series...

trendfunc <- function(time,slope,intercept){
  logint <- log(intercept)
  exp(logint + slope*(time-2012))
}

species= "Lcine" # "Mlucy" # "Efusc" # "Lbore" # "Mleib" #  "Msept" #  
acoust = T
plotSpecies <- function(species="Mlucy",BUGSoutput=Mod,trend=F,wns=F,svg=F){
  if(svg) setwd(imageDir)
  specnum <- which(allSpecies==species)
  expr <- sprintf("BUGSoutput$sims.list$abund.total[,%s,,%s]",specnum,1)
  dframe <- eval(parse(text=expr))
  mid <- apply(dframe,2,mean)
  up <- apply(dframe,2,function(t) quantile(t,0.975))
  low <- apply(dframe,2,function(t) quantile(t,0.025)) 
  if(trend){ 
    slope <- eval(parse(text=sprintf("mean(BUGSoutput$sims.list$trend[,%s])",specnum)))
    ndx <- ifelse(acoust,5,3)
    intercept <- eval(parse(text=sprintf("mean(BUGSoutput$sims.list$abund.total[,%s,ndx,2])",specnum)))
  }
  if(wns){
    ndx <- ifelse(acoust,6,4) 
    before <- eval(parse(text=sprintf("mean(BUGSoutput$sims.list$abund.total[,%s,2,3])",specnum)))
    after <- eval(parse(text=sprintf("mean(BUGSoutput$sims.list$abund.total[,%s,ndx,3])",specnum)))
  }

  filename <- paste(species,"_relAbund_withacoust.svg",sep="")
  mx1 = round(max(up)+0.05*max(up))
  mx2 = mx1/1000
  by1 = max(up)/10
  by2 = by1/1000
  if(svg) svg(filename=filename,width=4,height=4)
        if(!acoust) survYears = mistnetYears[-1]
	  errbar(x=survYears,y=mid,yplus=up,yminus=low,yaxt="n",main=species,
			xlab="Years", ylab="Abundance relative to Efusc 2012")
        abline(v=2008,lwd=3)
	  axis(2,at=seq(0,mx1,by=by1),labels=round(seq(0,mx2,by=by2),1))
        title(main=species)
        if(trend) curve(trendfunc(x,slope,intercept),from=2006, to=2013,lwd=1.5,lty=1,col=gray(0.3),add=T)
        if(wns){ lines(x=c(2006:2008),y=rep(before,times=3),lwd=2.5,lty=2,col=gray(0.5));lines(x=c(2008:2013),y=rep(after,times=6),lwd=2.5,lty=2,col=gray(0.5))}
  if(svg)dev.off()
}


### "Mlucy" # "Efusc" # "Lbore" # "Mleib" #  "Msept" # "Lcine" #

allSpecies
for(s in allSpecies){
  plotSpecies(s,Mod,trend=T,wns=T,svg=T)
}
 
s= "Lbore" # "Lcine" # "Mleib" # "Msept" # "Efusc" # "Mlucy" # "Lbore" # "Efusc" # "Lbore" # 
plotSpecies(s,Mod,trend=T,wns=T,svg=F)

 #str(m$sims.matrix)
 # ?as.mcmc

######################
 # head(summary(m),1000)  # review MCMC output: takes a while with so many parameters

varnames(m)

par(ask=T)
 #  plot(m) 



#######################
# Summarize catch per unit effort...


### MISTNET 
 
#### all bats...
t1 <- apply(counts.mistnet,2,function(t) sum(t,na.rm=T))  # total captures
t2 <- apply(neteffort,1,function(t) sum(t,na.rm=T))   # total effort
netCPUE <- t1/t2
netCPUE

 #### little brown...
t1 <- apply(counts.mistnet[1,,],1,function(t) sum(t,na.rm=T))  # total captures
t2 <- apply(neteffort,1,function(t) sum(t,na.rm=T))   # total effort
netCPUE <- t1/t2
netCPUE

cbind(survYears,netCPUE,micCPUE)


### ACOUSTIC 
 
#### all bats...
t1 <- apply(counts.acoust,2,function(t) sum(t,na.rm=T))  # total captures
t2 <- apply(micduration,1,function(t) sum(t,na.rm=T))   # total effort
micCPUE <- t1/t2
micCPUE

 #### little brown...
t1 <- apply(counts.acoust[1,,],1,function(t) sum(t,na.rm=T))  # total captures
t2 <- apply(micduration,1,function(t) sum(t,na.rm=T))   # total effort
micCPUE <- t1/t2
#t2
micCPUE
t1

sum(t1)
sum(t1)/sum(t2)

allSpecies

###################
####### VISUALIZE TREND ACROSS ALL SPECIES

trend <- data.frame(est=numeric(length(allSpecies)))
trend$group <- allSpecies
trend$upper <- 0
trend$lower <- 0
trend$group2 <- "NBAFS"

specNum <- c(1:length(allSpecies))
for(s in 1:length(allSpecies)){
  trend$est[s] <- median(as.vector(sapply(m[,sprintf("trend[%s]",specNum[s])], function(t) t)))
  trend$lower[s] <- quantile(as.vector(sapply(m[,sprintf("trend[%s]",specNum[s])], function(t) t)),0.05)
  trend$upper[s] <- quantile(as.vector(sapply(m[,sprintf("trend[%s]",specNum[s])], function(t) t)),0.95)
  #counter=counter+1
}
trend

graphics.off()
                   
setwd(imageDir)
rawname <- "trend_dotplot5_withacoust"
filename <- paste(rawname,"_",Sys.Date(),".svg")
svg(file=filename, width=3, height=5)        #Specify output graphics device:

#--Load extra library for plotting:
require(lattice)


#--Use "+" and filled square instead of "o" for plot symbols:
trellis.par.set(superpose.symbol=list(pch=c(3,15)))

#?reorder
dotplot.errors(x=trend,type.bar="CI",qlabel="Population growth trend")

#warnings()
#--Close plotting device:
dev.off()

###################################
####### VISUALIZE WNS EFFECT

wns <- data.frame(est=numeric(length(allSpecies)))
wns$group <- allSpecies
wns$upper <- 0
wns$lower <- 0
wns$group2 <- "NBAFS"

specNum <- c(1:length(allSpecies))
for(s in 1:length(allSpecies)){
  wns$est[s] <- median(as.vector(sapply(m[,sprintf("WNSEff[%s]",specNum[s])], function(t) t)))
  wns$lower[s] <- quantile(as.vector(sapply(m[,sprintf("WNSEff[%s]",specNum[s])], function(t) t)),0.05)
  wns$upper[s] <- quantile(as.vector(sapply(m[,sprintf("WNSEff[%s]",specNum[s])], function(t) t)),0.95)
  #counter=counter+1
}
wns

graphics.off()
                   
setwd(imageDir)
rawname <- "wns_dotplot5"
filename <- paste(rawname,"_",Sys.Date(),".svg")
svg(file=filename, width=3, height=5)        #Specify output graphics device:

#--Load extra library for plotting:
require(lattice)


#--Use "+" and filled square instead of "o" for plot symbols:
trellis.par.set(superpose.symbol=list(pch=c(3,15)))

#?reorder
dotplot.errors(x=wns,type.bar="CI",qlabel="WNS effect")

#warnings()
#--Close plotting device:
dev.off()



###############################
#  WRITE TABLE OF ALL MODEL PARAMETERS (for base model- fully time-dependent)
###############################

GetQuantiles <- function(vec,r){
  q <- numeric(3)
  q[2] <- round(quantile(vec,0.025),r)
  q[1] <- round(mean(vec),r)
  q[3] <- round(quantile(vec,0.975),r)
  return(q)
}


ParamNames <- c(
	"Abundance*",
	"P0, Mistnet",
	"P0, Acoustic",
  "Tree-bat Effect (logit)",
	"Temp Effect (logit)",
	"Precip Effect (logit)",
	"Wind Effect (logit)",
	"Habitat Effect (logit)",
	"Dist Water Effect (logit)",
	"Prob of zero captures",
  "Prob of zero acoust detections"
      #"Process variance (random int)"
)


ParamsTable <- data.frame(
    Params = ParamNames,
    PriorMean = 0,
    PriorLow95 = 0,
    PriorHigh95 = 0,
    PosteriorMean = 0,
    PosteriorLow95 = 0,
    PosteriorHigh95 = 0
)

### NAs for those params that are not yet included in the model
cols <- c(2:7)
rows <- c(8,9)
prcols <- c(2:4)
ParamsTable[rows,cols] <- NA
     
############
### base capture probability, mist net: p0.mistnet

hist(Mod$sims.list$p0.mistnet)
prior <- GetQuantiles(rbeta(1000000,1,1),3)
posterior <- GetQuantiles(Mod$sims.list$p0.mistnet,4)

ParamsTable[2,cols] <- c(prior,posterior)

############
### base capture probability, acoustic

hist(Mod$sims.list$p0.acoust)
prior <- GetQuantiles(rbeta(1000000,1,1),3)
posterior <- GetQuantiles(Mod$sims.list$p0.acoust,4)

ParamsTable[3,cols] <- c(prior,posterior) 

############
### tree bat effect

hist(Mod$sims.list$treebateff)
prior <- GetQuantiles(runif(1000000,0,10),2)
posterior <- GetQuantiles(Mod$sims.list$treebateff,2)

ParamsTable[4,cols] <- c(prior,posterior)


############
### temperature effect

hist(Mod$sims.list$tempeff)
prior <- GetQuantiles(rnorm(1000000,0,sqrt(10)),2)
posterior <- GetQuantiles(Mod$sims.list$tempeff,2)

ParamsTable[5,cols] <- c(prior,posterior)

############
### precip effect

hist(Mod$sims.list$precipeff)
prior <- GetQuantiles(rnorm(1000000,0,sqrt(10)),2)
posterior <- GetQuantiles(Mod$sims.list$precipeff,2)

ParamsTable[6,cols] <- c(prior,posterior)

############
### wind effect

hist(Mod$sims.list$windeff)
prior <- GetQuantiles(rnorm(1000000,0,sqrt(10)),2)
posterior <- GetQuantiles(Mod$sims.list$windeff,2)

ParamsTable[7,cols] <- c(prior,posterior)

############
### probability of zero captures, mist net

hist(Mod$sims.list$pZero[,1])
prior <- GetQuantiles(rbeta(1000000,1,1),3)
posterior <- GetQuantiles(Mod$sims.list$pZero[,1],4)

ParamsTable[10,cols] <- c(prior,posterior)

############
### prob zero captures, acoustic

hist(Mod$sims.list$pZero[,2])
prior <- GetQuantiles(rbeta(1000000,1,1),3)
posterior <- GetQuantiles(Mod$sims.list$pZero[,2],4)

ParamsTable[11,cols] <- c(prior,posterior)


# abundance

ParamsTable[1,prcols] <- GetQuantiles(runif(1000000,100,20000),0)

################
# SAVE TABLE
################

setwd(resultsDir)
filename <- "ParamsTable_baseBUGSmodel2.csv"
write.csv(ParamsTable,row.names=F,file=filename)


#############################
#   LOAD DOTPLOT ERROR BAR FUNCTION: from dotplots.errors, a new R function to ease the pain of creating dotplots 
#     written by Marcin Kozak
##############################
dotplot.errors <- function(x, myTheme = simpleTheme(pch = 19, col = 1), 
                           qlabel = "Estimate", add.text.to.qlabel = "", type.bar = "SE", 
                           conf.level = .95, end.length = .05, reorder.groups = TRUE, 
                           reordering = "decreasing", label.define = list(), reference.line = NULL, 
                           bar.color = 1, horizontal = TRUE, ...) 
{ require(lattice) 
  o <- c(which(colnames(x) == "group"), which(colnames(x) == "lower"), 
         which(colnames(x) == "est"), which(colnames(x) == "upper")) 
  if (length(o) != 4) stop("Error: Incorrect data frame") 
  x <- x[, o] 
  x$group <- factor(x$group) 
  if (horizontal == T) hor <- 1 else hor <- -1 
  if (reordering == "decreasing") 
    FUN.to.reorder <- function(x) mean(x) * hor 
  else FUN.to.reorder <- function(x) -mean(x) * hor 
  if (reorder.groups) x$group <- 
    with(x, reorder(group, est, FUN = FUN.to.reorder)) 
  if (type.bar == "CI") xlab <- substitute(expression(lab %+-% CL), 
                                           list(lab = qlabel, CL = paste(" ", as.character(conf.level * 100), 
                                                                         "% CI", add.text.to.qlabel, sep = ""))) else 
                                                                           if (type.bar == "SE") 
                                                                             xlab <- substitute(expression(lab %+-% " SE" ~ ~ AT), 
                                                                                                list(lab = qlabel, AT = add.text.to.qlabel)) else 
                                                                                                  xlab <- substitute(expression(lab %+-% Type.bar ~ ~ AT), 
                                                                                                                     list(lab = qlabel, Type.bar = type.bar, AT = add.text.to.qlabel)) 
  
  if (horizontal == T) 
    p <- stripplot(group ~ est, data = x, 
                   lower = x$lower, upper = x$upper, 
                   xlim = range(x[,2:4])+c(-.04,.04)*diff(range(x[,2:4])), 
                   xlab = c(list(xlab), label.define), 
                   par.settings = myTheme, 
                   panel = function(x, y, lower, upper, ..., subscripts) { 
                     if (is.null(reference.line) == F) 
                       panel.abline(v = reference.line, col = "grey") 
                     panel.abline(h = y, col = "grey", lty = "dashed") 
                     panel.arrows(x0 = lower[subscripts], y0 = y, 
                                  x1 = upper[subscripts], y1 = y, angle = 90, code = 3, 
                                  length = end.length, col = bar.color) 
                     panel.stripplot(x, y, ...) }, ...) else 
                       p <- stripplot(est ~ group, data = x, lower = x$lower, upper = x$upper, 
                                      ylim = range(x[,2:4])+c(-.04,.04)*diff(range(x[,2:4])), 
                                      ylab = c(list(xlab), label.define), 
                                      par.settings = myTheme, 
                                      panel = function(x, y, lower, upper, ..., subscripts) { 
                                        if (is.null(reference.line) == F) 
                                          panel.abline(h = reference.line, col = "grey") 
                                        panel.abline(v = x, col = "grey", lty = "dashed") 
                                        panel.arrows(y0 = lower[subscripts], x0 = x, 
                                                     y1 = upper[subscripts], x1 = x, angle = 90, code = 3, 
                                                     length = end.length, col = bar.color) 
                                        panel.stripplot(x, y, ...) }, ...) 
  
  print(p) 
  invisible(p) 
} 








