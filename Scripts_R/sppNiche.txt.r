#Creator: R.Loiselle


# This script evolved from NicheMod2.r for the analysis of P/A data from the NOAA RACE database.
# The script was built for the purpose of simplifying access and extraction of desired species throught the use of Matrices and their Indices

#Please skip ahead to data analysis section to view syntax related to building and accessing matrix data...   And additional Notes....

#===========================================================================================================================================================

## Install Packages necessary for script functioning
library(RODBC)
library(reshape)
library(Hmisc)
library(cluster)
library(gtools)
library(MASS)
library(vegan)
library(labdsv)
library(mgcv)
library(splines)
library(quantreg)
library (fBasics)   #matrices
library(maps)
library(fUtilities) #color ramps
library(RColorBrewer)

#======================================================================================================================================================
## Data Prep

#===========================================================================================================================================================

##**** AK DATA ****

## Set the working directory for the files to  be processed
arctic <- setwd("D:\\Projects\\Arctic Fishes\\Archive\\Niche_Modeling\\AK")
## Save the working directory into a variable called arctic.dir
arctic.dir <- getwd()
## Display the value stored in the arctic.dir variable to validate it is the right directory
arctic.dir
## Create a character vector listing the files contained in the working directory
files.arctic <-dir(pattern='.csv')
## List the contents of the character vector
files.arctic
## Define the columns in the files--> identifies the mode of the data found in columns, so that R understands how the data can be used
col.defs <- c(rep("numeric", 2),"character",rep("numeric", 2),rep("character", 2),rep("numeric",2),rep("character",2),rep("numeric", 7))
    col.defs
# Loop through all the files in the directory reading them into objects named with the file name
##(assigning the column definitions to each column in each table for all files (elements) in vector files.arctic
for (i in files.arctic) {
    assign(i, read.csv(i, colClasses = col.defs))         
    }

# Combining all datasets into one
## set up a list of file names to save into
ak <- NULL   #aknames <- files.arctic

## Combining all datasets into one
for(i in files.arctic) {
    tmp <- get(i)
    ak <- rbind(ak,tmp)
}

# check work
# for(i in files.arctic) {
    tmp3 <- get(i)
    print(dim(tmp3))
# }
    #36066+46648+28230+51086+49294+54726+62278+39921+23111+59186+49377+17871

## get list of column names
aknm <- names(ak)

#Loop through column elements, assigning NA to the digits
for(i in aknm) {
    ak[[i]][ak[[i]] %in% -9999] <- NA
    }

#making sure all longitudinal values are neg. variables (proper spatial locations)
    ak$LONGITUDE <- (-1) * abs(ak$LONGITUDE)   

## **** WC DATA ****
## Data Prep

## Set the working directory for the files to  be processed
setwd("D:\\Projects\\Arctic Fishes\\Archive\\Niche_Modeling\\WC")
## Save the working directory into a variable called arctic.dir
wc.dir <- getwd()
## Display the value stored in the arctic.dir variable to validate it is the right directory
wc.dir
## Create a character vector listing the files contained in the working directory
files.wc <-dir(pattern='.csv')
## List the contents of the character vector
files.wc
## Define the columns in the files--> identifies the mode of the data found in columns, so that R understands how the data can be used
col.defs <- c(rep("numeric", 2),"character","numeric",rep("character", 2),rep("numeric",4),rep("character",2),rep("numeric", 8))
    col.defs
# Loop through all the files in the directory reading them into objects named with the file name
##(assigning the column definitions to each column in each table for all files (elements) in vector files.arctic
for (i in files.wc) {
    assign(i,read.csv(i, colClasses = col.defs))         
    }

## Combining all datasets into one
## set up a list of file names to
wc <- NULL
## Combining all datasets into one
for(i in files.wc) {
    tmp3 <- get(i)
        #tmp4 <- dim(tmp3)
    wc <- rbind(wc,tmp3)
}

#check work
#for(i in files.wc) {
#    tmp3 <- get(i)
#    print(dim(tmp3))
#}
    #9346+9855+2394+1678+1288+3174+2614+2340+3951+4576+9774+7972+9304+14557+11306+8725+5952+9663+10128+10211+10503+11669+14277

## Assign -9999 as missing data records for all columns
wcnm <- names(wc)
for(i in wcnm) {
    wc[[i]][wc[[i]] %in% -9999] <- NA
}

## Correct coordinates
wc$LONGITUDE <- (-1) * abs(wc$LONGITUDE)
#Formatting date field
wc$DATE <- as.Date(wc$DATE, format='%m/%d/%Y')
## extracting year field from data
wc$YEAR <- substring(wc$DATE,1,4)

## Combine Cruise-Haul fields for unique identifier
wc$crha <- paste(wc$CRUISE,wc$HAUL,sep='.')
## Give station name using Cruise-Haul Identifier where Station ID is missing
for(i in 1:length(wc$STATION))  {
    if(wc$STATION[i] %in% NA) {
        wc$STATION[i] <- wc$crha[i]
    }
}

## remove columns which won't match AK data
wc2 <- subset(wc, select=c(-WEIGHT, -NUM_FISH, -HAULJOIN, -crha))
## reorder columns to match AK data
wc2 <- wc[,aknm]

## Merge both datasets and Prep for anlaysis
data2 <- rbind(ak, wc2)

## Removing names from the Scientific column which are too general
    #uniqspp <- sort(unique(data2$SCIENTIFIC))
GenusSppTMP = grep("^[a-zA-Z]*$|^.* (sp\\.).*$|^.* unid.*$|^.*sponge|^.*[Ee]gg|^.* hybrid|^.* species.*$|^[0-9].*$", data2$SCIENTIFIC, value=TRUE)
    data2 <- data2[!data2$SCIENTIFIC %in% GenusSppTMP,]
## Testing data for correctness
uniqspp <- sort(unique(data2$SCIENTIFIC))
    #testuniq <- unique(data2[,11])

#===========================================================================================================================================================

##Analysis

#===========================================================================================================================================================


  #NOTES:

     #   1.) These Matrices are designed to allow access to Species Of Interest (6) at a given location for a given year through the use of Placement Indices.
                  #  a.) Ie. not generalized!!
                  #  b.) Acesses only species of interest

     #     2.) Also used are hard coded color matrices for graphing purposes, but only interested in presences, not if a station was sampled in that year
     #     3.) Produces .pdf files of yearly PA locations and temp/depth profiles for each spp.


#===========================================================================================================================================================

#Testing data structure
#reformatting the data...getting unique data for each station/year combination
station <- subset(data2, select=c(LATITUDE, LONGITUDE, STATION, STRATUM, YEAR, BOT_DEPTH, BOT_TEMP, SURF_TEMP))
#melt to long format
meltstation <- melt(station, id=c("STATION", "YEAR", "STRATUM"), measured=c("LATITUDE", "LONGITUDE", "WTCPUE", "NUMCPUE", "BOT_DEPTH", "BOT_TEMP", "SURF_TEMP"), na.rm=TRUE)

#cast to wide format, averaging variables
station <- cast(meltstation, STATION + YEAR + STRATUM~ variable, mean)

#show counts of records per month
table(months(station$DATE))
#extract month from date and assign data class
habitat$month <- as.numeric(format(habitat$DATE, '%m'))

## These only had one year of data, check in original data
#data2 <- data2[!data2$STATION %in% c("S-32","T-31", "U-30", "U-24", "V-29"),]#these only had one year of data, check in original data
## These data are clustered for other collection purposes
#tmp <- grep("^[a-zA-Z][a-zA-Z][0-9]{4}$", data2$STATION, value=TRUE)
# Testing correctness
# unique(tmp)
## Subsetting for use in analysis
#data2 <- data2[!data2$STATION %in% tmp,]

## Get rid of unnecessary columns   ** FINAL DATAFRAME **
data2 <- subset(data2, select=c(-VESSEL, -CRUISE, -HAUL))

#===========================================================================================================================================================

## Remove ls() objects, except data2, for RAM
  rm(list=(ls(pattern="^[a-zA-Z]{1,3}[0-9]*$|^[a-zA-Z]{5,15}[0-9]*$|^.*csv$")))

#===========================================================================================================================================================

## Prepare station Presence/Absence

## Create an indices of Species/Year/Station numbers for access in Presence/Absence (PA) matrices
spp <- sort(unique(data2$SCIENTIFIC))
#spp2 <- make.names(spp, unique=TRUE)
sta <- sort(unique(data2$STATION))
year <- sort(unique(data2$YEAR))

## Station PA from sampling year
StaBD <-  xtabs(~ STATION + YEAR, data=data2)  # this one also conveys species richness
stapa <- ifelse(StaBD[1:11389,] > 0, 1, 0)
#change class
stapa <- as.array(stapa)

## Make syntactically valid names for calling
data2$SCIENTIFIC <- make.names(data2$SCIENTIFIC)

## list species of interest and gather PA data using this list
spplist <- make.names(c("Clupea pallasi", "Hippoglossoides elassodon", "Limanda aspera", "Microgadus proximus", "Pleuronectes quadrituberculatus", "Theragra chalcogramma"))

## Get PA for spp
#create empty object
spp.data <- NULL
#loop through species list,extract data for each species and bind by row for long format (melt-cast wasn't working due to size)
for(i in spplist) {
    tmp <-NULL
    tmp <-  data2[data2$SCIENTIFIC %in% i,]
    #print(dim(tmp))
    spp.data <- rbind(spp.data, tmp)
}

#Spp PA per year per station
PA <-  xtabs(~ STATION + YEAR + SCIENTIFIC, data=spp.data) 
#Number species present at a station per year
#staPA <-  xtabs(~ YEAR + STATION, data=data)

#if counts are greater than 0, assign 1 (present)
PA <- ifelse(PA[,,] > 0, 1, 0)
          
## Checking work... Duplicate sampling up to 5 times....
  #which(PA==2)
  #x <- c(7550, 26,6)
    #tmp <- array(1:1177800, dim=x)
  #for(i in spplist) {
  #for( st in sta) {
  #tmp <- data2[data2$SCIENTIFIC %in% i,]
  #tmp2 <- tmp[tmp$STATION %in% st,]
  #duplicated(tmp2[tmp2$YEAR,])
  #}
  #}
  
unique(spp.data$STATION)      # should = 7550
any(PA[,,]>1)     #should be false

#===========================================================================================================================================================

#Only interested in presences

#creating dimensions
x <- c(7550, 6)
#create array with the dimensions
PAf <- array(dim=x)

#loop to hard code colors into a matrix where species are absent
#looking at only the first six species
for(sp in 1:6)  {
    #All 26 years in data
    for(y in 1:26) {
        #all sampling stations
        for(st in 1:7550) {
            #if absent, assign "blue"
            if (PA[st,y,sp] > 0) {
            PAf[st,sp] ="blue"
            }
        }
    }
}

#===========================================================================================================================================================

## Reformatting the data...getting unique data for each station/year combination
station <- subset(data2, select=c(LATITUDE, LONGITUDE, STATION, STRATUM, YEAR))
meltstation <- melt(station, id=c("STATION", "YEAR", "STRATUM"), measured=c("LATITUDE", "LONGITUDE"), na.rm=TRUE)
stacoor <- cast(meltstation, STATION ~ variable, mean)

##Get variables into array form as with PA data
## Environmental data..... IF station was sampled (ie.PA=1), then
enviro <- subset(data2, select=c(STATION, YEAR, BOT_DEPTH, BOT_TEMP, SURF_TEMP))
    meltenviro <- melt(enviro, id=c("STATION", "YEAR"), measured=c("BOT_DEPTH", "BOT_TEMP", "SURF_TEMP"), na.rm=TRUE)

#dataframe form
enviro <- cast(meltenviro, YEAR + STATION ~ variable, mean)

#Array form
env2 <- cast(meltenviro, STATION ~ YEAR ~variable, mean)

##Get Variables for spp only data
spp.data2 <- subset(spp.data, select=c(STATION, YEAR, BOT_DEPTH, BOT_TEMP, SURF_TEMP))
meltspp.data <- melt(spp.data2, id=c("STATION", "YEAR"), measured=c("BOT_DEPTH", "BOT_TEMP", "SURF_TEMP"), na.rm=TRUE)
#dataframe form
#spp.enviro <- cast(meltspp.data, YEAR + STATION ~ variable, mean)
#Array form
spp.env2 <- cast(meltspp.data, STATION ~ YEAR ~variable, mean)


#===========================================================================================================================================================

#Plotting

#===========================================================================================================================================================

#setting up the plot frame
pdf(file="sppniche",8.5,11, paper="USr")
  layout(matrix(c(1,1,1,1,1,1,2,4,6,3,5,7), ncol=4))#,
       #widths=matrix(c(rep(lcm(6),6), rep(lcm(5),6), ncol=4)))#,
        # heights=rep(lcm(5),12)),
         #respect=TRUE)

#layout.show(7)

par(mai=c(1,1,1,1),mar=c(4,4,2,2))

#Start of the loop that goes through all the species in the dataset

x <- c(6184,26,3)
#sp=106
for(sp in 1:6) {
    tmp <- array(dim=x)
    for(y in 1:26) {
        for(v in 1:3){
            for(st in 1:6184) {
                if(PA[st,y,sp] > 0)  {
                    tmp[st,y,v] = env2[st,y,v]
                }
            }
        }
    }
    #map plot
    plot(stacoor$LATITUDE~stacoor$LONGITUDE, main=c(spplist[sp]), pch=16, cex=.3, col=as.character(PAf[,sp]), ylim=c(30, 65), xlim=c(-180, -115), xlab="Longitude", ylab="Latitude", cex.main=1.6)
    map('world', xlim = -c(179.2, 119.0), ylim = c(33.0, 65.0), add=TRUE) # To get entire boundary to plot with single set of data points, might need to add plot=FALSE, then send to vector and call plot after legend (eg. richness)
    legend("bottomleft", pch=c(16,16), col="blue", legend= paste("Spp Present, N=", sum(PAf[,sp] %in% "blue"), sep=""))
    
    ## Sample Depth
    length1 <- length(env2[,,1][!is.nan(env2[,,1])])
    if(length1 > 0) {
        h1 <- hist(env2[,,1], col="lightgray", xlab="Sample Depth, m", ylab="Number of samples", main = "", xlim=c(0,1600))
        # text(1000, (max(h1$counts)+(0.15*max(h1$counts))), "Sample depth", cex=1.5, pos=4)
        text(1000, (max(h1$counts)-(0.10*max(h1$counts))), "All Samples")
        text(1000, (max(h1$counts)-(0.20*max(h1$counts))), paste("N=", length1, sep=""), cex=.90)
        #text(1000, 0-(0.4*max(h1$counts)), "sample depth, m")
    } else { plot(1, type="n", axes=T, xlab="", ylab="")
        text(1,1, labels=paste("N=", length1, sep=""))
    }
    ## Spp Depth
    length1 <- length(tmp[,,1][!is.na(tmp[,,1])])
  if(length1 > 0) {
    h1 <- hist(tmp[,,1], col="lightgray", xlab="Species Depth, m", ylab="", main="", xlim=c(0,1600))
      #text(1000, (max(h1$counts)-(0.10*max(h1$counts))))
      text(1000, (max(h1$counts)-(0.20*max(h1$counts))), paste("N=", length1, sep=""), cex=.90)
  } else { plot(1, type="n", axes=T, xlab="", ylab="")
            text(1,1, labels=paste("N=", length1, sep=""))
          }
    ## Bottom Temp
    length1 <- length(env2[,,2][!is.nan(env2[,,2])])
    if(length1 > 0) {
        h1 <- hist(env2[,,2], col="lightgray", xlab="Bottom Temperature, C", ylab="Number of samples", main="", xlim=c(-5,15))
        #text(0, (max(h1$counts)+(0.15*max(h1$counts))), "Bottom Temperature", cex=1.5, pos=4)
        text(12, (max(h1$counts)-(0.10*max(h1$counts))), "All\nSamples", cex=.90)
        text(12, (max(h1$counts)-(0.35*max(h1$counts))), paste("N=", length1, sep="\n"), cex=.90)
        #text(10, 0-(0.4*max(h1$counts)), "temperature, C")
    } else { plot(1, type="n", axes=T, xlab="", ylab="")
            text(1,1, labels=paste("N=", length1, sep=""))
    }
    ## Spp Bottom Temp
    length1 <- length(tmp[,,2][!is.na(tmp[,,2])])
    if(length1 > 0) {
        h1 <- hist(tmp[,,2], col="lightgray", xlab="Species\nBottom Temperature, C", ylab="", main="", breaks=10, xlim=c(-5,13))
        #text(-1, (max(h1$counts)-(0.10*max(h1$counts))))
        text(-1, (max(h1$counts)-(0.20*max(h1$counts))), paste("N=", length1, sep=""), cex=.90)
    } else { plot(1, type="n", axes=T, xlab="", ylab="")
            text(1,1, labels=paste("N=", length1, sep=""))
    }
    ## Surface Temp
    length1 <- length(env2[,,3][!is.nan(env2[,,3])])
    if(length1 > 0) {
        h1 <- hist(env2[,y,3], col="lightgray", xlab="Surface Temperature, C", ylab="Number of samples", main="", xlim=c(-5,20))
        #text(0, (max(h1$counts)+(0.15*max(h1$counts))), "Surface Temperature", cex=1.5, pos=4)
        text(0, (max(h1$counts)-(0.10*max(h1$counts))), "All\nSamples", cex=.90)
        text(0, (max(h1$counts)-(0.35*max(h1$counts))), paste("N=", length1, sep="\n"), cex=.90)
        #text(125, 0-(0.4*max(h1$counts)), "temperature, C")
    } else { plot(1, type="n", axes=T, xlab="", ylab="")
            text(1,1, labels=paste("N=", length1, sep=""))
    }
    ## Spp Surface Temp
    length1 <- length(tmp[,,3][!is.na(tmp[,,3])])
    if(length1 > 0) {
        h1 <- hist(tmp[,,3], col="lightgray", xlab="Species\nSurface Temperature, C", ylab="", main="", breaks=10, xlim=c(-5,20))
        #text(1, (max(h1$counts)-(0.10*max(h1$counts))))
        text(1, (max(h1$counts)-(0.20*max(h1$counts))), paste("N=", length1, sep=""), cex=.90)
    } else { plot(1, type="n", axes=T, xlab="", ylab="")
            text(1,1, labels=paste("N=", length1, sep=""))
    }
    #savePlot(filename=paste(spplist[sp], sep=' '), type="pdf")
}

#turn of plotting device
dev.off()