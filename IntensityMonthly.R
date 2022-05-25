## Sensitivity Analysis of cluster param choices for travel bands
## Author: Rachel Cardell-Oliver
## Version: July 2021 Revised May 2022

#libraries and constants
{
  source("./GlobalConstants.R") 
}


#given a Jny matrix of M passengers and N times and 
# lower (>= B1) and upper (< B2) bounds for number of journeys
# vector of tokens per psg and tkmain vector of token names (strings) to track (eg top 20)
# return a list with
# B1 and B2 
# size of the band (number of months)
# token breakdown percentage of riders for each band (except Abs) 
# history of number of journeys per month and number of passengers per month
# this is simplified - see more detailed fn see below
{

#returns B1,B2,size,NJY,token counts,kld (10 items) for journey count in range >=B1 <B2
getTravelBandMetrics <- function(B1, B2, 
                                 journeyMat=journeyMatrix, 
                                 tokenGroupMat=tokenGroupMatrixMonthly, #needs monthly tokens
                                 tokexpected=tokenProbabilities) #populatio probabilities
{
  inband = c((journeyMat >= B1 & journeyMat < B2)*1)
  jnyband = c(journeyMat) * inband
  tokinband = c(tokenGroupMat)[which(inband==1)]
  size = sum(inband==1) #total number of months in the band - size of the cluster group
  M = dim(journeyMat)[2] #number of months in the dataset
  NJY = sum(jnyband) #total number of journeys taken in this band
  
  if (size>200) 
  {
    tokct = table(tokinband)
    tokct = tokct[order(names(tokct))] #ensure still alpha order
    tokprob = (tokct/sum(tokct))
    kld = KLdiv(tokprob, tokexpected)
    ## check name order OK for small groups
    if (sum(names(tokprob)==names(tokexpected))!=5) {
      print(paste(c("WARNING unexpected token name order",names(tokprob),
                    " != ",names(tokexpected)),collapse=" "))
    }
  
    #return metrics used for band selection
    bandMetrics = c(B1,B2,size,NJY,tokct,kld) 
  } else { #insufficient size
    bandMetrics = c(B1,B2,size,NJY, 0,0,0,0,0, 0)
  }
  return (bandMetrics)
}
}

## calculate travel band metrics (slow so only once)
## outputs "./results/bandSelectionData.csv"
{

#top travel band - cut high intensity monthly counts after 60
quantile(c(journeyMatrix[which(journeyMatrix>0)]),c(90,95,98,seq(99,100,by=0.05))/100) 
#90%    95%    98%    99% 99.05%  99.1% 99.15%  99.2% 99.25%  99.3% 99.35%  99.4% 99.45%  99.5% 99.55%  99.6% 99.65%  99.7% 99.75%  99.8% 99.85% 99.9% 99.95%   100% 
#32     37     42     47     47     47     48     48     49     49     50     51     51     52     53     54     55    56     58      59    62    65     71     209 
# 99.8-th percentile falling at 60

bandSelectionData = matrix(nrow=60,ncol=10) 
for (c in 1:60) { 
  bandSelectionData[c,] = getTravelBandMetrics(c,c+1)
}
colnames(bandSelectionData) = c("B1","B2","Freq","NJY",names(tokenProbabilities),"KLD")

#save data in results
write.csv(bandSelectionData,file="./results/bandSelectionData.csv",row.names=FALSE)

}



## FIG TO DISPLAY BANDS RESULTS
{
  
  pdf(file = "./results/travelBandsAnalysis.pdf",
      width = 12,height = 4)
  rr=1:60 #86 is 99.99 quantile but 60 enough to show
  {
    #journey counts
    par(mfrow=c(1,3))
    plot(bandSelectionData[rr,3]/1000000,
         type="h",lwd=2,main="Frequency All Passengers",
         ylab="Total Months (millions)", # (log scale)",
         xlab="Journeys per Month")
    abline(v=freqbounds-0.5,col="blue") #bounds from globalConstants
    
    ## per-Token frequencies
    toks = bandSelectionData[rr,5:9]
    plot(log(toks[,4],base=10),
         type="l",lwd=2,col=tokencols[order(tokennames)][4],
         main="Frequency per Passenger Type",
         ylab="Total Months (log scale)",
         xlab="Journeys per Month",
         ylim=range(log(toks,base=10)))
    for (c in 1:5) {
      lines(log(toks[,c],base=10),lwd=2,col=tokencols[order(tokennames)][c])
    }
    abline(v=freqbounds-0.5,col="blue") 
    legend("topright",legend=tokennames,col=tokencols,lty=c(1,1,1,1,1),lwd=c(2,2,2,2,2))
    
    # KLD
    plot(bandSelectionData[rr,10],
         type="l",lwd=2,main="Mix of Passenger Types",
         ylab="KLD",xlab="Journeys per Month")
    abline(v=freqbounds-0.5,col="blue") #lty="dotted"
    par(mfrow=c(1,1))
  }
  dev.off()
  
  
}


## TABLE of band significance for the chosen bands for paper
{

  ccallmetrics = getTravelBandMetrics(freqbounds[1],freqbounds[length(freqbounds)])

  bandtable = c()
  for (i in 1:5) {
    ccmetrics = getTravelBandMetrics(freqbounds[i],freqbounds[i+1])
    bandtable = rbind(bandtable, ccmetrics)
  }
  bandtable = rbind(bandtable,ccallmetrics)
  
  to = c(4,2,3,5,1)  #token order Std,Sch,Snr,Ter,Con
  bttok = bandtable[,5:9][,to]
  colnames(bttok) = tokenfullnames

  perbandtotal = rowSums(bttok)
  perbandpct = round(bttok*100 / perbandtotal)
  
  bandtable = cbind(bandtable[,1:4],perbandpct ,bandtable[,10])
  bandtable[,2]= bandtable[,2]-1 #report inclusive upper and lower band eg [3,16] not [3,17)s
  bandtable[5:6,2] = 60 #60+ is name for upper bound
  bandtable[,3:4] = round(bandtable[,3:4]/1000000,digits=1)    #cluster size million months #journey size million journeys
  bandtable[,10] = round(bandtable[,10],digits=3) #kldiv
  colnames(bandtable)[1:4] = c("B1","B2","Band size (months)","Journeys")
  colnames(bandtable)[10] = c("KLdiv")
  rownames(bandtable) = c(paste("I",1:5,sep=""),"Total")
  
  write.csv(bandtable,file="./results/intensityBandSummary.csv",row.names=TRUE)
  
  xtable(
    cbind(
      c(paste("I",1:5,sep=""),""),
      c("Rare (R)","Low (L)","Medium (M)","High (H)",
        "Very High (V)","Population"),
      bandtable),
    digits=c(0, 0,0, 0,0, 1,1, 0,0,0,0,0, 2),
    caption="Bands for monthly intensity of travel",
    label="tab:freqbands")

}


## prepare passenger intensity data for AnnualEngagment clustering
# converts journey counts into passenger intensity bands
# outputs ./data/passengerClassesMonthlyByYear.rds matrix to be used for AnnualEngagement clustering
{
NM =  dim(journeyMatrix)[2] #number of months
passengerClassesMonthly = journeyMatrix * 0 #initialise 
for (c in 1:NM) 
{
  passengerClassesMonthly[,c] = getTravelBand(journeyMatrix[,c])
}

saveRDS(passengerClassesMonthly,"./data/passengerClassesMonthly.rds") 

{
  passengerClassesMonthlyByYear = 
    rbind(passengerClassesMonthly[,1:12],
          passengerClassesMonthly[,12+1:12],
          passengerClassesMonthly[,24+1:12],
          passengerClassesMonthly[,36+1:12],
          passengerClassesMonthly[,48+1:12])
}

saveRDS(passengerClassesMonthlyByYear,"./data/passengerClassesMonthlyByYear.rds") 
# in AE we will convert these to histogram of class distributions per year
}

