## Investigate interesting subgroups of passengers based on 
## CIAM churn, intensity, annual engage and evolution categories
## also generate example passenger figure for motivation
## Author: Rachel Cardell-Oliver
## Version: 2020 Revised July 2021, May 2022


{
  source("./GlobalConstants.R") #load datasets, define utility fns and journey band definitions
  passengerClassesMonthly = readRDS("./data/passengerClassesMonthly.rds")
  engagePerYear = readRDS("./data/engagePerYear.rds")
}


## count number of months active and min and max engagement per passenger
{
  
  firstmonth = function(xx) { nxx = which(xx>0); return(min(nxx)) }
  whenstart = apply(journeyMatrix,1,firstmonth)
  lastmonth = function(xx) { nxx = which(xx>0); return(max(nxx)) }
  whenend = apply(journeyMatrix,1,lastmonth)
  numactivemonths = whenend - whenstart + 1

  maxintensity = apply(passengerClassesMonthly,1,max) 

  numactiveyears = rowSums(engagePerYear > 0)
  minAE = apply(engagePerYear,1,min)
  maxAE = apply(engagePerYear,1,max)
  AEinfirstyear = function(ts) { return( ts[which(ts!=0)[1]] ) }
  firstAE = apply(engagePerYear,1,AEinfirstyear)
}


## define interesting subgroups of passengers
{
  
  #SHORT BURST active <=12 months and intensity>=I3 at least once
  isshortburst = (numactivemonths<12 & maxintensity>=3)
  ShortburstPsg = which(isshortburst)
  
  #OCCASIONAL anything low = non-commuters, never have high engagement
  isoccasional = (numactiveyears>=1 & maxAE<=3)
  OccasionalPsg = which(isoccasional)
  
  # DRIFTING COMMUTER steps down from high start
  isdrifting = numactiveyears>=3 & firstAE>=4 &  minAE<4 
  CommuterDrift = which(isdrifting)
  
  # LOYAL COMMUTER stays high
  isloyal = numactiveyears>=4 & firstAE==5 & minAE>=4
  CommuterLoyal = which(isloyal)
  
  # OTHER all the rest
  NP = dim(journeyMatrix)[1]
  whichcaptured = sort(c(ShortburstPsg,OccasionalPsg,CommuterDrift,CommuterLoyal))
  Other = setdiff(1:NP,whichcaptured)
  
}

## table of properties of interesting subgroups
{
  #given a list of passenger IDs, find characteristics of that group
  psggroupsummary <- function(whichpassengers) 
  {
    numpsg = length(whichpassengers)
    numjny = sum(colSums(journeyMatrix[whichpassengers,]))
    
    tokct= table(tokenGroupMatrixMonthly[whichpassengers,])
  
    groupProb = tokct[2:6]/sum(tokct[2:6]) #all except Abs
    if (sum(names(groupProb)==names(tokenProbabilities))!=5) {
      print(paste(c("WARNING unexpected token name order nc=",nc,"cl=",ae,names(groupProb),
                    " != ",names(tokenProbabilities)),collapse=" "))
      print(groupProb)
    }
    kld = round(KLdiv(groupProb,tokenProbabilities),digits=3)
    
    #reordering of token names from alpha "Con" "Sch" "Snr" "Std" "Ter" 
    to = c(4,2,3,5,1)  #token order Std,Sch,Snr,Ter,Con
   
    summary = c(numpsg/1000000, numjny/1000000,  round(groupProb[to]*100,digits=1), kld)
    names(summary)[1:2] = c("Psg(million)","Jny(million)")
    names(summary)[8] = "KLD"
    return ( summary )
  }
}

{
    psggrouptable = rbind(
      psggroupsummary(ShortburstPsg),
      psggroupsummary(OccasionalPsg),
      psggroupsummary(CommuterDrift),
      psggroupsummary(CommuterLoyal),
      psggroupsummary(Other), 
      psggroupsummary (1:NP))
    
    psggrouptable=as.data.table(psggrouptable)

    rownames(psggrouptable)=c("Short Burst","Occasional",
                              "Drifting Commuter","Loyal Commuter",
                              "Others",
                              "Total")
    
    write.csv(psggrouptable,"./results/passengerGroupTable.csv",row.names=TRUE)

    #for paper with rounding
    psggrouptable[,1:2] = round(psggrouptable[,1:2],digits=3)
    psggrouptable[,3:7] = round(psggrouptable[,3:7])

    
    xtable(psggrouptable,
           digits=c(0,3,1, 0,0,0,0,0, 3),
           label="tab:psggroups",
           caption="Passenger subgroups metrics")
}

## create sample passenger motivation figure for paper 
## using randomly selected passengers from the 1.72 million
{
  showpassengers <- function(psglist)
  {
    np=length(psglist)
    par(mfrow=c(np,1),mai=c(0.4,0.8,0.3,0.3))
    for (i in 1:np) {
      p = psglist[i]
      tklist = setdiff(unique(tokenGroupMatrixMonthly[p,]),"Abs")
      pname = paste(paste("P",i,sep=""),paste(sapply(tklist,tokengetfullname),collapse="-"),sep=" ") #multi cards
      pj = journeyMatrix[p,]
      tk = unique(tokenGroupMatrixMonthly[p,])
      tk = tk[which(tk!="Abs")][1]
      tkcol = tokencols[which(tokennames==tk)]
      bb=barplot(pj,col=tkcol,ylab="Journeys/month",ylim=c(0,30)) #col=freqbandcols[pc+1]
      abline(v=bb[(1:5)*12]+0.5,col="gray") #show year boundaries
      #axis(4)
      legend("topleft",pname, bg="white") 
    }
    #add date axis
    axis(1,at=bb[(1:5)*12-6],labels=2015:2019,tick=FALSE)
    par(mfrow=c(1,1))
  }
  
  pdf(file="./results/examplepassengers.pdf",width=9,height=9)
  showpassengers(c(815691,889945,1347177,305802,1056134))
  dev.off()
}
