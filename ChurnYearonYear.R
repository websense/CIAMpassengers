## Passenger Churn model
## Author: Rachel Cardell-Oliver
## Version: V1 2021 V2 May 2022

#libraries and constants
{
  source("./Globalconstants.R") #load datasets, define utility fns and journey band definitions
}

{
peryearjourneys = 
  cbind(rowSums(journeyMatrix[,1:12]),
        rowSums(journeyMatrix[,12+1:12]),
        rowSums(journeyMatrix[,24+1:12]),
        rowSums(journeyMatrix[,36+1:12]),
        rowSums(journeyMatrix[,48+1:12]))

}

## active psgs each year candidate passengers for churn 2015-17, 2016-18, 2016-19
{
  isactive2015 = which(peryearjourneys[,1]>0)
  isactive2016 = which(peryearjourneys[,2]>0)
  isactive2017 = which(peryearjourneys[,3]>0)
  isactive2018 = which(peryearjourneys[,4]>0)
  isactive2019 = which(peryearjourneys[,5]>0)
  
  datasetChurn <- as.data.table(rbind(
    cbind(2016, peryearjourneys[isactive2016,1:3]),
    cbind(2017, peryearjourneys[isactive2017,2:4]),
    cbind(2018, peryearjourneys[isactive2018,3:5])))
  setnames(datasetChurn,c("ReferenceYear","JnyLast","JnyThis","JnyNext"))
}

{
  #Q: what proportion of all psgs appear at least once in the churn dataset?   A: 77.52
  uniquepsgChurn = length(Reduce(union, list(isactive2016,isactive2017,isactive2018)))
  uniquepsgChurnpct = round(uniquepsgChurn*100/dim(journeyMatrix)[1],digits=2) 
}

## Churn calculations
## candidate passengers for churn 2015-17, 2016-18, 2016-19 see CurateAnalysisDataSets
{
  

  churnMatrix <- datasetChurn
  churnMatrix[JnyLast == 0 & JnyThis > 0 & JnyNext > 0, gained:=TRUE]
  churnMatrix[JnyLast > 0 & JnyThis > 0 & JnyNext > 0, retained:=TRUE]
  churnMatrix[JnyLast > 0 & JnyThis > 0 & JnyNext == 0, lost:=TRUE]
  churnMatrix[JnyLast == 0 & JnyThis > 0 & JnyNext == 0, temporary:=TRUE]

}

## Q: does churn change over time? A: not much
{
  churnCountsByRefYear =
  rbind(colSums(churnMatrix[ReferenceYear==2016,.(gained,retained,lost, temporary)],na.rm=TRUE),
  colSums(churnMatrix[ReferenceYear==2017,.(gained,retained,lost, temporary)],na.rm=TRUE),
  colSums(churnMatrix[ReferenceYear==2018,.(gained,retained,lost, temporary)],na.rm=TRUE),
  colSums(churnMatrix[,.(gained,retained,lost, temporary)],na.rm=TRUE))
  
  ## very close all years, so just report overall figures for all years
  round(churnCountsByRefYear *100/rowSums(churnCountsByRefYear ),digits=1)

}

## Tokens and Card counts 
{
  #returns number of passengers c(gained, retained, lost, temporary, inactive) relative to reference year from start month sm:(sm+12-1)
  getchurncounts <- function(sm, journeyMatrix, tokenGroupMatrixMonthly, token="All")
  {
    
    if ((sm-12 < 1)|(sm+24-1 > 60))
    {
      print("ERROR: sm-12-1 to sm+24-1 are not included in 1:60. Must have all 3 years included")
      return(NULL)
    }
    
    ## ranges for the reference year and succ and pred
    baseyear = sm:(sm+12-1)
    prevyear = (sm-12):(sm-1)
    nextyear = (sm+12):(sm+24-1)
    #get records for token type, counting all passg who used that token in the base year
    #there will be some double counting
    if (token != "All") { 
      #some double counting when passg token changes
      tks = which(apply(tokenGroupMatrixMonthly[,baseyear],1,
                        function(ss) { return(is.element(token,ss))} ))
    } else {
      tks = 1:dim(tokenGroupMatrixMonthly)[1] #all
    }
    
    #number of journeys made per token passenger for the 3 churn ref years
    lastyear = rowSums(journeyMatrix[tks,prevyear])
    thisyear = rowSums(journeyMatrix[tks,baseyear])
    nextyear = rowSums(journeyMatrix[tks,nextyear])
    
    gained = sum(lastyear == 0 & thisyear>0 & nextyear>0)
    retained = sum(lastyear>0 & thisyear>0 & nextyear>0)
    lost = sum(lastyear>0 & thisyear>0 & nextyear == 0)
    temporary = sum(lastyear == 0 & thisyear>0 & nextyear == 0)
    inactive = sum(thisyear == 0) #not active in reference year should be 0 for tokens
    allactive = sum(c(gained, retained, lost, temporary))
    return( c(gained, retained, lost, temporary, inactive) ) 
    
  }
  
  res = c()
  for (tk in c(tokennames,"All")) {
    restk = colSums(
      rbind(
      getchurncounts(13, journeyMatrix, tokenGroupMatrixMonthly, tk),
      getchurncounts(25, journeyMatrix, tokenGroupMatrixMonthly, tk),
      getchurncounts(37, journeyMatrix, tokenGroupMatrixMonthly, tk)))
    res = rbind(res, restk)
  }
  
  rownames(res) = c(tokenfullnames,"Population")
  colnames(res) = c("gained", "retained", "lost", "temporary","inactive")
  churnPercentByToken = res[,1:4]*100/rowSums(res[,1:4])
  churnPercentByToken = cbind(churnPercentByToken,100)
  
  ## use KLD to check for surprising churn categories
  round(c(
    KLdiv(churnPercentByToken[1,1:4]/100,churnPercentByToken[6,1:4]/100), #diff churn class to popn
    KLdiv(churnPercentByToken[2,1:4]/100,churnPercentByToken[6,1:4]/100), 
    KLdiv(churnPercentByToken[3,1:4]/100,churnPercentByToken[6,1:4]/100), 
    KLdiv(churnPercentByToken[4,1:4]/100,churnPercentByToken[6,1:4]/100),
    KLdiv(churnPercentByToken[5,1:4]/100,churnPercentByToken[6,1:4]/100),
    KLdiv(churnPercentByToken[6,1:4]/100,churnPercentByToken[6,1:4]/100)),digits=3)
  #0.003 0.049 0.016 0.027 0.019 0.000  nothing too unexpected (in token mix)
  
  #Table 2 in paper: churn by token type
  write.csv( round(churnPercentByToken,digits=3), "./results/churnPercentByToken.csv")
}


## END 
