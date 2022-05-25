## Multi-Year Engagement Results
## Author: Rachel Cardell-Oliver
## Version: V1 2020, Sep 2021, Revised May 2022

#libraries and constants
{
  source("./GlobalConstants.R")
}

## Prepare data with AE cluster per year and token per year
{
  BESTCL=5
  entropiesZ <- readRDS(file=sprintf("./clustering/entropies%dAll.rds",BESTCL))
  NY = dim(entropiesZ)[1]/5 #number of passenger records for 1 year

  # 5 year passenger records
  y1range = 1:NY
  y2range = (NY*1) + y1range 
  y3range = (NY*2) + y1range 
  y4range = (NY*3) + y1range 
  y5range = (NY*4) + y1range 
  
  #get AE class per year
  engagePerYear = cbind(entropiesZ[y1range,9],
                        entropiesZ[y2range,9],
                        entropiesZ[y3range,9],
                        entropiesZ[y4range,9],
                        entropiesZ[y5range,9])
  saveRDS(engagePerYear,"./data/engagePerYear.rds")

  
  ## tokenperyear[a,b,c] is TRUE if tokennames[a] for passenger b in year 2014+c is active
  tokenperyear = array(data=FALSE, dim=c(5, NY, 5),dimnames=list(tokennames, NULL, 2015:2019))
  
  for (i in 1:5 ) {
    tk = tokennames[i]
    for (y in 1:5) {
       tokenperyear[i, ,y] = rowSums(tokenGroupMatrixMonthly[,(y-1)*12 + 1:12]==tk)>0
    }
  }

}


## GENERATE Markov year to year transitions between AE clusters
{
  # input yearStep offset from Base year to next (eg +2 is 2 years later)
  # uses engageperyear and sampleTokenGroupMatrixYearly
  # returns engageChangeProbs = array(dim=c(5,6,6) Ei, Ej, prob by each card
  learnMarkov <- function(yearStep)
  {
    engageChangeProbs = array(dim=c(5,6,5),dimnames=list(paste("E",1:5,sep=""),
                                                         paste("E",0:5,sep=""),
                                                         tokennames) )
    ## calculate transition probabilities from year i to year i+1
    for (ti in 1:5) {
      tok = tokennames[ti]
      for (e1 in 1:5) {
        for (e2 in 0:5) {
          ## excluding 0 to 0: must be active in y1
          ny12=0
          for (y in 1:(5-yearStep)) { #get number of passengers of type tok changing from e1 to e2 engagement from year A to B
            ny12 = ny12+sum(engagePerYear[,y]==e1 & engagePerYear[,(y+yearStep)]==e2 & tokenperyear[ti,,y])
          }
          ny1=0
          for (y in 1:(5-yearStep)) { #get base count in year A
            ny1 = ny1+sum(engagePerYear[,y]==e1 & tokenperyear[ti,,y])
          }
          #Pr of e2 next year give e1 this year
          engageChangeProbs[e1,e2+1,ti] = ny12 / ny1
        }
      }
    }
    
    return(engageChangeProbs)
  }
  
  # plot markov results
  #String for year change (eg Plus1)
  #Strings for y and x axes eg "Base Year" "Next Year"
  #tofile saves to figdir eg ./markovfigs/
  plotMarkov <- function(engageChangeProbs,titleY,titleY1,titleY2,
                         tofile=FALSE,figdir="./markovfigs/")
  {
    for (ti in 1:5) 
    {
      tkname = tokenfullnames[ti]
      if (tofile) {
        pdf(file = sprintf("%s%s%sMarkovHeatMap.pdf",figdir,tkname,titleY),
            width = 5,height = 5)
      } #else plot to screen
      {
        heatmap.2(
          engageChangeProbs[, , ti], 
          cellnote = round(engageChangeProbs[, , ti]*100), #show as pct
          Rowv = NULL,
          Colv = NULL,
          dendrogram = "none", #don't draw dendogram or reorder columns
          margins = c(5,5), #reduce margins
          col = brewer.pal(6, "YlOrRd"),
          notecol = "black",
          notecex=1.5, #1.25, #larger cell prob labels
          density.info = "none", # turns off density plot inside color legend
          trace = "none",
          key = FALSE, #don't show legend for colours
          ylab = titleY1,
          xlab = titleY2)
        mtext(tkname,3,line=-2) #title at lines counted from centre 
      }
      if (tofile) { dev.off() }
    }
  }
}

## GENERATE FIGS ggplot labelled heat maps
{

  oneyearChangeProbs = learnMarkov(1)
  twoyearChangeProbs = learnMarkov(2)  
  threeyearChangeProbs = learnMarkov(3)
  fouryearChangeProbs = learnMarkov(4)
  
  
  plotMarkov(oneyearChangeProbs,"OneYear","Engagement in base year","Engagement in base + 1 year",
             tofile=TRUE,figdir="./markovfigs/")
  plotMarkov(twoyearChangeProbs,"TwoYear","Engagement in base year","Engagement in base + 2 years",
             tofile=TRUE,figdir="./markovfigs/")
  plotMarkov(threeyearChangeProbs,"ThreeYear","Engagement in base year","Engagement in base + 3 years",
             tofile=TRUE,figdir="./markovfigs/")
  plotMarkov(fouryearChangeProbs,"FourYear","Engagement in base year","Engagement in base + 4 years",
             tofile=TRUE,figdir="./markovfigs/")
  
}


## END


