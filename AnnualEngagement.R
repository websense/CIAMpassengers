## Annual Engagement Clustering
## and analysis of clustering quality
## Author: Rachel Cardell-Oliver
## Versions: V1 2020, Sep 2021, Revised May 2022


#libraries and constants
{
  source("./GlobalConstants.R") 
  
  TESTCL=2:12  #number of clusters to test
  BESTCL=5  #best number of clusters (for full dataset) chosen by cluster quality analysis
}

## get data (created by travelBands) and organise as per Intensity group for clustering
{
peryearjourneys = 
        c(rowSums(journeyMatrix[,1:12]),
          rowSums(journeyMatrix[,12+1:12]),
          rowSums(journeyMatrix[,24+1:12]),
          rowSums(journeyMatrix[,36+1:12]),
          rowSums(journeyMatrix[,48+1:12]))
  
peryeartokens = 
    rbind(tokenGroupMatrixMonthly[,1:12],
      tokenGroupMatrixMonthly[,12+1:12],
      tokenGroupMatrixMonthly[,24+1:12],
      tokenGroupMatrixMonthly[,36+1:12],
      tokenGroupMatrixMonthly[,48+1:12])

passengerClassesMonthlyByYear = readRDS("./data/passengerClassesMonthlyByYear.rds") #created by IntensityMonthly.R

#only cluster on nonZero years
nonZeroyears = which(peryearjourneys > 0) #no need to cluster inactive years (they are E0)
length(nonZeroyears)/length(peryearjourneys) #0.5139591 51% of all possible from 5 years
    
#make entropy vectors of distribution of classes for clustering
Nyears = dim(passengerClassesMonthlyByYear)[1] #9838285
maxF = max(passengerClassesMonthlyByYear) #max band label = 5
entropiesY = matrix(data=0,ncol=length(freqbandnames),nrow=Nyears) 
for (c in 0:maxF) {
  entropiesY[,c+1] <- rowSums(passengerClassesMonthlyByYear==c,na.rm=TRUE)
}
colnames(entropiesY) = freqbandnames

saveRDS(entropiesY,file="./data/entropiesY.rds")
}

## function to cluster the yearly data with chosen NC=BESTCL
#clusters sample profiles in entropiesY into NC clusters and saves centres and labels
#outputs the files clustering/AEDclustering5.rds clustering/entropy5CentersSize.rds and clustering/entropies5All.rds
{
makeAndSaveAEclustering <- function(entropiesY, NC=BESTCL)
{
    set.seed(9876) #seed for repeatable results
    nonZeroyears = which(peryearjourneys > 0)
    kmY = kmeans(entropiesY[nonZeroyears,],NC,nstart=10) 
    saveRDS(kmY,file=sprintf("./clustering/AEDclustering%d.rds",NC))

    #order engagement clusters by size of cluster
    oo=order(kmY$size,decreasing=TRUE)
    xx=cbind(kmY$centers, kmY$size)[oo,]
    colnames(xx)[7]=c("ClusterSize(years)") 
    rownames(xx)=paste("E",1:NC,sep="") #clusterID
    entropyCentersSize = xx #already ordered
    

    #cluster IDs rename so they are in size order
    entropyCluster = kmY$cluster*10
    for (cid in 1:NC) {
      entropyCluster[which(kmY$cluster==oo[cid])] = cid 
    }
    #save labels for all years (used later for transition probs)
    entropiesZ = cbind(entropiesY,rowSums(passengerClassesMonthlyByYear!=0),peryearjourneys,0) 
    entropiesZ[nonZeroyears,9] = entropyCluster #assign learned clusters to nonZero years
    colnames(entropiesZ)[7:9]=c("ActiveMonths","Journeys","entropyCluster")

    saveRDS(entropyCentersSize,file=sprintf("./clustering/entropy%dCentersSize.rds",NC))
    saveRDS(entropiesZ,file=sprintf("./clustering/entropies%dAll.rds",NC))
  }
}


# SAVE FILE show clusters as bar charts for BESTCL clusters
{
  makeAndSaveAEclustering(entropiesY,BESTCL)
  entropyCentersSize = readRDS(sprintf("./clustering/entropy%dCentersSize.rds",BESTCL))
  
  pdf(file=sprintf("./results/Entropy%dClustersBarChart.pdf",BESTCL),height=5.5,width=6)  
  {
  par(mfrow=c(1,1))
  nb=6 #number of bands
  bp=barplot(t(entropyCentersSize[,1:nb]),
             yaxt="n",xaxt="n",
             ylab="Travel Months per Year",
             ylim=c(0,16), #legend space
             xlab = "Annual Engagement Distribution Clusters", 
             col=freqbandcols)

  axis(1,at=bp,labels=paste("E",1:BESTCL,sep=""),
       tick=FALSE,line=-.5)
  axis(2,at=0:12,labels=TRUE)  #months in year
  axis(4,at=0:12,labels=TRUE) 
   legend("top",title="Travel Frequency Bands (see Table 3)",
          ncol=3,bg="white",
          legend=freqbandnames,
          fill=freqbandcols)
  par(mfrow=c(1,1))
}
  dev.off()

}


## Table of AE class metrics
## Get token data and generate table of properties of the clusters
## Saves entropyTable5.csv in clustering directory
{
  makeEntropyTable <- function(nc=BESTCL)
  {
    entropiesZ <- readRDS(file=sprintf("./clustering/entropies%dAll.rds",nc))
    entropyCluster = entropiesZ[,9]
    
    #reordering of token names from alpha "Con" "Sch" "Snr" "Std" "Ter" 
    to = c(4,2,3,5,1)  #token order Std,Sch,Snr,Ter,Con
    entropyTable = matrix(0,nrow=(nc+1),ncol=8)
    rownames(entropyTable)=c(paste("E",1:nc,sep=""),"Total")
    colnames(entropyTable) = c("ClusterSize","AvgJnyPerYear",tokenfullnames,"KLD")

    # get per token probabilities per cluster
    for (ae in 1:(nc+1)) 
    {
       if (ae<=NC) {
        thiscluster = which(entropyCluster==ae)
       } else {
        thiscluster = which(entropyCluster!=0) #all engaged
      }
       
       clsize = round(length(thiscluster)/1000000,digits=2) #cl size in millions
       avgjnyperyear = round(mean(entropiesZ[thiscluster,8]))
       gp = table(peryeartokens[thiscluster,])
       groupProb = gp[2:6]/sum(gp[2:6]) #all except Abs
       if (sum(names(groupProb)==names(tokenProbabilities))!=5) {
         print(paste(c("WARNING unexpected token name order nc=",nc,"cl=",ae,names(groupProb),
                       " != ",names(tokenProbabilities)),collapse=" "))
         print(groupProb)
       }
       kld = round(KLdiv(groupProb,tokenProbabilities),digits=3)
  
       entropyTable[ae,] = c(clsize, avgjnyperyear, round(groupProb[to]*100,digits=3), kld)
    }
  
    write.csv(entropyTable,sprintf("./clustering/entropyTable%d.csv",nc),row.names=TRUE)
  }
}



### Analysis of clustering quality
{
## generate multiple versions with different numbers of clusters testing which number of clusters is best
  for (nc in TESTCL) {
    #makeAndSaveAEclustering(entropiesY,nc)
    makeEntropyTable(nc)
  }

  ## stability test
  # samplepct in 0.1 to 1.0 how much of the full dataset to select
  # this is pretty slow - reclustering + ARI is slow to calc
  # given nc eg=5 read previously saved entropies5All which contains engagement profiles and assigned cluster ID
  # cluster output for this number of clusters
  clusteringRobustness <- function(nc,repeats=5,samplepct=0.5) 
  {
    entropiesZ = readRDS(sprintf("./clustering/entropies%dAll.rds",nc))
    nonZero = which(entropiesZ[,9]>0) #ignore E0
    clustering = entropiesZ[nonZero,9]
    size = length(clustering)
    samplesize = round(size*samplepct)
    similarity=vector(length=repeats)
    for (i in 1:repeats) { #repeat reclustering 5 times
      #take a (new) sample each time
      ss = sample(1:size,samplesize) 
      cl1 = clustering[ss] #assigned clusters at ss positions
      #re-cluster test with % of data
      kmTest = kmeans(entropiesZ[nonZero[ss],1:6],nc,nstart=5)  #cluster the profiles (only)
      cl2 = kmTest$cluster  #cluster assignments should be similar for same observations
      similarity[i] = adjustedRandIndex(cl1,cl2)
    }
    return(similarity)
  }
}


## run cluster quality tests
  {
## TEST 1: measure CLUSTERING  QUALITY results for variance fig
  # creates and saves clusteringquality matrix
  {
    clusteringquality = matrix(NA,nrow=length(TESTCL),ncol=6)
    clusteringquality[,1] = TESTCL #number of clusters
    for (nc in TESTCL) {
      row = nc-1
      kmY = readRDS(file=sprintf("./clustering/AEDclustering%d.rds",nc))
      clusteringquality[row,2] = round(kmY$betweenss / kmY$totss,digits=3)
      ## is between cluster sum of squares (separation)
      #A slight variation of this method plots the curvature of the within group variance
      #kmY$tot.withinss - always decreases as NC inc, but gains get smaller 
      
      entropyTable =  read.csv(sprintf("./clustering/entropyTable%d.csv",nc))[,2:9]
      
      ## range of passenger years per cluster (cluster size)
      psgrr = range( entropyTable[1:nc,1]/sum(entropyTable[1:nc,1]))
      clusteringquality[row,3:4] = psgrr

      jnyrr = range( entropyTable[1:nc,2]/sum(entropyTable[1:nc,2]))
      clusteringquality[row,5:6] = jnyrr
      
    }
    colnames(clusteringquality)=c("NC","VarExplained","MinClSize","MaxClSize","MinClJny","MaxClJny")
    write.csv(clusteringquality,file="./results/clusteringquality.csv",row.names=FALSE)
  }

  ## TEST 2: Robustness - same clusters under multiple clusterings
  # very slow so save result as we go
  # 8 did not converge warnings out of 12*3*5*5 clusterings
  {
    repeats = 5 
    robustness=c()
    for (samplepct in c(0.25,0.5,1.0)) 
      {
      for (nc in TESTCL) {
        sampleres =   clusteringRobustness(nc,repeats,samplepct) 
        robustness = rbind(robustness,
                           c(nc,repeats,samplepct,sampleres))
        #slow so save on the way
        colnames(robustness)=c("NC","Repeats","Sample%",paste("ARI",1:repeats,sep=""))
        write.csv(robustness,file="./results/robustness.csv",row.names=FALSE)
      }
    }
 }
   

  pdf(file="./results/robustnessToSampling.pdf",width=10,height=10)
  {
    par(mfrow=c(2,2))
    boxplot(t(robustness[which(robustness[,3]==0.25),4:8]),main="Sample 0.25",names=TESTCL,ylim=c(0.6,1))
    boxplot(t(robustness[which(robustness[,3]==0.5),4:8]),main="Sample 0.5",names=TESTCL,ylim=c(0.6,1))
    boxplot(t(robustness[which(robustness[,3]==1.0),4:8]),main="Sample all",names=TESTCL,ylim=c(0.6,1))
    boxplot(c(),main="All samples",xlim=c(2,12),ylim=c(0.6,1),xlab="Number of Clusters",ylab="Adjusted Rand Index")
    for (c in TESTCL) {
      boxplot(c(robustness[which(robustness[,1]==c),4:8]),at=c,add=TRUE)
    }
    axis(1,at=TESTCL,TESTCL)
    par(mfrow=c(1,1))
  }
  dev.off()
  
  ## TEST 3: Distinctiveness: distribution of KLD over all clusters
  {  
     distinctiveness = matrix(NA,nrow=length(TESTCL),ncol=max(TESTCL)+1)
     distinctiveness[,1] = TESTCL #number of clusters
     for (nc in TESTCL) {
      entropyTable =  read.csv(sprintf("./clustering/entropyTable%d.csv",nc))
      distinctiveness[nc-1,2:(nc+1)] = entropyTable[1:nc,9] #KLD column
     }
     write.csv(distinctiveness,file="./results/distinctiveness.csv",row.names=FALSE)
  }
}



# FIGURE three metrics used to determine best number of clusters for AE
# BEST is trade-off between variance explained and distinctiveness, 
# which both improve with more clusters
# and stabilty of the clusterings (to diff data and repeated cls) which starts failing after 5
{
  pdf(file="./results/clusteringEvaluation.pdf",width=15,height=5)
  par(mfrow=c(1,3))
  {
    {
      rr=TESTCL
      #var explained higher is better
      plot(rr,clusteringquality[rr-1,2], #var explained
           type="b",pch=1, lwd=2,col="blue",
           ylim=c(0,1),
           xlab="Number of Clusters",ylab="Variance, Cluster Size (max,min)",
           main="Variance Explained by Clustering")
      grid()
      # cluster size, lower max and higher min are better - avoid extreme unbalance
      lines(rr,clusteringquality[rr-1,3],lty="dotted",type="b",pch="P") #psg
      lines(rr,clusteringquality[rr-1,4],lty="dashed",type="b",pch="P") #psg
      lines(rr,clusteringquality[rr-1,5],lty="dotted",type="b",pch="J") #jny
      lines(rr,clusteringquality[rr-1,5],lty="dashed",type="b",pch="J") #jny
      abline(v=BESTCL,lty=1,col="red")
    }
    
    # ARI similarity for different training data, higher is better (more robust)
    # show all as crosses (boxplots misleading)
    ## rand index compared to base is fine
    ## but need to normalise for K by comparing with a random sample?
    ## I think ARI does that
    {
      boxplot(c(), main="Stability of Clusterings",xlim=c(2,12),ylim=c(0.6,1),xlab="Number of Clusters",ylab="Adjusted Rand Index")
      for (c in TESTCL) {
        allaris = c(robustness[which(robustness[,1]==c),4:8])
        boxplot(allaris,at=c,add=TRUE,border="darkgray",range=0,width=1)
        lines(rep(c,times=length(allaris)),allaris,yaxt="n",
              type="p",pch=4,col="blue")
      }
      axis(1,at=TESTCL,TESTCL)
      grid(nx=0,ny=NULL)
      abline(v=BESTCL,lty=c(1),col="red")
    }
    
  
    #KLdiv divergence from main population increases with number of clusters, higher is better
    # show all as crosses (boxplots misleading)
    {
      rr=TESTCL
      klds = distinctiveness[,2:13]
      boxplot(t(klds),ylim=c(0,max(klds,na.rm=TRUE)),border="darkgray",
              at=TESTCL,names=TESTCL,range=0,
              xlab="Number of Clusters",ylab="KLdiv clusters vs population",
              main="Distinctiveness of Categories per Cluster")
      for (nc in TESTCL) {
        lines(rep(nc,times=nc),klds[nc-1,1:nc],type="p",pch=4,col="blue") #crosses
      }
      grid(nx=0,ny=NULL)
      abline(v=BESTCL,lty=c(1),col="red")
    }

  }
  par(mfrow=c(1,1))
  dev.off()
} 

## END

