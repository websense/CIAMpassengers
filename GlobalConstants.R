## Definitions, constants and functions required for analysis
## Load data files, set common names and colour constants and useful metrics
## Author: Rachel Cardell-Oliver
## Version: 2020, 2021 Revised May 2022

{
  require(RColorBrewer)
  require(xtable)
  require(data.table)
  require(mclust) #for adjustedRandIndex
  require(gplots) #for heatmap.2
}

## pre-processed data used for all analyses
## the datasets for the paper are used under a privacy agreement so can not be published
## please see below for the format required and substitute other data sets
{
  journeyMatrix = readRDS("./data/journeyMatrix.rds")
  #Format: matrix N=1,722,120 passengers x 60 cells journey counts
  # sample row: 0     2     0     0     3    22     4     6    13 ... 0
  
  tokenGroupMatrixMonthly = readRDS("./data/tokenGroupMatrixMonthly.rds")
  #Format: matrix N passengers x 60 cells  Abs  Con  Day  Sch  Snr  Std  Ter 
  # sample row: "Std" "Std" "Std" "Abs" "Std" ... "Ter"
  
  ## token probabilities in the full dataset - use this in other files
  tokentable = table(tokenGroupMatrixMonthly)
  tokenProbabilities = tokentable[2:6]/sum(tokentable[2:6])
  write.csv(tokenProbabilities,"./results/tokenProbabilities.csv")
}

## Freq Bands info
{
  ## from IntensityMonthly.R experiments 
  freqbounds =  c(1, 3, 17, 33,  41, 210)
  freqbandnamesdescriptive = c("Zero","Rare","Low","Medium","High","VeryHigh")
  freqbandnames = paste("I",0:5,sep="")
  #per month journeys for 4 week months
  freqbandshort = c("0..0","1..2","3..16","17..32","34..40","41+")
  freqbandcols =  brewer.pal(6,"RdYlBu")[6:1] #same cols as in engagement clusters
  
  #given a journeys per month count, return intensity band
  getTravelBand <- function(v)  
  {
    b = freqbounds
    c <- rep(NA,times=length(v))   #NA for inactive years
    c[which(v==0)] = 0 #v==0 no travel
    for (i in 1:(length(b)-1)) {
      c[which(v>=b[i] & v<b[i+1])] = i
    }
    return(c)
  }

  engagementcols =   brewer.pal(6,"YlOrBr") #E0 to E5, like a heat map
}

## Passenger category (token) information
{

  tokencols = brewer.pal(5,"Set1")[c(1,3,2,5,4)]
  tokennames = c("Std", "Sch", "Snr","Ter", "Con" ) # (not Abs)
  tokencols = brewer.pal(5,"Set1")
  
  tokenfullnames = c("Standard","School", "Senior", "Tertiary", "Concession") #order for paper tables
  
  #used in plotpassenger for showing card changes
  tokengetfullname <- function(shortname) {
    return(tokenfullnames[which(tokennames==shortname)])
  }
  
}




## Functions for calculating useful clustering and other metrics
{
  #shannon entropy H function
  entropy <- function(v) {
    pi = v / sum(v)  #probability of each symbol
    logpi = log(pi,base=2)
    #%In the case of P(xi) = 0 for some i, the value 
    #of the corresponding summand 0 logb(0) is taken to be 0 
    #https://en.wikipedia.org/wiki/Entropy_(information_theory)
    logpi[which(v==0)]=0
    return (-1 * sum(pi*logpi))
  }
  
  #KL divergence Kullback-Leibler Divergence for small prob vectors
  # https://www.countbayesie.com/blog/2017/5/9/kullback-leibler-divergence-explained
  # eg x is prob of token type for one gp and y for another
  # KLdiv is not symmetric (divergence not a distance) 
  ## KLdiv (Observed ∣∣ ApproximatedorExpected) 
  ## so y should be the "base" eg population probs 
  ## compared with observed in some subgroup
  KLdiv <- function(observed,expected)
  {
    if (sum(is.nan(observed))>0 | sum(is.nan(expected))>0) { 
      print("ERROR: KLdiv observed and y must be probability vectors not NaNs")
      return(NA) 
    }
    #http://uc-r.github.io/comparing_numeric_values/
    if (!all.equal(sum(observed),1.0) | !all.equal(sum(expected),1.0) | length(observed)!=length(expected) ) {
      print("ERROR: x and y must be probability vectors (sum to 1)  and have same length")
      return (NA) 
    }
    lxy = log2(observed/expected)
    lxy[which(lxy==-Inf)]=0  #map 0 probs to 0 in sum
    return( sum(observed*lxy) )
  }
  
  
## Other cluster quality metrics, not included in paper 
## but were used for testing and maybe useful for future analysis
{ 
  ## kulczynski measure from Han p270
  # K(P,T) = 0.5 ( p(P|T) + P(T|P) )
  # 0 means independence, 1 means A and B always occur together
  #input: confusion matrix for all passenger group and token counts
  # return all Kulc measures for each pi and ti
  #pi indicates which passenger group (col eg 1:9)
  #ti indicates which token (row eg 1:4 or 1:7)
  #number of digits returned, default 2
  kulc <- function( confmat, dig=2 ) 
  {
    ncols = dim(confmat)[2]
    nrows = dim(confmat)[1]
    km = confmat*0.0 #matrix of kulc values
    for (c in 1:ncols) {
      countB = sum(confmat[,c])
      for (r in 1:nrows) {
        countA = sum(confmat[r,])
        countAB = confmat[r,c]
        
        probAgivenB = countAB / countB
        probBgivenA = countAB / countA
        
        km[r,c] =  0.5 * (probAgivenB + probBgivenA)
        #max conf max(probPgivenT , probTgivenP)
        #cosine sqrt(probPgivenT * probTgivenP)
        #kulc 0.5 * (probPgivenT + probTgivenP)
      }
    }
    return(round(km,digits=dig))
  }
  
  
  ## imbalance ratio
  ## range [0,1] 0 is perfectly balanced and 1 is imbalanced
  #number of digits returned, default 2
  ### error: IR(a,b) = IR(b,a)
  IR <- function(confmat, dig=2)
  {
    ncols = dim(confmat)[2]
    nrows = dim(confmat)[1]
    ir = confmat*0.0 #matrix of kulc values
    for (c in 1:ncols) {
      countB = sum(confmat[,c])
      for (r in 1:nrows) {
        countA = sum(confmat[r,])
        countAB = confmat[r,c] #both
        
        #print(c(c,r,countA,countB,countAB))
        ir[r,c] = (abs(countA - countB) / (countA + countB - countAB))
      }
    }
    return (round(ir,digits=2))
  }
  
  
  #Jensen-Shannon Divergence 
  # distance between 2 discrete prob distributions
  # https://stats.stackexchange.com/questions/208118/distance-measure-of-two-discrete-probability-histograms-distance-between-two-ve
  # Jensen-Shannon distance is the 1st thing I'd consider. If you don't insist on having a "distance function", you can directly use Jensen–Shannon divergence, from which this distance is derived.
  # JSD is based on the Kullback–Leibler divergence, with some notable (and useful) differences, including that it is symmetric and it always has a finite value. 
  ### The square root of the Jensen–Shannon divergence is a metric often referred to as Jensen-Shannon distance
  # https://en.wikipedia.org/wiki/Statistical_distance
  # Distance must be 1) d(x,y)>=0 2) d(x,x)=0 3) dxy=dyx 4) dxz<=dxy+dyz triangle 
  # just 1 and 2 is a divegence
  # JS divergence is widely used to measure the difference between two probability distributions. 
  # It fits your case, as the inputs are two probability vectors. JS divergence is a straightforward modification of the well-known Kullback–Leibler divergence.
  # Generally, KL and JS divergence require the input vectors have nonzero entries. 
  # In case of zeros in the input, many people simply choose to throw out those values. 
  # Check https://mathoverflow.net/a/72672 for more details on this issue.
  # 
  #https://cran.r-project.org/web/packages/philentropy/vignettes/Information_Theory.html
  # philentropy
  #This function computes the Jensen-Shannon Divergence JSD(P || Q) between two probability distributions P and Q with equal weights π_1 = π_2 = 1/2.
  #The Jensen-Shannon Divergence JSD(P || Q) between two probability distributions P and Q is defined as:
  #   
  #   JSD(P||Q)=0.5∗(KL(P||R)+KL(Q||R))
  # 
  # where R = 0.5 * (P + Q) denotes the mid-point of the probability vectors P and Q, 
  # and KL(P || R), KL(Q || R) denote the Kullback-Leibler Divergence of P and R, as well as Q and R.
  # The square root of the Jensen–Shannon divergence is a metric often referred to as Jensen-Shannon distance
  # TODO THIS VERSION IS Not symmetric
  # JS distance for small prob vectors
  # eg x is prob of token type for one gp and y for another
  JSD <- function(x,y)
  {
    z = 0.5 * (x+y) 
    jsd = sqrt(0.5 * (KLdiv(x,z) + KLdiv(y,z)))
    return(jsd)
  }
  
  intermittency <- function(ts) {
    rr = rle(sort(ts))
    mm = as.integer(rr$values[which.max(rr$lengths)[1]])
    ## probability that t_i is the mode value
    return( sum(ts==mm) / length(ts) )
  }
  
  #Histogram distance
  #earth mover distance: 12 0 0 0 to 0 12 0 0 is closer than 0 0 0 12
  #non-normalised version - number of bins must be the same in both
  #hc1 and hc2 are AED vectors
  ## Bin size b $$MDPA(A,B,b) = cumSum ( abs(A_n  - B_n) \times b ) / N$$
  ## From Cha2002
  ## Histograms A and B (normalised) with b bins
  # Algorithm: Distance-ordinal-histogram(int∗A; int∗B)
  # 1 prefixsum = 0
  # 2 hdist = 0 
  # 3 for i =0 to b − 1 
  # 4     prefixsum += A[i] − B[i]
  # 5     hdist += |prefixsum|
  # 6 return(hdist)
  # since bins are uneven?? OR do this for raw journey counts
  # option allow bin medians and then weight the transfers
  # default null is not to do it, else take into account the sum
  # for travel bands binmedians = c(0,2,7,23,35,43)
  ## STOP binmedians does not quite work - use jny cnt if this is the real aim
  MDPAdist <- function(h1c,h2c) #,binmedians=NULL)
  {
    # n1 = sum(h1c) #sample size for hist 1,2
    # n2 = sum(h2c)
    # if (abs(n1-n2)>1) { #allow for leap years +- 1 day OK
    #   print(sprintf("ERROR: number of obs %d and %d must be same for MDPAdist",n1,n2))
    # }
    nb = length(h1c) #number of bins
    # if (!is.null(binmedians)) {
    #   bindists = (binmedians[2:nb] - binmedians[1:(nb-1)])
    # } else {
    #   bindists = c(1,rep(1,times=nb-1))
    # }
    prefixsum = 0
    hdist = 0 
    for (i in 1:nb) { 
       prefixsum=prefixsum+(h1c[i]-h2c[i])  #*bindists[i]
       hdist = hdist+abs(prefixsum) 
    } 
    
    ##print(c( hdist, sum(abs(cumsum(h1c-h2c))) ))
    return(hdist)
  }
}

 

#given a list of passengerID indices, 
#show their journey time series as bar charts
#used for motivation figure



## ROBUSTNESS METRICS FOR CLUSTER CENTRE DISTANCE METRIC FOR PAPER
## Lin meaningless metrics
{
  #centres are matrix of cluster centres from eg km$centers, assume at least 2 clusters
  #can be 2 clusterings for same dataset or clusters from a different dataset
  #what about failed clusterings - all 0s
  clusterdistance <- function(centersA, centersB)
  {
    k = dim(centersA)[1] #number of clusters
    sum = 0
    for (i in 1:k) 
    {
      mind = dist(rbind(centersA[i,],centersB[1,]))[1] #initial est
      #find min dist from center Ai to any centre in B
      for (j in 2:k) {
        minij = dist(rbind(centersA[i,],centersB[j,]))[1]
        if (minij < mind) {
          mind = minij
        } 
        #print(c(minij,mind))
      }
      sum = sum+mind
      #print(c(mind,sum))
    }
    return(sum)
  }
  
  #year is 1:5, compare centres for a given dataset index datai (1..15) for 1:reps repeats
  withinset <- function(Xcenters,datai,reps)
  {
    sum = 0
    ns = 0
    for (repi in 1:reps) {
      for (repj in 1:reps) {
        #only include valid clusterings in the sum
        if (sum(Xcenters[repi,datai,,])>0 & sum(Xcenters[repj,datai,,])>0) {
          sum = sum + clusterdistance(Xcenters[repi,datai,,],
                                      Xcenters[repj,datai,,])
          ns = ns+1
        }
      }
    }
    #average of valid comparisons
    return ( c(sum,ns))
  }
  
  #between set
  #compare clusters found by different methods
  # for datasets set eg 1:nd or 1:5
  # for repeats whichreps
  betweenset <- function(Xcenters,whichds = 1:ds, whichreps = 1:reps)
  {
    sum = 0
    ns = 0
    for (ri in whichreps) {
      for (rj in whichreps) {
        for (di in whichds) {
          for (dj in whichds) {
            if (sum(Xcenters[ri,di,,])>0 & sum(Xcenters[rj,dj,,])>0) {
              sum = sum + clusterdistance(Xcenters[ri,di,,],Xcenters[rj,dj,,])
              ns = ns+1
            }
          }
        }
      }
    }
    return ( c(sum , ns))
  }
}



## MORE ROBUSTNESS METRICS (see compareTSclustering.R and compareSnapshots.R)
## clustering proximity measure
{
  #proximity for ordinal attributes
  # cl1 and cl2 are cluster labels for the same data: assume in our case 0 to 5, already ranked
  # or convert to 1 to max, then convert to (0,1.0]
  # return the cluster proximity similarity value
  ## OR POAsimilarity <- function( cl1, cl2 ) but this assumes both cluster labels are ranked
  ## which doesn't quite work for our AED clusterings
  POAsimilarity <- function( cl1, cl2 ) {
    if (length(cl1) != length(cl2)) {
      print(sprintf("ERROR: cluster vectors must have same length %d != %d",length(cl1),length(cl2)))
      return (NULL)
    }
    minr = min(c(cl1,cl2),na.rm=TRUE)
    if (minr==0) { cl1 = cl1+1; cl2 = cl2+1 } #rank should start from 1
    maxc1 = max(cl1)
    maxc2 = max(cl2) #scale by this for 0,1 scale
    scaled1 = cl1/maxc1
    scaled2 = cl2/maxc2
    dissim = abs(scaled1 - scaled2)
    sim = 1 - dissim
    avgsim = mean(sim)
    return(avgsim)
  }
}


## RAND index for comparing 2 clusterings
## matches by position even if number of clusters and IDs differ
{
## In R 
#require(mclust)
#adjustedRandIndex(x, y) for class label vectors x and y
# 1 is close 0 is random
## can have diff number of clusters and names in each - match for match by pos
## eg adjustedRandIndex(c(1,1,1,3,3,3,3,3), c(2,2,2,9,8,4,6,7)) = 0.2432432
## adjustedRandIndex(c(1,1,1,3,3,3,3,3), c(4,4,4,2,2,2,2,2))  = 1
## adjustedRandIndex(c(1,1,1,3,3,3,3,3), c(2,2,2,5,5,4,4,4)) = 0.5555556
## adjustedRandIndex(c(1,1,1,3,3,3,3,3), c(2,2,2,2,2,4,4,4)) = 0.1384615
## adjustedRandIndex(c(1,1,1,3,3,3,3,3), c(9,8,7,9,8,4,6,7)) = -0.2108108 worse than random
## classError(cl1,cl2) is similar
## Error for a given classification relative to a known truth. 
## Location of errors in a given classification relative to a known truth.
}

## Weighted Sequential Instability (quoted by Villard19 from 
  ## Leskovec, Rajaraman and Ullman, 2014 or 2020 book Mining of Massive Datasets, CUP)
  # cl seq is a vector of clusters in sequence labels 1:NC
  # distmatrix is and NCxNC matrix with Euc distance between cluster centres
{
  WSI <- function( clseq, distmatrix ) {
    n = length(clseq)
    dd = 0 #distance measure
    for (i in 1:(n-1)) {
      cldist = distmatrix[clseq[i],clseq[i+1]]
      dd = dd + cldist
    }
    dd = dd/n #shoudl this be n-1?
    return(dd)
  }
}  
}

## END 