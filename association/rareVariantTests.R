#####################################################
# This script implements the VT test for pooled association of rare variants with a phenotype.
# See Price et al. AJHG 2010
# The script also includes implementations of T1, T5, and WE (Madsen-Browning) tests, optionally weighted with PolyPhen scores.
# (see the paper for details).
#
# For speed, the script supports three modes of running: local on single CPU, multicore, and cluster (see options below).
#
# NOTE: currently, the script is configured to run the tests on a single gene.
# 
# Usage Rscript rareVariantTests.R -p <permutations> -n <nodes> -a <phenotypeFile> -b <snpWeightFile> -c <genotypeFile> [--multicore] [--seed seed]
#   <permutations> is an integer number of permutations to perform
#   <nodes> is an integer number of nodes in the SunGridEngine to use (set 0 to run the script locally, without a cluster)
#   <phenotypeFile> is the name of a file with lines of format: individualID phenotypeValue
#   <snpWeightFile> is the name of a file with lines of format (weight is between 0 and 1): snpid weight
#   <genotypeFile> is the name of a file with lines of format (genotype is the number of rare alleles of the snp in the individual, typically one of {0,1,2}): individualID snpid genotype
#   --multicore is an optional flag that indicates that mutliple CPUs of the machine can be used for computations (using this flag implies that the cluster is not to be used)
#   <seed> is an optional random seed value
# 
# Example 1 (run on multicore):
#  Rscript --slave rareVariantTests.R -p 100000 -n 0 -a AHITUV/data1.ind0 -b AHITUV/data1.csnp0 -c AHITUV/data1.cgeno0 --multicore
#
# Example 2 (run on a single CPU):
#  Rscript --slave rareVariantTests.R -p 100000 -n 0 -a AHITUV/data1.ind0 -b AHITUV/data1.csnp0 -c AHITUV/data1.cgeno0
#
# Example 3 (run on Sun Grid Engine cluster - split the permutations into 50 parts):
#  Rscript --slave rareVariantTests.R -p 100000 -n 50 -a AHITUV/data1.ind0 -b AHITUV/data1.csnp0 -c AHITUV/data1.cgeno0
#
# Output: p-values for
# score1, score1P - test T1 and T1P (see paper)
# score2, score2P - test T5 and T5P (see paper)
# score3, score3P - test WE and WEP (see paper)
# score4, score4P - test VT and VTP (see paper)
#
# This R implementation by Adam Kiezun, based on reference implementation in C by Alkes Price.
#######################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(Rsge))
suppressPackageStartupMessages(library(doMC))

#Read the command-line options
opt <- getopt(matrix(c(
  'indFile', 'a', 1, 'character', "file with individuals and phenotypes",
  'snpFile', 'b', 1, 'character', "file with SNP weights",
  'genoFile', 'c', 1, 'character', "file with genotypes",
  'flipPhenotype', 'f', 0, 'logical', "should the phenotype be multiplied by -1",
  'permutations', 'p', 1, "integer", "number of permutations to perform",
  'nodes', 'n', 1, "integer", "number of nodes in the cluster (if 0 script will run locally)",
  'profile', 'r', 1, "character", "optional name of profile log ",
  'multicore', 'm', 0, "logical", "is it mutlicore (temporary)",
  'seed', 's', 1, "integer", "random seed"
), byrow=T, nrow=9));

if (is.null(opt$permutations)) { opt$permutations = 1000 }
if (is.null(opt$nodes)) { opt$nodes = 48 }
if (is.null(opt$seed)) { opt$seed = 0 }
if (is.null(opt$flipPhenotype)) {opt$flipPhenotype = FALSE}

#Move the random seed to avoid overlaps when we run this script many times
opt$seed <- 99991 * opt$seed 

ind <- read.table(opt$indFile, col.names=c("indid", "pheno"))
csnp <- read.table(opt$snpFile, col.names=c("snpid", "polyphen"))
cgeno <- read.table(opt$genoFile, col.names=c("indid", "snpid", "count"))

#Sometimes phenotypes are annotated the opposite way of what we're expecting. If yes, then flip.
ind <- ind[,c("indid", "pheno")]
if (opt$flipPhenotype){
  cat("flipping phenotypes\n")
  ind$pheno <- -1.0*ind$pheno
}
meanpheno <- mean(ind$pheno)

#For each SNP, how many times it is seen
csnp$counts <- sapply(csnp$snpid, function(x){ sum(cgeno[cgeno$snpid == x,]$count) })

#For a SNP with total count c, how many counts are lower than c? (ie what is the rank or c in the order of counts)
csnp$countg <- sapply(csnp$counts, function(x){ length(unique(csnp[csnp$counts < x,]$counts)) })

#Sample size
N <- dim(ind)[1]

m1 <- merge(cgeno, csnp, by=c("snpid"))

#adjust polyphen scores
if (length(m1[(m1$counts >= N/50) & (m1$polyphen < 1.0),]$polyphen) > 0){
 m1[(m1$counts >= N/50) & (m1$polyphen < 1.0),]$polyphen <- 0.5
}

#pre-compute metrics that are independent of permutations
m1$countSquare <- m1$count*m1$count
m1$countPolyphen <- m1$count*m1$polyphen
m1$countSquarePolyphenSquare <- m1$countPolyphen*m1$countPolyphen
f <- (1+m1$counts)/(2+2*N)
m1$weight <- 1/sqrt(f*(1.0-f))
m1$countWeight <- m1$count*m1$weight

#Create a single table by joining SNPs and genotypes by SNPid, and joining individuals by individual ID
m <- merge(m1, ind, by=c("indid"))
m <- m[m$counts < N,] #ignore common variants

#Compute sum for subsets of indices (the subsets are pre-computed)
mysum <- function(X, range, arr, whiches){
  for (i in range) { arr[i] <- sum(X[whiches[[i]]])}
  arr
}

#To improve speed, pre-compute everything that is independent of permutations
ctg <- m$countg
mCount <- m$count;
mPolyphen <- m$polyphen;
mCountWeight <- m$countWeight
indPheno <- ind$pheno

nx <- length(ctg)
fctg <- as.factor(list(ctg)[[1]])
index <- fctg
one <- 1L
group <- rep.int(one, nx) + one * (as.integer(index) - one)
len <- length(unique(group))
arr <- double(len)
oneToLen <- 1:len
Msize <- dim(m)[1]

whiches <- vector("list", len)
for (i in oneToLen) { whiches[[i]] <- which(group == i) }

count <- mysum(mCount, range=oneToLen, arr=arr, whiches=whiches)
countSquare <- mysum(m$countSquare, range=oneToLen, arr=arr, whiches=whiches)
countSquarePolyphenSquare <- mysum(m$countSquarePolyphenSquare, range=oneToLen, arr=arr, whiches=whiches)
countPolyphen <- mysum(m$countPolyphen, range=oneToLen, arr=arr, whiches=whiches)

#Indices of variants below frequency thresholds
mBelow50 <- which(m$counts < N/50)
mBelow10 <- which(m$counts < N/10)

#Pre-compute cumulative sums, for VT test
csCount <- cumsum(count)
csCountSquare <- cumsum(countSquare)
csCountPolyphenMeanpheno <- cumsum(countPolyphen * meanpheno)
csCountSquarePolyphenSquare <- cumsum(countSquarePolyphenSquare)
csCountMeanpheno <- csCount*meanpheno
sqrtCsCountSquare <- sqrt(csCountSquare)
sqrtCsCountSquarePolyphenSquare <- sqrt(csCountSquarePolyphenSquare)

#Indices of individuals from m in ind (may be duplicate)
matchIds <- match(m$indid, ind$indid)

#Compute the test scores for many tests, for 1 permutation
getScores <- function(permute){
  
  if (permute){
    pheno <- sample(ind$pheno)[matchIds]
  } else {
    pheno <- ind$pheno[matchIds]
  }
  
  phenoCount         <- pheno * mCount
  phenoCountPolyphen <- phenoCount * mPolyphen
  phenoCountWeight   <- pheno * mCountWeight

  #Scores that count only rare variants, optionally weighted
  score1 <- sum(phenoCount[mBelow50])
  score2 <- sum(phenoCount[mBelow10])
  score1P <- sum(phenoCountPolyphen[mBelow50])
  score2P <- sum(phenoCountPolyphen[mBelow10])

  #Madsen-Browning score, optionally weighted
  score3 <- sum(phenoCountWeight)
  score3P <- sum(phenoCountWeight * mPolyphen)

  #VT test, optionally weighted
  #Aggregate for each count, to find the optimal threshold for VT test
  csPhenoCount <- cumsum(mysum(phenoCount, range=oneToLen, arr=arr, whiches=whiches))
  csPhenoCountPolyphen <- cumsum(mysum(phenoCountPolyphen, range=oneToLen, arr=arr, whiches=whiches))
  
  score4 <- max((csPhenoCount-csCountMeanpheno)/sqrtCsCountSquare)
  score4P <- max((csPhenoCountPolyphen-csCountPolyphenMeanpheno)/sqrtCsCountSquarePolyphenSquare)

  c(score1=score1, score1P=score1P, score2=score2, score2P=score2P, score3=score3, score3P=score3P, score4=score4, score4P=score4P)
}

#Unpermuted data for which we're looking for pvalues
unpermutedScores <- as.data.frame(t(getScores(permute=FALSE)))
print(unpermutedScores)

#For a specific score, returns how often permuted data has a higher score than unpermuted data.
permwins <- function(scores.df, scorename) {
  unpermuted <- (unpermutedScores[c(scorename)])[[1]]
  ceiling(sum(scores.df[,c(scorename)] > unpermuted) + 0.5*sum(scores.df[,c(scorename)] == unpermuted))
}

#For all scores, returns how often permuted data has a higher score than unpermuted data.
getPermWins <- function(subrange, seed){
  set.seed(seed)
  scores <- sapply(X = subrange, simplify = T, USE.NAMES = T, FUN = function(x) { getScores(permute=TRUE) })
  scores.df <- as.data.frame(t(scores))
  pw1  <- permwins(scores.df, "score1")
  pw1P <- permwins(scores.df, "score1P")
  pw2  <- permwins(scores.df, "score2")
  pw2P <- permwins(scores.df, "score2P")
  pw3  <- permwins(scores.df, "score3")
  pw3P <- permwins(scores.df, "score3P")
  pw4  <- permwins(scores.df, "score4")
  pw4P <- permwins(scores.df, "score4P")
  c(score1=pw1, score1P=pw1P, score2=pw2, score2P=pw2P, score3=pw3, score3P=pw3P, score4=pw4, score4P=pw4P)
}

#P-values
pval <- function(permwins){ (permwins+1)/(opt$permutations+1) }

#Splits the range of all permutations to k parts
splitRanges <- function(permutations, k) {
  i <- 1:permutations
  if (k == 0)
    list(i)
  else
    structure(split(i, cut(i, k)), names = NULL)
}

#Start the profiler
if (!is.null(opt$profile)) { Rprof(opt$profile) }

if (is.null(opt$multicore)){
   #Running on the cluster
   cluster <- opt$nodes != 0
   cat("running ", opt$permutations, " permutations on ", opt$nodes, "nodes " , " cluster=", cluster, "\n")
   options(sge.user.options = paste(getOption("sge.user.options"), " -p -1", sep="")) #SGE job priorities

   subranges <- splitRanges(opt$permutations, opt$nodes)
   pw <- sge.parSapply(cluster=cluster, global.savelist=c("matchIds", "opt", "csCount", "csCountSquare", "csCountPolyphenMeanpheno", "csCountSquarePolyphenSquare", "csCountMeanpheno", "sqrtCsCountSquare", "sqrtCsCountSquarePolyphenSquare","indPheno", "Msize", "mCount", "mPolyphen", "mCountWeight", "mBelow50", "mBelow10", "mysum", "oneToLen", "whiches", "arr",  "m", "meanpheno", "N", "csnp", "ind", "cgeno", "subranges", "getScores", "permwins", "unpermutedScores"), function.savelist=c("getPermWins"), njobs=length(subranges), 1:length(subranges), function(i) { getPermWins(subranges[[i]], seed=opt$seed+i)})
   
   pw.df <- as.data.frame(t(pw))
   pw.df <- as.data.frame(t(apply(pw.df, 2, sum)))
} else {
  #Running multicore
  #XXX: detectCores is an internal function of the multicores package
  cores <- multicore:::detectCores(all.tests=TRUE)
  cat("running ", opt$permutations, " permutations on ", cores, " cores ", "\n")

  registerDoMC()
  subranges <- splitRanges(opt$permutations, cores)  
  pw <- foreach(i=1:cores, .combine='+') %dopar% getPermWins(subranges[[i]], seed=opt$seed+i)
  pw.df <- as.data.frame(t(pw))[1,]
}

if (!is.null(opt$profile)) { Rprof(NULL) }   

cat("counts of how often permuted data has higher value than unpermuted data\n")
print(pw.df)
cat("p-values\n")
print(pval(pw.df))
