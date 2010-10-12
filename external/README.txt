Downloaded from: http://genetics.bwh.harvard.edu/rare_variants/

-------

This script implements the VT test for pooled association of rare variants with a phenotype.
See Price et al. AJHG 2010
The script also includes implementations of T1, T5, and WE (Madsen-Browning) tests, optionally weighted with PolyPhen scores.
(see the paper for details).

For speed, the script supports three modes of running: local on single CPU, multicore, and cluster (see options below).

NOTE: currently, the script is configured to run the tests on a single gene.
 
Usage Rscript rareVariantTests.R -p <permutations> -n <nodes> -a <phenotypeFile> -b <snpWeightFile> -c <genotypeFile> [--multicore] [--seed seed]
   <permutations> is an integer number of permutations to perform 
   <nodes> is an integer number of nodes in the SunGridEngine to use (set 0 to run the script locally, without a cluster)
   <phenotypeFile> is the name of a file with lines of format: individualID phenotypeValue
   <snpWeightFile> is the name of a file with lines of format (weight is between 0 and 1): snpid weight
   <genotypeFile> is the name of a file with lines of format (genotype is the number of rare alleles of the snp in the individual, typically one o\
f {0,1,2}): individualID snpid genotype
   --multicore is an optional flag that indicates that mutliple CPUs of the machine can be used for computations (using this flag implies that the\
 cluster is not to be used)
   <seed> is an optional random seed value
                                                                                                                                                 
Example 1 (run on multicore):                                                                                                                     
 Rscript --slave rareVariantTests.R -p 100000 -n 0 -a AHITUV/data1.ind0 -b AHITUV/data1.csnp0 -c AHITUV/data1.cgeno0 --multicore                 
                                                                                                                                                  
Example 2 (run on a single CPU):                                                                                                                  
 Rscript --slave rareVariantTests.R -p 100000 -n 0 -a AHITUV/data1.ind0 -b AHITUV/data1.csnp0 -c AHITUV/data1.cgeno0                              

Example 3 (run on Sun Grid Engine cluster - split the permutations into 50 parts):                                                              
 Rscript --slave rareVariantTests.R -p 100000 -n 50 -a AHITUV/data1.ind0 -b AHITUV/data1.csnp0 -c AHITUV/data1.cgeno0                             

Output: p-values for
 score1, score1P - test T1 and T1P (see paper)
 score2, score2P - test T5 and T5P (see paper)
 score3, score3P - test WE and WEP (see paper)
 score4, score4P - test VT and VTP (see paper)

This R implementation by Adam Kiezun, based on reference implementation in C by Alkes Price.

----------------
Example: Rscript --slave rareVariantTests.R -p 100000 -n 0 -a AHITUV/data1.ind0 -b AHITUV/data1.csnp0 -c AHITUV/data1.cgeno0 --multicore

Welcome to Rsge
    Version: 0.6.3 

  score1  score1P score2  score2P   score3  score3P   score4  score4P
1     24 16.50323     22 15.50323 619.7794 421.4031 2.800978 3.258393
running  100000  permutations on  8  cores  
counts of how often permuted data has higher value than unpermuted data
  score1 score1P score2 score2P score3 score3P score4 score4P
1   3197    1305   5391    2396   1052     460    924     188
p-values
      score1    score1P     score2    score2P     score3     score3P
1 0.03197968 0.01305987 0.05391946 0.02396976 0.01052989 0.004609954
       score4     score4P
1 0.009249908 0.001889981


