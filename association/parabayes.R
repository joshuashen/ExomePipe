## a parametric Bayesian method for rare variants association
# Author: Yufeng "Joshua" Shen, Columbia University

# idea: 
# 1. null hypothesis: all variants in a gene are neutral; 
#    alternative: at least one variants in a gene increase risk
# 2. Denote the total number of variants as V, the number of risk variants as D,
#    define p(D,V) as p(data|D out of V are risk variants). then 
#   we can compute likelihood of p(d,V) for d=0,1,...,V using dynamic programming
# 3. p(data| variant i is causal) = Integral {non-central hypergeometric with OR}{d_OR}
# 4. p(data | alternative) = Integral {p(d,V)*prior(d)}{d_d}, with prior(d) 
#   determined by assumptions about PAR of the gene and PAR of the variants

## PAR calculation:  for rare variants, OR ~ RR = PAR / [Px * (1 - PAR)] + 1
##  PAR = Px / [ Px + 1/(OR - 1)] , where Px is the prevalence, OR must be larger than 1.
## if OR is smaller than 1, then transform to its inverse. 


## MCMCpack library for non-central hypergeometric density function
#  dnoncenhypergeom()
library("MCMCpack")


### Note: assume additive model
function <- dplikelihood(genotable) 
{
	gt <- read.table(genotable, sep="\t", header=T)
	attach(gt)
### format of genotable: 
###  varID, c1, c0, n1, n0, score etc
## c1: number of alleles in cases; c0: number of alleles in controls; 
## n1: number of cases; n0: number of controls
## score: prior information. could be useful later
	
	# initial condition
	numVar = nrow(gt)
	
	p = matrix(rep(1,(numVar+1)*(numVar+2)), nrow=numVar+2, ncol=numVar+1)
	# loop through variants
	p[1,] = 0 ## p(-1,) = 0
	for (v in 1:numVar) 
	{
		vindex = v + 1
		for (d in 0:v) 
		{	
			dindex = d + 2 
			p[dindex, vindex]=p[dindex - 1, pvindex -1]*palter(gt[v,]) + p[dindex,vindex-1]*pnull(gt[v,])
			
		}
	}
	
	## p[numVar+2,] is the final likelihood of each v

	### compute the posterior probability of data under alternative hypothesis (d>0)
	

}

function <- pnull(gvector)
{
	return dnoncenhypergeom(gvector[2], gvector[4], gvector[5], gvector[2]+gvector[3], 1) # OR=1 for null
}

function <- palter(gvector)
{	## question: do we want to make OR depends on MAF?  intuitively we should
	# prior probability of OR:  Gaussian  (Stephens and Balding), with variance 
	# depends on allele frequency: N(0,t),  t = 1/(p*(1-p))^0.5 / C, C is a scaling factor (C=20 seems reasonable)
	
	# make a vector containing OR and prior
	c = 10  ## scaling factor for prior variantce of effect. 
	### c=15 -> t(freq=0.1) = 0.22, which is close to the default value (0.2) of SNPTEST.
	### c=10 -> t(freq=0.1) = 0.33
  
	### We COULD take care of uncertainties in genotype calling or missing data by integral over other possible 
 	### genotypes . 

	freq = (gvector[2] + gvector[3]) / (gvector[4] + gvector[5])	##combined frequency in cases and controls
	t = 1/(freq * (1-freq))^0.5 / c

	# only calculate approximate posterior likelihood.	
	## divide into bins of OR, then 
	# orarray = exp(rnorm(100000, 0, t))
	binmax = 10 * t
	binmin = -10 * t
	bins = 500
	barray = c(0:bins) * (binmax - binmin)/bins + binmin
	density = dnorm(barray, 0, t)
	totalp = 0
	for (i in 1:length(barray))
	{
		or = exp(barray[i])
		lh = dnoncenhypergeom(gvector[2], gvector[4], gvector[5], gvector[2]+gvector[3],or)
		totalp = totalp + lh * density[i]
	}	
	totalp = totalp / sum(density)   
	return totalp

}


