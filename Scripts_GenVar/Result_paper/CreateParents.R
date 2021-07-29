# Generate initial haplotypes
FOUNDERPOP = runMacs(nInd=70, 
                     nChr=21, 
                     segSites=600,
                     inbred=TRUE, 
                     species="WHEAT")

# Set simulation parameters
SP = SimParam$new(FOUNDERPOP)
SP$restrSegSites(100,500)   # maxQtl = 100 and maxSnp = 500 
SP$addSnpChip(500)          # nSnpPerChr
SP$addTraitAG(100,mean=4,var=0.1,varGxE=0.2)  # nQtlPerChr = 100
SP$setVarE(varE=0.4)
SP$nThreads = 1L

# Create parents
Parents = newPop(FOUNDERPOP)
rm(FOUNDERPOP)
