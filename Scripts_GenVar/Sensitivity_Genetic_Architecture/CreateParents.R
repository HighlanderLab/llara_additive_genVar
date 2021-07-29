
# Generate initial haplotypes
FOUNDERPOP = runMacs(nInd = 70, 
                     nChr = 21, 
                     segSites = (nSNPs + nQTLs),
                     inbred = TRUE, 
                     species = "WHEAT")


# Set simulation parameters
SP = SimParam$new(FOUNDERPOP)
SP$restrSegSites(minQtlPerChr = nQTLs, minSnpPerChr = nSNPs, overlap = FALSE) 
SP$addSnpChip(nSNPs)                # nSnpPerChr
SP$addTraitAG(nQTLs, mean = 4, var = 0.1, varGxE = varGE)       # nQtlPerChr = nQTLs
SP$setVarE(varE = 0.4)
SP$nThreads = 1L

# Create parents
Parents = newPop(FOUNDERPOP)
rm(FOUNDERPOP)

