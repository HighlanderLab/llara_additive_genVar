# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data
# pValueGxY is created outside the script

# EYT2 and select variety
EYT2 = setPheno(EYT1, reps = repEYT, fixEff = year, p = pValueGxY[year])

# EYT1
EYT1 = selectInd(AYT, nEYT, use = "pheno")
EYT1 = setPheno(EYT1, reps = repEYT, fixEff = year, p = pValueGxY[year])

# AYT
AYT = selectInd(PYT, nAYT, use = "pheno")
AYT = setPheno(AYT, reps = repAYT, fixEff = year, p = pValueGxY[year])

# PYT
if(ncol(ebv(HDRW)) == 0){
  PYT = selectWithinFam(HDRW, famMax, use = "pheno")
  PYT = selectInd(PYT, nPYT, use = "pheno")
}else{
  PYT = HDRW
}
PYT = setPheno(PYT, reps = repPYT, fixEff = year, p = pValueGxY[year])

# HDRW
DH = setEBV(DH, gsModel)
output$accDH[year] = cor(ebv(DH), gv(DH))
DH = selectWithinFam(DH, famMax, use = "ebv")
HDRW = selectInd(DH, nPYT, use = "ebv")

# DH
DH = makeDH(F1, nDH)


