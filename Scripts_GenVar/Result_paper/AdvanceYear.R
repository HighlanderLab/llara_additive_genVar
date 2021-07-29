# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data
# p is created outside the script

# EYT2 and select variety
EYT2 = setPheno(EYT1,reps=repEYT,fixEff=year)

# EYT1
EYT1 = selectInd(AYT,nEYT)
EYT1 = setPheno(EYT1,reps=repEYT,fixEff=year)

# AYT
AYT = selectInd(PYT,nAYT)
AYT = setPheno(AYT,reps=repAYT,fixEff=year)

# PYT
PYT = selectWithinFam(HDRW,famMax)
PYT = selectInd(PYT,nPYT)
PYT = setPheno(PYT,reps=repPYT,fixEff=year)

# HDRW
HDRW = setPheno(DH,reps=repHDRW,fixEff=year)

# DH
DH = makeDH(F1,nDH)


