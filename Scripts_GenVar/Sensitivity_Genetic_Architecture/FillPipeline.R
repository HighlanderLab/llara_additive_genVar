# Set initial yield trials with unique individuals  
for(year in 1:7){
  # Crossing
  F1 = randCross(Parents, nCrosses)
  if(year < 7){
    # DH
    DH = makeDH(F1, nDH)
  }
  if(year < 6){
    # HDRW
    HDRW = setPheno(DH, reps = repHDRW, p = pValueGxY[year])
  }
  if(year < 5){
    # PYT
    PYT = selectWithinFam(HDRW, famMax, use = "pheno")
    PYT = selectInd(PYT, nPYT, use = "pheno")
    PYT = setPheno(PYT, reps = repPYT, p = pValueGxY[year + 1L])
  }
  if(year < 4){
    # AYT
    AYT = selectInd(PYT, nAYT, use = "pheno")
    AYT = setPheno(AYT, reps = repAYT, p = pValueGxY[year + 2L])
  }
  if(year < 3){
    # EYT1
    EYT1 = selectInd(AYT, nEYT, use = "pheno")
    EYT1 = setPheno(EYT1, reps = repEYT, p = pValueGxY[year + 3L])
  }
  if(year < 2){
    # EYT2 and select variety
    EYT2 = setPheno(EYT1, reps = repEYT, p = pValueGxY[year + 4L])
  }
}

