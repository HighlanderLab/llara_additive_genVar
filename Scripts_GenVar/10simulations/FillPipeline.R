#Set initial yield trials with unique individuals
for(year in 1:7){
  # Crossing
  F1 = randCross(Parents,nCrosses)
  if(year<7){
    # DH
    DH = makeDH(F1,nDH)
  }
  if(year<6){
    # HDRW
    HDRW = setPheno(DH,reps=repHDRW)
  }
  if(year<5){
    # PYT
    PYT = selectWithinFam(HDRW,famMax)
    PYT = selectInd(PYT,nPYT)
    PYT = setPheno(PYT,reps=repPYT)
  }
  if(year<4){
    # AYT
    AYT = selectInd(PYT,nAYT)
    AYT = setPheno(AYT,reps=repAYT)
  }
  if(year<3){
    # EYT1
    EYT1 = selectInd(AYT,nEYT)
    EYT1 = setPheno(EYT1,reps=repEYT)
  }
  if(year<2){
    # EYT2 and select variety
    EYT2 = setPheno(EYT1,reps=repEYT)
  }
}
