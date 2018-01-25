#=============================
calculate_hdcza = function(angle, k =60, perc = 0.1, inbedthreshold = 15, bedblocksize = 30, outofbedsize = 60, ws3 = 5) {
  # exploratory function 27/7/2017
  medabsdi = function(angle) {
    angvar = stats::median(abs(diff(angle))) #50th percentile, do not use mean because that will be outlier dependent
    return(angvar)
  }
  x = zoo::rollapply(angle, k, medabsdi) # 5 minute rolling median of the absolute difference
  nomov = rep(0,length(x)) # no movement
  inbedtime = rep(NA,length(x))
  pp = quantile(x,probs=c(perc)) * inbedthreshold 
  # if (pp == 0) pp = 7
  if (pp < 0.13) pp = 0.13 # needed because dummy data is inserted
  if (pp > 0.50) pp = 0.50 # needed because dummy data is inserted
  nomov[which(x < pp)] = 1
  nomov = c(0,nomov,0)
  s1 = which(diff(nomov) == 1) #start of blocks in bed
  e1 = which(diff(nomov) == -1) #end of blocks in bed
  bedblock = which((e1 - s1) > ((60/ws3)*bedblocksize*1)) #which are the blocks longer than bedblocksize in minutes?
  if (length(bedblock) > 0) { #
    s2 = s1[bedblock] # only keep the bedblocks that are long enough
    e2 = e1[bedblock] # only keep the bedblocks that are long enough
    for (j in 1:length(s2)){
      inbedtime[ s2[j]:e2[j]] = 1 #record these blocks in the inbedtime vector
    }
    # fill up gaps in time between bed blocks
    outofbed = rep(0,length(inbedtime))
    outofbed[which(is.na(inbedtime) == TRUE)] = 1
    outofbed = c(0,outofbed,0)
    s3 = which(diff(outofbed) == 1) #start of blocks out of bed?
    e3 = which(diff(outofbed) == -1) #end blocks out of bed?
    outofbedblock = which((e3 - s3) < ((60/ws3)*outofbedsize*1))
    if (length(outofbedblock) > 0) { # only fill up gap if there are gaps
      s4 = s3[outofbedblock]
      e4 = e3[outofbedblock]
      if (length(s4) > 0) {
        for (j in 1:length(s4)){
          inbedtime[ s4[j]:e4[j]] = 1
        }
      }
    }
    if (length(inbedtime) == (length(x)+1)) inbedtime = inbedtime[1:(length(inbedtime)-1)]
    # keep indices for longest in bed block:
    inbedtime2 = rep(1,length(inbedtime))
    inbedtime2[which(is.na(inbedtime) == TRUE)] = 0
    s5 = which(diff(c(0,inbedtime2,0)) == 1) #start of blocks out of bed?
    e5 = which(diff(c(0,inbedtime2,0)) == -1) #end blocks out of bed?
    inbeddurations = e5 - s5
    longestinbed = which(inbeddurations == max(inbeddurations))
    lightsout = s5[longestinbed] - 1
    lightson = e5[longestinbed] - 1
  } else {
    lightson = c()
    lightsout = c()
    tib.threshold = c()
  }
  #
  # Sleep detection:
  angle[which(is.na(angle) == T)] = 0
  cnt = 1
  sdl1 = rep(0,length(angle)) # time => angle
  ws3 = 5
  timethreshold = 5
  anglethreshold = 5
  postch = which(abs(diff(angle)) > anglethreshold) #posture change of at least j degrees
  # count posture changes that happen less than once per ten minutes
  q1 = c()
  if (length(postch) > 1) {
    
    q1 = which(diff(postch) > (timethreshold*(60/ws3))) #less than once per i minutes
  }
  if (length(q1) > 0) {
    for (gi in 1:length(q1)) {
      sdl1[postch[q1[gi]]:postch[q1[gi]+1]] = 1 #periods with no posture change
    }
  } else { #possibly a day without wearing
    if (length(postch) < 10) {  #possibly a day without wearing
      sdl1[1:length(sdl1)] = 1 #periods with no posture change
    } else {  #possibly a day with constantly posture changes
      sdl1[1:length(sdl1)] = 0 #periodsposture change
    }
  }
  sleep = sdl1
  sleep = as.data.frame(sleep)
    
  
  
  #
  tib.threshold = pp
  invisible(list(lightsout=lightsout,lightson=lightson,tib.threshold=tib.threshold, sleep=sleep))
}
