rm(list=ls())
graphics.off()

pathpsg = "/media/windows-share/Exeter/psg_study/cleaned_psg"
pathacc = "/media/windows-share/Exeter/psg_study/cleaned_acc"
source("~/GGIR/psg-newcastle/tib.R")

namespsg = dir(pathpsg, full.names = TRUE)
namesacc = dir(pathacc, full.names = TRUE)
NN = length(namespsg)+length(namesacc)
comp = data.frame(fname=rep(" ",NN),method=rep(" ",NN),time=rep(" ",NN),id=rep(0,NN),stringsAsFactors = FALSE)
cnt = 1
for (i in 1:length(namespsg)) {
  A = read.csv(namespsg[i],nrow=2)
  comp$fname[cnt] = namespsg[i]
  comp$method[cnt] = "psg"
  comp$time[cnt] = as.character(A$time[1])
  comp$id[cnt] = unlist(strsplit(unlist(strsplit(namespsg[i],"pant"))[2],".cs"))[1]
  cnt = cnt + 1
}
for (j in 1:length(namesacc)) {
  B = read.csv(namesacc[j],nrow=2)
  comp$fname[cnt] = namesacc[j]
  comp$method[cnt] = "acc"
  comp$time[cnt] = as.character(B$timestampPOSIX[1])
  comp$id[cnt] = unlist(strsplit(unlist(strsplit(namesacc[j],"pant"))[2],"_loca"))[1]
  cnt = cnt + 1
}

# find the pairs based on id and date....
compacc = comp[comp$method=="acc",]
comppsg = comp[comp$method=="psg",]
compaccR = compacc[grep("right wrist",compacc$fname),]
compaccL = compacc[grep("left wrist",compacc$fname),]

psgR = merge(comppsg,compaccR,by="id")
psgR$time.psg = as.numeric(as.POSIXlt(psgR$time.x))
psgR$time.accR = as.numeric(as.POSIXlt(psgR$time.y))
diffR = (psgR$time.accR - psgR$time.psg) / 3600

psgL = merge(comppsg,compaccL,by="id")
psgL$time.psg = as.numeric(as.POSIXlt(psgL$time.x))
psgL$time.accL = as.numeric(as.POSIXlt(psgL$time.y))
diffL = (psgL$time.accL - psgL$time.psg) / 3600

# verify by looking at signal...
# psgdata = psgR
psgdata = psgL
output = data.frame(id=rep(0,nrow(psgdata)))
output$Sens1 = 0 #output$tib.threshold =  0
#output$row = output$error_onset = output$error_wakeup = 0 
# output$perc1 = output$perc2 = output$perc3 = output$perc4 = 0

pdf(file=paste0("/media/windows-share/Exeter/psg_left.pdf"),width = 7,height = 3.5)
for (h in 1:nrow(psgdata)) { # c(6,7,9,10,20,21,22,26)){ #1:nrow(psgdata)
  # print(h)
  id = as.numeric(unlist(strsplit(unlist(strsplit(psgdata$fname.x[h],"pant"))[2],".cs"))[1])
  PSG = read.csv(psgdata$fname.x[h])
  ACC = read.csv(psgdata$fname.y[h])
  PSG$time = as.POSIXlt(PSG$time,origin="1970-1-1",tz = "Europe/London")
  ACC$timestampPOSIX = as.POSIXlt(ACC$timestampPOSIX,origin="1970-1-1",tz = "Europe/London")
  # repeat 6 times to match resolution of ACC
  PSG = PSG[rep(seq_len(nrow(PSG)), each=6),]
  starttime = as.numeric(PSG$time[1])
  PSG$time = as.POSIXlt(seq(starttime,starttime+(nrow(PSG)*5)-1,by=5),origin="1970-1-1",tz = "Europe/London")
  PSG$timenum = as.numeric(PSG$time)
  ACC$timenum = as.numeric(ACC$timestampPOSIX)
  # merge psg and acc data into one dataframe based on timestamps
  psgacc = data.frame()
  cnt = 0
  cnt2 = 0
  while(nrow(psgacc) == 0) {
    PSG$timenum = PSG$timenum + 1 # psg times are not rounded to 5 seconds
    psgacc = merge(PSG,ACC,by="timenum")
    if (cnt == 10) {
      if (PSG$timenum[1] - ACC$timenum[1] > 12*3600) {
        print(paste0("timeshift backward",id))
        PSG$timenum = PSG$timenum - 24*3600
      } else if (PSG$timenum[1] - ACC$timenum[1] < -12*3600) {
        print(paste0("timeshifted forward ",id))
        PSG$timenum = PSG$timenum + 24*3600 
      }
      cnt = 0
      cnt2 = cnt2 + 1
      PSG$time = as.POSIXlt(PSG$timenum,origin="1970-1-1",tz = "Europe/London")
      if (cnt2 == 3) psgacc = 1
    }
    cnt = cnt + 1
  }
  # when pairs are clear, look for optimal shift
  radi = 120
  shift = rep(0,(2*radi)+1)
  radilist = -radi:radi
  for (g in radilist) {
    t2 = nrow(psgacc)
    shift[g+radi+1] = length(which(psgacc$stagescore[(radi+g+1):(t2-radi+g)] == 0 & psgacc$ENMO[(radi+1):(t2-radi)] > 0.02))
  }
  finalshift = radilist[median(which.min(shift))]
  # print(paste0(finalshift," ",min(shift)))
  
  # shift the timestamps for psg
  newtime = psgacc$time[(radi+finalshift+1):(t2-radi+finalshift)]
  psgacc = psgacc[(radi+1):(t2-radi),]
  psgacc$time = newtime
  psgacc = merge(PSG,ACC,by="timenum")

  # apply timeinbed detection and compare
    
  # first expand data with dummy data to create realistic conditions
  blocksize = 12 * 60 * 1.5
  angleblock = rnorm(n = blocksize,mean = 0,sd = 10) + sin((1:blocksize)/pi*0.1) * 40
  psgacc_expand = psgacc[1:blocksize,]
  psgacc_expand$anglez = angleblock
  psgacc_expand$timenum = seq(psgacc$timenum[1]-(5+blocksize*5)+1,psgacc$timenum[1]-5,by=5)
  psgacc_expand$time = as.POSIXlt(psgacc_expand$timenum,origin="1970-1-1",tz = "Europe/London")
  psgacc_expand$stagescore = 0
  psgacc_expand$ENMO = 0.1

  psgacc = rbind(psgacc_expand,psgacc)
  psgacc_expand$timenum = seq(psgacc$timenum[nrow(psgacc)]+5,psgacc$timenum[nrow(psgacc)]+(5*blocksize),by=5)
  psgacc_expand$time = as.POSIXlt(psgacc_expand$timenum,origin="1970-1-1",tz = "Europe/London")
  psgacc_expand$stagescore = 0
  psgacc_expand$ENMO = 0.1
  psgacc = rbind(psgacc,psgacc_expand)

  timeinbed = inbed(psgacc$anglez, k =60, perc = 0.1, inbedthreshold = 15, bedblocksize = 30, outofbedsize = 60, ws3 = 5)
  
  #===================================================================
  # plot for visual check
  show = blocksize:(nrow(psgacc)-blocksize+1)
  # x11()
  par(mfrow=c(2,1),mar=c(2,5,0,2))
  CL =0.9
  CA = 0.8
  plot(psgacc$time[show],psgacc$stagescore[show],pch=20,ylim=c(-0.5,4.5),type="l",main="",bty="l",axes=FALSE,xlab="",ylab="psg sleep stage", cex.lab = CL)
  axis(side = 2,at = 0:4,labels = c("Wake","REM","N1","N2","N3"),tick = TRUE,cex.axis=0.6,las=1)
  text(x = psgacc$time[show][2],y=1,labels = paste0(id),cex = 2,pos = 3,col="green")
  timeinbedt0t1 = c(timeinbed$lightsout,timeinbed$lightson)
  plot(psgacc$time[show],psgacc$anglez[show],pch=20,ylim=c(-90,90),type="l", main = "",bty="l",ylab="angle (degrees)",xlab="Time",cex.lab=CL,cex.axis=CA)
  if (length(timeinbedt0t1) == 2) lines(psgacc$time[timeinbedt0t1],rep(89,2),col="gray",lwd=3,lty=1,lend=2)
  
  psgacc$timeinbed = 0
  if (length(timeinbedt0t1) == 2) psgacc$timeinbed[timeinbed$lightsout:timeinbed$lightson] = 1
  psgacc = psgacc[(blocksize+1):(nrow(psgacc)-blocksize),] # remove blocks of dummy data
  
  #=================================================
  # Define sleep and wake
  definesleep = c(1,2,3,4) # option to redefine which PSG labels are counted towards sleep
  definewake = c(0) # option to redefine which PSG labels are counted towards wake
  
  #=====================================================
  # Sens 1: Percentage of sleep stage time that was correctly classified as time in bed
  totalsleepstage = length(which(psgacc$stagescore %in% definesleep == TRUE))
  perc_psg_eq_accinbed = round(((totalsleepstage -
                                   length(which(psgacc$timeinbed == 0 &
                                                                  psgacc$stagescore %in% definesleep ==TRUE)))
                                / totalsleepstage) * 100,digits=2)
  output$Sens1[h] = perc_psg_eq_accinbed
  # Negative predictive value
  totaloutofbed = length(which(psgacc$timeinbed == 0))
  
  if (totaloutofbed > (12*60)) {
    output$NPV[h] = round(((totaloutofbed -
                               length(which(psgacc$timeinbed == 0 &
                                              psgacc$stagescore %in% definesleep ==TRUE)))
                            / totaloutofbed) * 100,digits=2)
  } else {
    output$NPV[h] = NA
  }
  output$totaloutofbed[h] = totaloutofbed/12 # in minutes
  # #===================================
  # # Prepare for calculation of NPV and Sens2
  # tmp0 = psgacc$stagescore
  # notwakeorsleepi = which(tmp0 %in% definewake == FALSE & tmp0 %in% definesleep == FALSE)
  # if (length(notwakeorsleepi) > 0) tmp0[notwakeorsleepi] = 0.5
  # wakei = which(tmp0 %in% definewake ==TRUE)
  # if (length(wakei) > 0) tmp0[wakei] = 0
  # sleepi = which(tmp0 %in% definesleep) 
  # if (length(sleepi) > 0) tmp0[sleepi] = 1
  # tmp1 = cumsum(c(0,tmp0))
  # SBS = 60 # minutes of block size to downsample to
  # stagescore_downsampled = diff(tmp1[seq(1,length(tmp1),by=12*SBS)]) / (12*SBS)
  # tmp3 = cumsum(c(0,psgacc$timeinbed))
  # timeinbeddownsampled = diff(tmp3[seq(1,length(tmp3),by=12*SBS)]) / (12*SBS)
  # coverage = 0.2
  # output$Nblocks_outofbed[h] = length(which(timeinbeddownsampled < coverage))
  # min_Nblocks_outofbed= 1
  # #===================================
  # # nega_predic_value based on 30 minute data blocks
  # if (length(which(timeinbeddownsampled < coverage)) >= min_Nblocks_outofbed) {
  #   A1 = length(which(timeinbeddownsampled < coverage & stagescore_downsampled < coverage)  )
  #   B1 = length(which(timeinbeddownsampled < coverage))
  #   ratio = A1 / B1
  #   output$nega_predic_value[h] = ratio * 100
  # } else {
  #   output$nega_predic_value[h] = NA
  # }
  # #===================================
  # # Sens 2: sensitivity based on 30 minute data blocks
  # if (length(which(timeinbeddownsampled > (1-coverage))) >= min_Nblocks_outofbed) {
  #   ratio = length(which(timeinbeddownsampled > (1-coverage) & stagescore_downsampled > (1-coverage))  ) / length(which(stagescore_downsampled > (1-coverage)))
  #   output$Sens2[h] = ratio * 100
  # } else {
  #   output$Sens2[h] = NA
  # }
  output$id[h] = id
  output$tib.threshold[h] = round(timeinbed$tib.threshold,digits=2)
}
dev.off()

# prepare some variables for manual investigation
npv = output$NP[which(is.na(output$NPV) == FALSE)]
sen1 = output$Sens1[which(is.na(output$Sens1) == FALSE)]
# npv_nb = output$Nblocks_outofbed[which(is.na(output$nega_predic_value) == FALSE)]
sen1_nb = output$Nblocks_outofbed[which(is.na(output$Sens1) == FALSE)]

print("--------------------------")
print(output[order(output$Sens1,decreasing = TRUE),])
# print(output)