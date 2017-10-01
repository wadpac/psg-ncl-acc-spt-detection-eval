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
output$row = output$perc_sleep_intimeinbed = output$error_waketime_min = output$error_onset = output$error_wakeup = 0
output$perc1 = output$perc2 = output$perc3 = output$perc4 = 0


for (h in 1:6) { # c(6,7,9,10,20,21,22,26)){ #1:nrow(psgdata)
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

  show = blocksize:(nrow(psgacc)-blocksize+1)
  # x11()
  jpeg(file=paste0("/media/windows-share/Exeter/examplepsg_",h,".jpeg"),width = 7,height = 7,res=400,un="in")
  par(mfrow=c(3,1),mar=c(2,4,4,3))
  plot(psgacc$time[show],psgacc$stagescore[show],pch=20,ylim=c(-0.5,4.5),type="l",main="sleep stages",bty="l",axes=FALSE,xlab="",ylab="")
  axis(side = 2,at = 0:4,labels = c("W","R","N1","N2","N3"),tick = TRUE)
  # axis(side = 2,at = 0:4,labels = c("W","R","N1","N2","N3"),tick = TRUE)
  # lines(psgacc$time[c(timeinbed$lightsout,timeinbed$lightson)],rep(0.5,2),col="blue",lwd=3)
  plot(psgacc$time[show],psgacc$ENMO[show],ylim=c(0,0.2),type="l",main="magnitude of acceleration (gravitational acceleration)",bty="l",axes=FALSE,ylab="")
  axis(side = 2,at = c(0,0.1),labels = c("0","0.1"),tick = TRUE)
  # lines(psgacc$time[c(timeinbed$lightsout,timeinbed$lightson)],rep(0.01,2),col="blue",lwd=3)
  plot(psgacc$time[show],psgacc$anglez[show],pch=20,ylim=c(-90,90),type="l", main = "angle-z (degrees)",bty="l",ylab="")
  # lines(psgacc$time[c(timeinbed$lightsout,timeinbed$lightson)],rep(50,2),col="blue",lwd=3)
  dev.off()
  # Percentage of sleep stage time that was correctly classified as time in bed
  psgacc$timeinbed = 0
  psgacc$timeinbed[timeinbed$lightsout:timeinbed$lightson] = 1
  psgacc = psgacc[(blocksize+1):(nrow(psgacc)-blocksize),] # remove blocks of dummy data
  
  
  # plot for visual check
  totalsleepstage = length(which(psgacc$stagescore > 0))
  # percentage of sleep periods that falls within the time in bed
  perc_sleep_intimeinbed = round(((totalsleepstage - length(which(psgacc$timeinbed == 0 & psgacc$stagescore > 0))) / totalsleepstage) * 100,digits=2)
  
  perc_stage = function(x,stage) {
    Tstage = length(which(x$stagescore == stage))
    perc_stage = round(((Tstage - length(which(x$timeinbed == 0 & x$stagescore == stage))) / Tstage) * 100,digits=2)
    return(perc_stage)
  }
  perc1 = perc_stage(psgacc,stage=1)
  perc2 = perc_stage(psgacc,stage=2)
  perc3 = perc_stage(psgacc,stage=3)
  perc4 = perc_stage(psgacc,stage=4)
  output$perc1[h] = perc1
  output$perc2[h] = perc2
  output$perc3[h] = perc3
  output$perc4[h] = perc4
  
  
  output$perc_sleep_intimeinbed[h] = perc_sleep_intimeinbed
  
  
  # Time (minutes) in sleep during the timewindow classified as not in bed
  allsleep = which(psgacc$stagescore >0 & psgacc$timeinbed == 1)
  onset = allsleep[1]
  wakeup = allsleep[length(allsleep)]
  sleepperiod = onset:wakeup
  data = data.frame(algo = psgacc$timeinbed[-sleepperiod],stagescore = psgacc$stagescore[-sleepperiod])
  error_waketime_min = round(length(which((data$algo == 0 & data$stagescore > 0))) / 12,digits=2)
  output$error_waketime_min[h] = error_waketime_min
  
  # onset and wakingup 'error'
  timewindowacc = range(which(psgacc$timeinbed == 1))
  onsetacc = timewindowacc[1]
  wakeupacc = timewindowacc[2]
  timewindowpsg = range(which(psgacc$stagescore > 0))
  onsetpsg = timewindowpsg[1]
  wakeuppsg = timewindowpsg[2]
  error_onset = (onsetacc - onsetpsg) / 12
  error_wakeup = (wakeupacc - wakeuppsg) / 12
  output$error_onset[h] = error_onset
  output$error_wakeup[h] =  error_wakeup
    output$id[h] = id
  output$row[h] = h
  # print(output[h,])
  # print performance
  # print(paste0("row: ",h," id: ",id," ",perc_sleep_intimeinbed," ",error_waketime_min))
}

print("--------------------------")
print(output[order(output$perc_sleep_intimeinbed,decreasing = TRUE),])
# print(output)