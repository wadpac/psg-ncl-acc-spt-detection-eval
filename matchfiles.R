rm(list=ls())
graphics.off()
#==================================================
# user input needed:
simulate24 = TRUE # whether to simulate 24 hours of data or not
# specify data directories
pathpsg = "/media/vincent/Exeter/psg_study/cleaned_psg"
pathacc = "/media/vincent/Exeter/psg_study/cleaned_acc"
pathparticipantinfo = "/media/vincent/Exeter/psg_study/data participants/participants_diagnosis.csv"
pathfigures = "/media/vincent/Exeter"
# load the function for hdcza algorithm
source("~/GGIR/psg-newcastle/tib.R") 
#==================================================
# extract filenames
namespsg = dir(pathpsg, full.names = TRUE)
namesacc = dir(pathacc, full.names = TRUE)
NN = length(namespsg)+length(namesacc)
# match psg and accelerometer files
# initialize data frame
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

# find the pairs based on id and date:
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

# now double check matched by visually looking at signal...
pdf(file=paste0("/media/vincent/Exeter/Supplement_psg_accelerometer.pdf"),width = 7,height = 3.5)
for (location in c("left","right")) {
  if (location == "right") psgdata = psgR
  if (location == "left") psgdata = psgL
  # initialize output matrix (this is where we will store the key information needed for the paper)
  output = data.frame(id=rep(0,nrow(psgdata)))
  output$Sens1 = 0 #initialize to 0
  for (h in 1:nrow(psgdata)) { # loop through all psg recordings
    id = as.numeric(unlist(strsplit(unlist(strsplit(psgdata$fname.x[h],"pant"))[2],".cs"))[1])
    PSG = read.csv(psgdata$fname.x[h])
    ACC = read.csv(psgdata$fname.y[h])
    PSG$time = as.POSIXlt(PSG$time,origin="1970-1-1",tz = "Europe/London")
    ACC$timestampPOSIX = as.POSIXlt(ACC$timestampPOSIX,origin="1970-1-1",tz = "Europe/London")
    # repeat psg scores 6 times to match 5 second resolution of ACC vclassifications
    PSG = PSG[rep(seq_len(nrow(PSG)), each=6),]
    starttime = as.numeric(PSG$time[1])
    PSG$time = as.POSIXlt(seq(starttime,starttime+(nrow(PSG)*5)-1,by=5),origin="1970-1-1",tz = "Europe/London")
    PSG$timenum = as.numeric(PSG$time)
    ACC$timenum = as.numeric(ACC$timestampPOSIX)
    # merge psg and acc data into one dataframe based on timestamps
    psgacc = data.frame()
    cnt = cnt2 = 0
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
    # ========================================================
    # We know that time blocks are not perfectly synced (just like in all PSG labs around the world)
    # Therefore, look for optimal shift (a bit similar to http://iopscience.iop.org/article/10.1088/2057-1976/aa985f/meta)
    radi = 120
    shift = rep(0,(2*radi)+1)
    radilist = -radi:radi
    for (g in radilist) {
      t2 = nrow(psgacc)
      shift[g+radi+1] = length(which(psgacc$stagescore[(radi+g+1):(t2-radi+g)] == 0 & psgacc$ENMO[(radi+1):(t2-radi)] > 0.02))
    }
    finalshift = radilist[median(which.min(shift))]
    # shift the timestamps for psg
    newtime = psgacc$time[(radi+finalshift+1):(t2-radi+finalshift)]
    psgacc = psgacc[(radi+1):(t2-radi),]
    psgacc$time = newtime
    psgacc = merge(PSG,ACC,by="timenum")
    
    #=======================================================
    # apply spt-window detection and compare estimates
    # first expand data with simulated data to create realistic conditions
    if (simulate24 == TRUE) {
      Naddedhours = floor((24 - (nrow(psgacc) / (12*60))) / 2) # Number of hours needed to simulate 24 hours
    } else {
      Naddedhours = 1.5 # notice that further on one hour of this is removed to make evaluation specific to recording + 30 minutes
    }
    blocksize = 12 * 60 * Naddedhours
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
    # now use the expanded data for spt window detection.
    sptwindow = inbed(psgacc$anglez, k =60, perc = 0.1, inbedthreshold = 15, bedblocksize = 30, outofbedsize = 60, ws3 = 5)
    # detect sleep episodes within the spt window based on previously described algorithm: journals.plos.org/plosone/article?id=10.1371/journal.pone.0142533
    postch = which(abs(diff(psgacc$anglez)) > 5) #posture change of at least j degrees
    # count posture changes that happen less than once per ten minutes
    sdl1 = rep(0,length(psgacc$anglez))
    q1 = c()
    if (length(postch) > 1) {
      q1 = which(diff(postch) > (5*(60/5))) #less than once per i minutes
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
    sleepepisodes = sdl1
    # Derive sleep efficiency and various other summary variables & Area under the cruve for spt window classification
    if (length(sptwindow$lightson) != 0 & length(sptwindow$lightsout) != 0) {
      # Estimated sleep from the algorithm
      sleepclas = sptwindow$sleep[sptwindow$lightsout:sptwindow$lightson,1]
      output$est_sle_eff[h] = (length(which(sleepclas == 1)) / length(sleepclas)) * 100
      output$est_PSTdur[h] = length(sleepclas) / 720
      tib_est_series = rep(0,nrow(sptwindow$sleep))
      tib_est_series[sptwindow$lightsout:sptwindow$lightson] = 1
      sleepepisodes[which(tib_est_series != 1)] = 0
      sleep_est = data.frame(est = tib_est_series,time=psgacc$timenum,sleepepisodes_est = sleepepisodes)
      output$est_PSTonset[h] = strftime(as.POSIXlt(sleep_est$time[sptwindow$lightsout], origin="1970-1-1", tz="Europe/London")[1], format='%H:%M:%S')
      output$est_PSTwake[h] = strftime(as.POSIXlt(sleep_est$time[sptwindow$lightson], origin="1970-1-1", tz="Europe/London")[1], format='%H:%M:%S')
      # True sleep according to PSG
      pss = psgacc$stagescore
      psstmp = pss
      psstmp[psstmp != 0] = 1
      i234 = which(pss != 0)[1]:(length(pss)-match(1,rev(psstmp))+1)
      tib_true_series = rep(0,length(pss))
      tib_true_series[i234] = 1
      tib_true_series_time = psgacc$timenum
      sleepscore = pss[i234]
      output$true_sle_eff[h] = (length(which(sleepscore != 0)) / length(sleepscore)) * 100
      output$true_PSTdur[h] = length(sleepscore) / 720
      output$true_PSTonset[h] = as.character(psgacc$Time..hh.mm.ss.[which(pss != 0)[1]])
      output$true_PSTwake[h] = as.character(psgacc$Time..hh.mm.ss.[(length(pss)-match(1,rev(psstmp))+1)])
      
      sleepscore[which(sleepscore != 0)] = 1
      sleep_psg = data.frame(psg = tib_true_series,time=tib_true_series_time,sleepepisodes_true=psstmp)
      # merge psg and accelerometer (algorithm) based estimates
      sleep_est_psg = merge(sleep_est,sleep_psg,by="time",all = FALSE)
      if (simulate24 == FALSE) {
        sleep_est_psg = sleep_est_psg[721:(nrow(sleep_est_psg)-720),] # this will remove 1 hour of the dummy data
      }
      # Calcultate AUC
      library(pROC)
      roccurve = roc(sleep_est_psg$psg ~ sleep_est_psg$est)
      output$auc[h] = auc(roccurve)
      output$accuracy[h] = length(which(sleep_est_psg$psg == sleep_est_psg$est)) / nrow(sleep_est_psg)

    }
    output$est_sleepdur[h] = length(which(sleep_est_psg$sleepepisodes_est == 1)) / (12*60)
    output$true_sleepdur[h] = length(which(sleep_est_psg$sleepepisodes_true == 1)) / (12*60)
    #===================================================================
    # plot for visual check
    show = blocksize:(nrow(psgacc)-blocksize+1)
    par(mfrow=c(2,1),mar=c(2,5,0,2))
    CL =0.9
    CA = 0.8
    plot(psgacc$time[show],psgacc$stagescore[show],pch=20,ylim=c(-0.5,5),type="l",main="",bty="l",axes=FALSE,xlab="",ylab="psg sleep stage", cex.lab = CL)
    axis(side = 2,at = 0:4,labels = c("Wake","REM","N1","N2","N3"),tick = TRUE,cex.axis=0.6,las=1)
    text(x = psgacc$time[show][1], y=4.5, labels = paste0("id ", id, " ", location, " wrist"), cex = 0.7, pos = 4, col="red")
    sptwindowt0t1 = c(sptwindow$lightsout,sptwindow$lightson)
    plot(psgacc$time[show], psgacc$anglez[show], pch=20, ylim=c(-90,90), type="l", main = "",bty="l", ylab="angle (degrees)",xlab="Time",cex.lab=CL,cex.axis=CA)
    if (length(sptwindowt0t1) == 2) lines(psgacc$time[sptwindowt0t1], rep(89,2), col="blue", lwd=5, lty=1, lend=2)
    
    psgacc$sptwindow = 0
    if (length(sptwindowt0t1) == 2) psgacc$sptwindow[sptwindow$lightsout:sptwindow$lightson] = 1
    psgacc = psgacc[(blocksize+1):(nrow(psgacc)-blocksize),] # remove blocks of simulated data
    
    #=================================================
    # Define sleep and wake
    definesleep = c(1,2,3,4) # option to redefine which PSG labels are counted towards sleep
    definewake = c(0) # option to redefine which PSG labels are counted towards wake
    
    #=====================================================
    # Sens 1: Percentage of sleep stage time that was correctly classified as time in bed
    totalsleepstage = length(which(psgacc$stagescore %in% definesleep == TRUE))
    perc_psg_eq_accinbed = round(((totalsleepstage -
                                     length(which(psgacc$sptwindow == 0 &
                                                    psgacc$stagescore %in% definesleep ==TRUE)))
                                  / totalsleepstage) * 100,digits=2)
    output$Sens1[h] = perc_psg_eq_accinbed
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
    output$id[h] = id
    output$tib.threshold[h] = round(sptwindow$tib.threshold,digits=2)
  }
  sen1 = output$Sens1[which(is.na(output$Sens1) == FALSE)]
  sen1_nb = output$Nblocks_outofbed[which(is.na(output$Sens1) == FALSE)]
  
  print("--------------------------")
  print(output[order(output$Sens1,decreasing = TRUE),])
  convchartime2number = function(x) {
    gettime = function(y) {
      if (is.na(y) == FALSE) {
        tmp = as.numeric(unlist(strsplit(y,":")))
        y = tmp[1]+(tmp[2]/60)+(tmp[3]/3600)
        if (y < 12) y = y + 24
      }
      return(y)
    }
    x = as.numeric(sapply(X=x,FUN=gettime))
    return(x)
  }
  output$true_PSTonset_nm = convchartime2number(output$true_PSTonset)
  output$true_PSTwake_nm = convchartime2number(output$true_PSTwake)
  output$est_PSTonset_nm = convchartime2number(output$est_PSTonset)
  output$est_PSTwake_nm = convchartime2number(output$est_PSTwake)
  
  output$error_PSTonset = output$est_PSTonset_nm - output$true_PSTonset_nm
  output$error_PSTwake = output$est_PSTwake_nm - output$true_PSTwake_nm

  diagn = read.csv(file=pathparticipantinfo)
  if (location == "right") output_right = merge(output,diagn,by="id")
  if (location == "left") output_left = merge(output,diagn,by="id")
}
dev.off()

print(table(abs(round((output_left$est_PSTdur - output_left$true_PSTdur)-0.5)+0.5)))
print(table(abs(round((output_right$est_PSTdur - output_right$true_PSTdur)-0.5)+0.5)))
print(table(abs(round((output_left$est_sleepdur - output_left$true_sleepdur)-0.5)+0.5)))

addLoAlines = function(meth_diff){
  abline(h=mean(meth_diff),lty=1)
  abline(h=mean(meth_diff)-(1.96*sd(meth_diff)),lty=3)
  abline(h=mean(meth_diff)+(1.96*sd(meth_diff)),lty=3)
}

# Figure 3 in paper
CX = 0.9
CXM = 0.9
CXL = 0.9
CXA = 0.8
FL = 1 # font label
CXLegend = 0.7
YLIM = c(-7,4.1)
jpeg(file=paste0(pathfigures,"/Figure3.jpeg"),width = 7,height = 7,units = "in",res=500)
par(mfrow=c(2,2),mar=c(4,5,5,2))

#----------------------------
meth_diff = output_left$est_PSTdur - output_left$true_PSTdur
print(range(meth_diff))
norm = which(output_left$Disorder==0)
diag = which(output_left$Disorder==1)
plot(output_left$true_PSTdur[norm],meth_diff[norm],xlab="PSG (hours)",
     ylab="HDCZA - PSG (hours)",main="SPT-window duration (left wrist)",
     pch=19,bty="l",ylim=YLIM,xlim=range(output$true_PSTdur),cex=CX,cex.lab=CXL,cex.axis=CXA,cex.main=CXM,font.lab=FL)
lines(output_left$true_PSTdur[diag],meth_diff[diag],pch=1,type="p",cex=CX)
addLoAlines(meth_diff)

meth_diff = output_right$est_PSTdur - output_right$true_PSTdur
print(range(meth_diff))
norm = which(output_right$Disorder==0)
diag = which(output_right$Disorder==1)
plot(output_right$true_PSTdur[norm],meth_diff[norm],xlab="PSG (hours)",
     ylab="HDCZA  - PSG (hours)",main="SPT-window duration (right wrist)",
     pch=19,bty="l",ylim=YLIM,xlim=range(output$true_PSTdur),cex=CX,cex.lab=CXL,cex.axis=CXA,cex.main=CXM,font.lab=FL)
lines(output_right$true_PSTdur[diag],meth_diff[diag],pch=1,type="p",cex=CX)
addLoAlines(meth_diff)
#-------------------------------------
meth_diff = output_left$est_sleepdur - output_left$true_sleepdur
norm = which(output_left$Disorder==0)
diag = which(output_left$Disorder==1)
plot(output_left$true_sleepdur[norm],meth_diff[norm],xlab="PSG (hours)",
     ylab="HDCZA - PSG (hours)",main="Sleep duration (left wrist)",
     pch=19,bty="l",ylim=YLIM,xlim=range(output$true_sleepdur),cex=CX,cex.lab=CXL,cex.axis=CXA,cex.main=CXM,font.lab=FL)
lines(output_left$true_sleepdur[diag],meth_diff[diag],pch=1,type="p",cex=CX)
addLoAlines(meth_diff)

meth_diff = output_right$est_sleepdur - output_right$true_sleepdur
norm = which(output_right$Disorder==0)
diag = which(output_right$Disorder==1)
plot(output_right$true_sleepdur[norm],meth_diff[norm],xlab="PSG (hours)",
     ylab="HDCZA  - PSG (hours)",main="Sleep duration (right wrist)",
     pch=19,bty="l",ylim=YLIM,xlim=range(output$true_sleepdur),cex=CX,cex.lab=CXL,cex.axis=CXA,cex.main=CXM,font.lab=FL)
lines(output_right$true_sleepdur[diag],meth_diff[diag],pch=1,type="p",cex=CX)
addLoAlines(meth_diff)
dev.off()


# MAE calculation
output_left$error_PSTonset_abs = abs(output_left$error_PSTonset) # calculate absolute error in wake time per night
output_left$error_PSTwake_abs = abs(output_left$error_PSTwake) # calculate absolute error in onset time per night
MAE = round(mean(c(output_left$error_PSTwake_abs,output_left$error_PSTonset_abs)),digits=3) # calcualte mean acros individuals.
print(paste0("MAE left = ", MAE * 60))

output_right$error_PSTonset_abs = abs(output_right$error_PSTonset) # calculate absolute error in wake time per night
output_right$error_PSTwake_abs = abs(output_right$error_PSTwake) # calculate absolute error in onset time per night
MAE = round(mean(c(output_right$error_PSTwake_abs,output_right$error_PSTonset_abs)),digits=3) # calcualte mean acros individuals.
print(paste0("MAE right = ", MAE * 60))
