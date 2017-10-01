rm(list=ls())
graphics.off()
path = "/media/windows-share/Exeter/psg_study/output_accdata_psg_nc/meta/basic/"

psglogs = "/media/windows-share/Exeter/psg_study/data psg/txt"

accnames = dir(path,full.names = TRUE)
lognames = dir(psglogs,full.names = TRUE)

# load psg data and put all of it in one dataframe with dummy wake data before and after
converttime = function(x,recdate) {
  splittime = as.numeric(unlist(strsplit(as.character(x),":")))
  hour = splittime[1]
  # min = splittime[2]
  # sec = splittime[3]
  splitday = as.numeric(unlist(strsplit(recdate,"/")))
  day = splitday[1]
  month = splitday[2]
  year = splitday[3]
  if (hour < 17) recdate = paste0(day+1,"/",month,"/",year)
  t = as.character(as.POSIXlt(paste0(recdate," ",x),format="%d/%m/%Y %H:%M:%S",tz="Europe/London"))
  return(t)
}
fixmidnight = function(x) {
  a = unlist(strsplit(x," "))
  if (length(a) == 1) x = paste0(a," 00:00:00")
  return(x)
}


repropsg = TRUE
if (repropsg == TRUE) {
  for (i in 1:length(lognames)) { #
    logheader = as.matrix(read.csv(lognames[i],nrow=15,sep = "\t"))
    logdata = read.csv(lognames[i],skip=16,sep = "\t")
    recdate = logheader["Recording Date:",]
    tmp = unlist(sapply(logdata$Time..hh.mm.ss.,converttime,recdate))
    tmp2 = unlist(sapply(tmp,fixmidnight))
    logdata$time = as.POSIXlt(tmp2,tz="Europe/London")
    # ignore files with psg ending before 6am or duration < 6 hours
    if (logdata$time[length(logdata$time)]$hour >= 3 & (nrow(logdata) / 120) > 4) {
      id = as.numeric(unlist(strsplit(unlist(strsplit(lognames[i],"ecsleep"))[2],"_psg"))[1])
      logdata$stagescore = 0
      logdata$stagescore[which(logdata$Sleep.Stage == "W")] = 0
      logdata$stagescore[which(logdata$Sleep.Stage == "R")] = 1
      logdata$stagescore[which(logdata$Sleep.Stage == "N1")] = 2
      logdata$stagescore[which(logdata$Sleep.Stage == "N2")] = 3
      logdata$stagescore[which(logdata$Sleep.Stage == "N3")] = 4
      # Now save as csv file with id number in the file name
      write.csv(logdata,file=paste0("/media/windows-share/Exeter/psg_study/cleaned_psg/psg_participant",id,".csv"),row.names = FALSE)
    }
  }
}


# load acc data: id, date, angle and enmo
# look for most plausible match with psgscores

for (j in 1:length(accnames)) {
  print(j)
  load(accnames[j])
  M$metalong$timestampPOSIX = as.POSIXlt( M$metalong$timestamp,format="%Y-%m-%dT%H:%M:%S%z",tz="Europe/London")
  M$metashort$timestampPOSIX = as.POSIXlt( M$metashort$timestamp,format="%Y-%m-%dT%H:%M:%S%z",tz="Europe/London")
  id = as.numeric(unlist(strsplit(unlist(strsplit(accnames[j],"SLEEP"))[2],"_"))[1])
  nonwear = which(M$metalong$nonwearscore >= 2)
  if (length(nonwear) > 4) {
    # print(paste0(j," ",id))
    if (id == 23) { #data removed based on visual inspection of non-wear period
      M$metashort = M$metashort[1:21000,] 
    } else if (id == 28) {
      M$metashort = M$metashort[1:14000,]
    } else if (id == 29) {
      M$metashort = M$metashort[1:16000,]
    } else if (id == 34) {
      M$metashort = M$metashort[8000:37000,]
    } else if (id == 35 | id  == 45 | id == 53) {
      M$metashort = M$metashort[1:16000,]
    } else if (id == 38) {
      M$metashort = M$metashort[1:11000,]
    } else if (id  == 39) {
      M$metashort = M$metashort[1:13000,]
    } else if (id == 43 | id == 46) { # minor non-wear
    } else if (id == 47) {
      M$metashort = M$metashort[1:11000,]
    } else if (id == 49) {
      M$metashort = M$metashort[1:8600,] 
    } else if (id == 57) {
      M$metashort = M$metashort[1:15000,]
    } else if (id == 60) {
      M$metashort = M$metashort[1:13500,]
    }
  }
  write.csv(M$metashort,file=
              paste0("/media/windows-share/Exeter/psg_study/cleaned_acc/acc_participant",id,
                     "_location",as.character(I$header['Device_Location_Code',]),".csv"),row.names = FALSE)
}