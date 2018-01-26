rm(list=ls())
graphics.off()
# path = "/media/vincent/Exeter/psg_study/output_accdata_psg_nc/meta/basic/"
path = "/media/vincent/Exeter/psg_study/output_psgaxivityfiles/meta/basic/"

# psglogs = "/media/vincent/Exeter/psg_study/data psg/txt"
psglogs = "/media/vincent/Exeter/Geneactiv validation data/Sleep Staging data"


accnames = dir(path,full.names = TRUE)
lognames = dir(psglogs,full.names = TRUE)

# load psg data and put all of it in one dataframe with dummy wake data before and after
converttime = function(x,recdate) {
  AMPM = as.character(unlist(strsplit(as.character(x)," "))[2])
  if (is.na(AMPM) == FALSE) {
    x = unlist(strsplit(as.character(x)," "))[1]
    splittime = as.numeric(unlist(strsplit(as.character(x),":")))
    hour = splittime[1]
    minute = splittime[2]
    second = splittime[3]
    if (AMPM == "PM") hour = hour + 12
    splitday = as.numeric(unlist(strsplit(recdate,"/")))
    month = splitday[1]
    day = splitday[2]
    year = splitday[3]
    if (hour < 17) recdate = paste0(month,"/",day+1,"/",year)
    
    t = as.character(as.POSIXlt(paste0(recdate," ",hour,":",minute,":",second),format="%m/%d/%Y %H:%M:%S",tz="America/New_York"))
    if (is.na(t) == TRUE) {
      if (day >= 28 & hour < 17) {
        recdate = paste0(month+1,"/",1,"/",year)
        t = as.character(as.POSIXlt(paste0(recdate," ",hour,":",minute,":",second),format="%m/%d/%Y %H:%M:%S",tz="America/New_York"))
      }
    }
  } else {
    t = NA
  }
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
    logheader = as.matrix(read.csv(lognames[i],nrow=11,sep = "\t")[,1:2])
    logdata = read.csv(lognames[i],skip=13,sep = "\t",header = TRUE) #read.csv(lognames[i],skip=12,sep = "\t")
    logheader = as.data.frame(logheader[,2],row.names = logheader[,1])
    recdate = unlist(strsplit(as.character(logheader["Study Date:",])," "))[1]
    tmp = unlist(sapply(logdata$Start.Time,converttime,recdate))
    NAvalues = which(is.na(tmp) == TRUE)
    if (length(NAvalues) > 0) {
      tmp = tmp[-NAvalues]
      logdata = logdata[-NAvalues,]
    }
    tmp2 = unlist(sapply(tmp,fixmidnight))
    logdata$time = as.POSIXlt(tmp2,tz="America/New_York")
    # ignore files with psg ending before 6am or duration < 6 hours
    if (logdata$time[length(logdata$time)]$hour >= 3 & (nrow(logdata) / 120) > 4) {
      id = as.numeric(unlist(strsplit(unlist(strsplit(lognames[i],"WIN"))[2],"_sta"))[1])
      logdata$stagescore = 0
      logdata$stagescore[which(logdata$Sleep.Stage == "W")] = 0
      logdata$stagescore[which(logdata$Sleep.Stage == "R")] = 1
      logdata$stagescore[which(logdata$Sleep.Stage == "N1")] = 2
      logdata$stagescore[which(logdata$Sleep.Stage == "N2")] = 3
      logdata$stagescore[which(logdata$Sleep.Stage == "N3")] = 4
      # Now save as csv file with id number in the file name
      write.csv(logdata,file=paste0("/media/vincent/Exeter/psg_study/cleaned_psg_penn/psg_participant",id,".csv"),row.names = FALSE)
    }
  }
}
print("psg done")

# load acc data: id, date, angle and enmo
# look for most plausible match with psgscores

for (j in 1:length(accnames)) {
  print(j)
  load(accnames[j])
  M$metalong$timestampPOSIX = as.POSIXlt( M$metalong$timestamp,format="%Y-%m-%dT%H:%M:%S%z",tz="Europe/London") # conversion to different timezone not needed (a bit strange)
  M$metashort$timestampPOSIX = as.POSIXlt( M$metashort$timestamp,format="%Y-%m-%dT%H:%M:%S%z",tz="Europe/London")
  id = as.numeric(unlist(strsplit(unlist(strsplit(accnames[j],"000000"))[2],"[.]cwa"))[1])
  M$metalong = M$metalong[1:72,]
  M$metashort = M$metashort[1:12960,]
  nonwear = which(M$metalong$nonwearscore >= 2)
  if (length(nonwear) > 4) {
    #& id %in% c(1009,1017,1027,1029,1005,1007,1022,1025,1003,1012,1013,1015,1018,1019,1026,1028) == FALSE
    print(paste0(j," ",id))
    if (id == 1029 | id == 1013 | id == 1025 | id == 1012 | id == 1028) { #data removed based on visual inspection of non-wear period
      M$metashort = M$metashort[1:10000,]
    } else if (id == 1022) {
      M$metashort = M$metashort[1:3000,]
    } else if (id == 1003) {
      M$metashort = M$metashort[1:9000,]
    } else if (id == 1015) {
      M$metashort = M$metashort[1000:12000,]
    } else if (id  == 1018) {
      M$metashort = M$metashort[1:12000,]
    } else if (id == 1019) { 
      M$metashort = M$metashort[1800:11000,]
    } else if (id == 1026) {
      M$metashort = M$metashort[1:11000,]
    }
    # x11()
    # par(mfrow=c(2,1))
    # plot(M$metashort$anglez,main="z")
    # plot(M$metalong$nonwearscore)
  }
  write.csv(M$metashort,file=
              paste0("/media/vincent/Exeter/psg_study/cleaned_acc_penn/acc_participant",id,
                     "_location",as.character(I$header['Device_Location_Code',]),".csv"),row.names = FALSE)
}
