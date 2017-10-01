rm(list=ls())

f0 = c() #file to start with if used in serial analyses
f1 = c() #file to end with if used in serial analyses (modify accordingly, if infinite then it will process until last file)
mode= c(1) #What part of the analysis needs to be done (options: 1,2,3,4 and 5)
#======================================================

datadir = "/media/windows-share/Exeter/psg_study/accdata_psg_nc"
# datadir = "/media/windows-share/London/testbinaries"
outputdir = "/media/windows-share/Exeter/psg_study"
studyname = "psg"
dirR = "/home/vincent/GGIR/mcs-acc/R"
loglocation = ""

#=====================================================================================
# library(GGIR)
# load functions directly from local clone of the R package repository
ffnames = dir(dirR) # creating list of filenames of scriptfiles to load
for (i in 1:length(ffnames)) {
  source(paste(dirR,"/",ffnames[i],sep="")) #loading scripts for reading geneactiv data
}

# library("MASS")
# library("mmap")
# library("GENEAread")
# library("bitops")
# library("matlab")
# library("signal")
# library("tuneR")
# library("zoo")
# library("data.table")

g.shell.GGIR(#=======================================
             # INPUT NEEDED:
             #-------------------------------
             # General parameters
             #-------------------------------
             mode=mode, #specify above
             datadir=datadir, #specify above
             outputdir=outputdir, #specify above
             studyname=studyname, #specify above
             f0=f0, #specify above
             f1=f1, #specify above
             overwrite = TRUE, #overwrite previous milestone data?
             do.imp=TRUE, # Do imputation? (recommended)
             idloc=1, #id location (1 = file header, 2 = filename)Rcpp::
             print.filename=TRUE,
             storefolderstructure = TRUE,
             #-------------------------------
             # Part 1 parameters:
             #-------------------------------
             # Key functions: reading file, auto-calibration, and extracting features
             windowsizes = c(5,900,3600), #Epoch length, non-wear detection resolution, non-wear detection evaluation window
             do.cal= TRUE, # Apply autocalibration? (recommended)
             do.enmo = TRUE, #Needed for physical activity analysis
             do.anglez=TRUE, #Needed for sleep detection
             do.angley=TRUE,
             do.anglex=TRUE,
             
             # do.roll_med_acc_x=TRUE,
             # do.roll_med_acc_y=TRUE,
             # do.roll_med_acc_z=TRUE,
             # do.dev_roll_med_acc_x=TRUE,
             # do.dev_roll_med_acc_y=TRUE,
             # do.dev_roll_med_acc_z=TRUE,
             
             chunksize=1, #size of data chunks to be read (value = 1 is maximum)
             printsummary=TRUE,
             #-------------------------------
             # Part 2 parameters:
             #-------------------------------
             # Key functions: Non-wear detection, imputation, and basic descriptives
             strategy = 1, #Strategy (see tutorial for explanation)
             ndayswindow=7, #only relevant when strategy = 3
             hrs.del.start = 0, # Only relevant when strategy = 2. How many HOURS need to be ignored at the START of the measurement?
             hrs.del.end = 0, # Only relevant when strategy = 2. How many HOURS need to be ignored at the END of the measurement?
             maxdur = 2, # How many DAYS of measurement do you maximumally expect?
             includedaycrit = 5, # number of minimum valid hours in a day to attempt physical activity analysis
             L5M5window = c(0,24), #window over which to calculate L5 and M5
             M5L5res = 10, #resolution in minutes of M5 and L5 calculation
             winhr = c(5,10), # size of M5 and L5 (5 hours by default)
             
             qlevels = c(c(1380/1440),c(1410/1440)), #quantiles to calculate, set value at c() if you do not want quantiles
             qwindow=c(0,24), #window over which to calculate quantiles
             ilevels = c(seq(0,400,by=50),8000), #acceleration values (metric ENMO) from which a frequency distribution needs to be derived, set value at c() if you do not want quantiles
             mvpathreshold =c(100,120), #MVPA (moderate and vigorous physical activity threshold
             bout.metric = 4,
             closedbout=FALSE,
             IVIS_windowsize_minutes = 60,
             IVIS_epochsize_seconds = 3600,
             #-----------------------------------
             # Report generation
             #-------------------------------
             # Key functions: Generating reports based on meta-data
             do.report=c(), #for what parts does and report need to be generated? (option: 2, 4 and 5)
             visualreport=FALSE,
             dofirstpage = TRUE, #first page of pdf-report with simple summary histograms
             viewingwindow=1) #viewingwindow of visual report: 1 centres at day and 2 centers at night
