# Code related to the evaluation of Sleep Period Time window detection algorithm in the Newcastle PSG study

## preprocess_acc_data.R
Use R package GGIR to pre-process the accelerometer binary files, e.g. feature extraction and autocalibration.

## transferdatatocsvfiles.R
Script to load the output from GGIR and store relevant data in separate csv files.

## newcastlepsg-acc-spt-detection-eval.R
Script to match PSG and accelerometer files and perform the comparison between PSG sleep scores and accelerometer derived sleep classifications. Function calculate_hdcza.R is used to detect the SPT-window with algorithm HDCZA.

See for more information [insert reference to paper: van Hees bioRxiv et al. 2018].