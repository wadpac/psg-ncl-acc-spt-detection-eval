# Code related to the analyses of psg and accelerometer data

## preprocess_acc_data.R
Uses GGIR to pre-process the accelerometer binary files.

## transferdatatocsvfiles.R
Script to load the output from GGIR and store relevant data in separate csv files.

## matchfiles.R
Script to match PSG and accelerometer files and perform the comparison between PSG sleep scores and accelerometer derived sleep classifications.

## tib.R
Function used by matchfiles.R to detect the SPT-window.
