#1/bin/bash

#Start in directory containing subset of wav files
# make dirs separating downsample handling and full wavs
mkdir full
mkdir ds
mv *.wav full
cd full
python ../../../scripts/decimator_ondirs.py --i "$PWD" --r 5000 #downsample all wavs for model

mv *ds.wav ../ds #move downsample output to ds folder
cd ../ds/
#set up dirs for matlab output
mkdir matlab
mkdir logs
mkdir temp
mv *.wav matlab
