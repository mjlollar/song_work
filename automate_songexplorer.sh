#1/bin/bash

# Generate directories, run downsample, run songxplorer on RIL line

# Start in directory containing subset of wav files.
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

#### Broken from here down, does not recognize songexplorer call right now. Try to fix later
#cd matlab
##Run Classify and Ethogram on downsampled wavs
#wavfiles=(*)
#for wavfile in ${wavfiles[@]} ; do
#	$SONGEXPLORER_BIN classify.py \
#	--context_ms=204.8 \
#	--shiftby_ms=0.0 \
#	--model=/home/matt/songexplorer_work/train_1r-20221014T200134Z-001/train_1r/frozen-graph.ckpt-16000.pb \
#	--model_labels=/home/matt/songexplorer_work/train_1r-20221014T200134Z-001/train_1r/labels.txt \
#	--wav=$PWD/${wavfile} \
#	--parallelize=8192 \
#	--audio_tic_rate=5000 \
#	--audio_nchannels=1 \
#	--video_findfile=same-basename \
#	--video_bkg_frames=0 \
#	--video_frame_rate=0 \
#	--video_frame_height=0 \
#	--video_frame_width=0 \
#	--video_channels=0 \
#	--deterministic=0 \
#	--labels= \
#	--prevalences=0.2,0.2,0.2,0.2,0.2
#	$SONGEXPLORER_BIN ethogram.py /home/matt/songexplorer_work/train_1r-20221014T200134Z-001 train_1r thresholds.ckpt-16000.csv $PWD/${ethofile} 5000
#done

#mv *.log ../logs
#mv *-*-*.wav ../temp
