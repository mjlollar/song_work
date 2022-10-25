#!/bin/bash

# e.g. ./classify_auto.sh FR48N

### Input directory
WAV_dir=$1

### Move to wavs, list in array
cd /home/matt/songexplorer_work/$WAV_dir
WavFiles=$( ls *.wav )

### Call classify for each wav file in input directory
### Paths relative to my system
for wavfile in $WavFiles ; do
	echo $wavfile
	$SONGEXPLORER_BIN classify.py \
		--context_ms=204.8 \
		--shiftby_ms=0.0 \
		--model=/home/matt/songexplorer_work/train_1r-20221014T200134Z-001/train_1r/frozen-graph.ckpt-16000.pb \
		--model_labels=/home/matt/songexplorer_work/train_1r-20221014T200134Z-001/train_1r/labels.txt \
		--wav=/home/matt/songexplorer_work/$WAV_dir/$wavfile \
		--parallelize=8192 \
		--audio_tic_rate=5000 \
		--audio_nchannels=1 \
		--video_findfile=same-basename \
		--video_bkg_frames=0 \
		--video_frame_rate=0 \
		--video_frame_height=0 \
		--video_frame_width=0 \
		--video_channels=0 \
		--deterministic=0 \
		--labels= \
		--prevalences=
 done
