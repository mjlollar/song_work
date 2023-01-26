#1/bin/bash

### Place into Main directory containing RILs divided in subdirectories
### My directory paths: /home/matt/RIL_song_wavs/<RIL#>/ds/matlab/[.wav and .predicted files]

### Bash command:
#for x in */ ; do 
#	cd "$x"; /home/matt/songexplorer_work/scripts/post_songexplorer_generate_dirs.sh;  
#	cd /home/matt/songexplorer_work/RIL_song_wavs/test; 
#done

cd ds/matlab #Move into dir containing predicted and wav files
# Move temp files around
mv *-*-*.wav ../temp/
mv *.log ../logs

#Run BatchSongAnalysis
/home/matt/MATLAB/R2022b/bin/matlab -nosplash -nodisplay -sd /home/matt/Documents/MATLAB/Post_DS_Analysis/Post_DS_Analysis -r "BatchDeepSongAnalysis_32rig_bygeno_MJL_forRILs({'$PWD'}, {'mel'});quit;"

#Run BatchExport
cd ../matlab_RESULTS
/home/matt/MATLAB/R2022b/bin/matlab -nosplash -nodisplay -sd /home/matt/Documents/MATLAB -r "BatchExport_MJL('$PWD');quit;"

#Output BatchExport as '<RIL>_rawdata.csv'
mv EXPORT.csv ../../"${x}"_rawdata.csv
