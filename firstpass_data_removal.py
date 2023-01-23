import sys
import os
import numpy as np
import pandas as pd

infile = sys.argv[1] #use full path if csv is not in current directory.
outfile = sys.argv[2] #use full path for output, else outputs to working directory

#Updated col names as of: Dec 8 2022 - MJL
df = pd.read_csv(infile, usecols=['maternal_geno','phenotype','chamber','SinesPerMinNew','PulsesPerMinNew','SongPerMinNew','Sine2PulseNew','Sine2PulseMJL','Sine2PulseMJLbase',
  'SinesPerMinOld','PulsesPerMinOld','SongPerMinOld','Sine2PulseOld','ModeSineMFFT','ModePulseMFFT','SineFrameTotal','PulseFrameTotal','PulsesInTrain',
  'AllPulsesPerMin','PulseTrainsPerMin','MedianSineTrainLength', 'MedianSineAmplitudes','MedianPulseTrainLength','MedianPulseAmplitudes','ModePeak2PeakIPI',
  'ModeEnd2PeakIPI','SineTrainsPerMin','PulseTrainsPerMin','Number_slow_pulses','Number_fast_pulses','Slow2FastandSlowPulses','AmbientPerMin','OtherPerMin'])
### Add/Drop cutoffs as desired
df.drop(df[df['SongPerMinNew'] <= float(1.0)].index, inplace=True) #No song,sine,pulse less than 1 sec per minute recording
df.drop(df[df['PulsesPerMinNew'] <= float(1.0)].index, inplace=True)
df.drop(df[df['SinesPerMinNew'] <= float(1.0)].index, inplace=True)
df.drop(df[df['SinesPerMinNew'] >= float(57.0)].index, inplace=True) #No song or sine recorded as more than 95% of recording
df.drop(df[df['SongPerMinNew'] >= float(57.0)].index, inplace=True)
df.drop(df[df['AmbientPerMin'] >= float(0.2)].index, inplace=True) # No recording with ambient/other greater than 20% of second recording
df.drop(df[df['OtherPerMin'] >= float(0.2)].index, inplace=True)

#Order and choose which columns to write to out
column_order = ['maternal_geno','phenotype','chamber','SinesPerMinNew','PulsesPerMinNew','SongPerMinNew','Sine2PulseNew','Sine2PulseMJL','Sine2PulseMJLbase',
  'SinesPerMinOld','PulsesPerMinOld','SongPerMinOld','Sine2PulseOld','ModeSineMFFT','ModePulseMFFT','SineFrameTotal','PulseFrameTotal','PulsesInTrain',
  'AllPulsesPerMin','PulseTrainsPerMin','MedianSineTrainLength', 'MedianSineAmplitudes','MedianPulseTrainLength','MedianPulseAmplitudes','ModePeak2PeakIPI',
  'ModeEnd2PeakIPI','SineTrainsPerMin','PulseTrainsPerMin','Number_slow_pulses','Number_fast_pulses','Slow2FastandSlowPulses','AmbientPerMin','OtherPerMin']

df.to_csv(outfile, na_rep='NaN', index=False, columns=column_order) #write to out
