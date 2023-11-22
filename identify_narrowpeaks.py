#MJL last edit 09/10/23

import argparse
import numpy as np
import pandas as pd

parser= argparse.ArgumentParser()
parser.add_argument('--i', type=str, help="all window input file", required=True)
parser.add_argument('--p', type=str, help="peak window input file", required=True)
args = parser.parse_args()

df_allwin = pd.read_csv(args.i, sep='\t')
df_peaks = pd.read_csv(args.p, sep='\t')

#Get Peak Conf start/stop sites, highest peak
conf_winstarts = df_peaks['lowerStartCI'].tolist()
conf_winstops = df_peaks['upperEndCI'].tolist()
peak_peak = df_peaks['Start'].tolist()
peak_chroms = df_peaks['Chr'].tolist()
assert len(conf_winstarts) == len(conf_winstops) == len(peak_peak) == len(peak_chroms) #sanity

#Initialize some lists
conf_finalstart = []
conf_finalend = []
narrow_CI_applied = []

#Loop through each peak identified in all_peaks output
for peak in range(0, len(conf_winstarts)):
	# Get window number length of basic conf
	peak_num = peak + 1
	index_start = df_allwin.query('Start == @conf_winstarts[@peak] and Chr == @peak_chroms[@peak]').index[0] #3738
	index_end = df_allwin.query('End == @conf_winstops[@peak] and Chr == @peak_chroms[@peak]').index[0] #3775
	simple_CI_length = index_end - index_start + 1 #Length of Max LOD - 1.5 interval  #38

	#Get window number of area above pvalue
	peak_start = df_allwin.query('Start == @peak_peak[@peak] and Chr == @peak_chroms[@peak]').index[0] #3755
	thresh_wins = 1 # Start at 1 to count peak win
	i = 1
	while True: #forward
		row_id = int(peak_start + i)
		wintype = df_allwin.loc[df_allwin.index[row_id], 'group']
		if wintype == 'P':
			thresh_wins = thresh_wins + 1
			i = i + 1
			continue
		elif wintype == 'V':
			thresh_CI_end = row_id - 1
			break
		else:
			print("something is wrong") #sanity
			break
	i = 1 #reverse
	while True:
		row_id = int(peak_start - i)
		wintype = df_allwin.loc[df_allwin.index[row_id], 'group']
		if wintype == 'P':
			thresh_wins = thresh_wins + 1
			i = i + 1
			continue
		elif wintype == 'V':
			thresh_CI_start = row_id + 1
			break
		else:
			print("something is wrong") #sanity
			break

	if simple_CI_length >= thresh_wins: #MaxLOD - 1.5 window
		print("MaxLOD - 1.5 C.I. most inclusive for peak " + str(peak_num) + " in set")
		conf_finalstart.append(conf_winstarts[peak])
		conf_finalend.append(conf_winstops[peak])
		narrow_CI_applied.append('No')

	else:
		print("Narrow Peaks Identified for peak " + str(peak_num) + " in set")
		narrow_CI_id = []
		CI_bases_start = []
		CI_bases_stop = []

		# Get MaxLOD - 1.5 pvalue of peak window
		lod_threshold = float(df_peaks.loc[df_peaks.index[peak], 'LODnum']) - 1.5 #LOD threshold for peak #1.36403...
		sig_interval = list(range(thresh_CI_start, thresh_CI_end+1)) #Total area above significance to evaluate (3753, 3755 + 1)
		start_bases = df_allwin.loc[df_allwin.index[thresh_CI_start],'Start']
		stop_bases = df_allwin.loc[df_allwin.index[thresh_CI_end],'End']
		conf_finalstart.append(start_bases)
		conf_finalend.append(stop_bases)
		narrow_CI_applied.append('Yes')

		#Id peak or valley
		for win_index in sig_interval:
			lod_value = float(df_allwin.loc[df_allwin.index[win_index], 'LODnum'])
			CI_bases_start.append(df_allwin.loc[df_allwin.index[win_index],'Start'])
			CI_bases_stop.append(df_allwin.loc[df_allwin.index[win_index],'End'])
			if lod_value >= lod_threshold:
				narrow_CI_id.append('P')
			else:
				narrow_CI_id.append('V')

		#Separate output file for criteria B matching results
		df_out = pd.DataFrame({'Start': CI_bases_start, 'Stop': CI_bases_stop, 'ID': narrow_CI_id})
		out_name = args.i + '_narrow_ID_peakNum_' + str(peak_num) + '.csv'
		df_out.to_csv(out_name, sep=',', index=False)

### output total results
df_peaks['FinalCIStart'] = conf_finalstart
df_peaks['FinalCIEnd'] = conf_finalend
df_peaks['Narrow_applied'] = narrow_CI_applied
outname = args.i + '_MJL_addition.csv'
df_peaks.to_csv(outname, sep=',', index=False)
