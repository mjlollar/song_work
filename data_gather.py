#!/usr/bin/env python


### run in directory containing RIL sub dirs
### Generates combined RIL in defined directory and renames/replaces 

### Required Python libraries
import argparse
import numpy as np
import pandas as pd
import os

### Input arguments
os.getenv('filename')
parser = argparse.ArgumentParser(description='agg_data post matlab')
parser.add_argument('--i', help='Input File', type=str,required=True)
args = parser.parse_args()

in_dir = args.i
dir_num = in_dir.split('/')[0]

in_file = '<DIR>' + dir_num + '/RIL__rawdata.csv'
out_file ='<DIR>' + dir_num + '/' + dir_num + '_rawdata.csv'

### Add to dataframe RIL genotype information
df = pd.read_csv(in_file)
df.insert(0,'genotype',dir_num)

### Output updated dataframe and remove old
df.to_csv(out_file, index=False)
os.remove(in_file)

### Add to combined RIL data file
with open('<DIR>', 'a') as out_comb:
	df.to_csv(out_comb, mode='a', header=out_comb.tell()==0, index=False)
out_comb.close()
