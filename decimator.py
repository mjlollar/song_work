import scipy.signal as sp
import scipy.io as spi
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='DECIMATOR')
parser.add_argument('--i', help='Input wav path', required=True)
parser.add_argument('--r', help='New sample rate', required=True, type=int)
args = parser.parse_args()

in_file = args.i
new_rate = str(args.r)
new_rate = int(new_rate, base=16)

### Read in Wav
sample_rate, data = spi.wavfile.read(in_file)

### Calculate decimate factor
new_rate = args.r
q = int(sample_rate / new_rate)

### Downsample data
new_sample = sp.decimate(data, q)

### Output new wav
out_file = in_file.split('.')[0] + 'ds.wav'
spi.wavfile.write(out_file, rate=new_rate, data=new_sample.astype(np.int16))
