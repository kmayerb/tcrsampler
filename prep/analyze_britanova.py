import os
import pandas as pd
from tcrsampler.sampler import TCRsampler


t = TCRsampler()
fn= os.path.join('britanova_chord_blood.csv')
t.ref_df = pd.read_csv(fn)
t.build_background( max_rows = 1000)
t.sample([['TRBV10-2*01', 'TRBV10-2*01*01', 1],['TRBV27*01', 'TRBV27*01*01', 4] ], depth = 10)

for k,v in t.ref_dict.items():
	print (k, v.shape[0])
	