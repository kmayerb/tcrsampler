
import os
import pandas as pd
import numpy as np
from progress.bar import IncrementalBar

def assemble_britanova_chord_blood(
	path = os.path.join('/','Volumes','Samsung_T5','kmayerbl','tcr_data','britanova'),
	chord_blood_samples = [ 'A5-S11.txt','A5-S12.txt','A5-S13.txt','A5-S14.txt',
							'A5-S15.txt','A5-S16.txt','A5-S17.txt','A5-S18.txt'],
	select_these_columns = ['v','j', 'cdr3aa','count', 'freq'],
	rename_these_columns = {'v':'v_reps',
				            'j':'j_reps', 
				            'cdr3aa':'cdr3',
				            'count':'count',
				            'freq':'freq'}):
		
	l = list()
	bar = IncrementalBar('Read Chord Blood', max = 2, suffix='%(percent)d%%')
	for f in chord_blood_samples:
		bar.next()
		assert os.path.isfile( os.path.join(path, f))
		fx =  os.path.join(path, f) 
		df = pd.read_csv(fx, sep = "\t")
		df = df[select_these_columns].\
			rename(columns = rename_these_columns)
		df['subject'] = os.path.basename(fx)
		l.append(df)
	bar.next();bar.finish()
	del df
	df = pd.concat(l)
	del l

	def _valid_cdr3(cdr3):
		amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
		valid = np.all([aa in amino_acids for aa in cdr3])
		return valid

	bar = IncrementalBar('Validate CDR3   ', max = 2, suffix='%(percent)d%%')
	bar.next()
	ind = df['cdr3'].apply(lambda x : _valid_cdr3(x))
	df = df[ind]
	bar.next()
	bar.finish()

	bar = IncrementalBar('Append *01 to V ', max = 2, suffix='%(percent)d%%')
	bar.next()
	v_reps = df['v_reps'].apply(lambda x: f"{x}*01")
	df = df.assign(v_reps = v_reps )
	bar.next();	bar.finish()

	bar = IncrementalBar('Append *01 to J ', max = 2, suffix='%(percent)d%%')
	bar.next()
	j_reps = df['j_reps'].apply(lambda x: f"{x}*01")
	df = df.assign(j_reps = j_reps )
	bar.next();	bar.finish()

	bar.finish()
	return df


if __name__ == '__main__':
	df = assemble_britanova_chord_blood()
	df.to_csv("britanova_chord_blood.csv")
	df[['v_reps', 'j_reps', 'cdr3', 'count', 'freq', 'subject']].to_csv("britanova_human_beta_t_cb.tsv.sampler.tsv", sep = "\t", index = False)
