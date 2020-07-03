
import os
import pandas as pd
import numpy as np
from progress.bar import IncrementalBar
import sys

def _valid_cdr3(cdr3):
	if not isinstance(cdr3, str):
		return False
	else:
		amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
		valid = np.all([aa in amino_acids for aa in cdr3])
		return valid

def fix_adaptive(dft):
	l = list()
	# Deal with ties, make a row for all possible combindations
	for i, row in dft.iterrows():
			vs = list()
			js = list()
			if row['v_gene'] == 'unresolved':
					vs =row['v_gene_ties'].split(",")
			else:
					vs = [row['v_gene']]
			if row['j_gene'] == 'unresolved':
					js = row['j_gene_ties'].split(",")
			else:
					js = [row['j_gene']]
			for v in vs:
					for j in js:
						 row['v_gene'] = v
						 row['j_gene'] = j
					l.append(row.to_list())

	df = pd.DataFrame(l, columns = dft.columns )

	# fix adaptive agravating TCR naming convention to match IMGT
	def _fix_v(v_gene):
		v_gene = v_gene.replace("TCRB", "TRB").\
			replace("-0", "-").replace("BV0","BV")
		if v_gene.find('*0') == -1:
			v_gene = v_gene + "*01"
		return v_gene
	 
	def _fix_j(j_gene):	
		j_gene = j_gene.replace("TCRB", "TRB").\
			replace("-0", "-").replace("BJ0","BJ")
		if j_gene.find('*0') == -1:
				j_gene = j_gene + "*01"
		return j_gene

	f_v = df['v_gene'].apply(lambda x : _fix_v(x))
	f_j = df['j_gene'].apply(lambda x : _fix_j(x))
	df = df.assign(v_gene = f_v)
	df = df.assign(j_gene = f_j)
	df.reset_index(drop = True)

	return df

def assemble_emerson_cmv_negative(
	path = os.path.join('/','Volumes','Samsung_T5','kmayerbl','tcr_data','emerson','cohort2'),
	samples = [ 'Keck0002_MC1.tsv',
				'Keck0014_MC1.tsv',
				'Keck0017_MC1.tsv',
				'Keck0019_MC1.tsv',
				'Keck0042_MC1.tsv',
				'Keck0114_MC1.tsv',
				'Keck0116_MC1.tsv'],
	select_these_columns =	['sample_name','v_gene','v_gene_ties','j_gene','j_gene_ties', 'cdr3_amino_acid','frequency', 'templates']):
	bar = IncrementalBar('Clean Emerson CMV-', max = len(samples), suffix='%(percent)d%%')
	lx = list()
	for f in samples:
		bar.next()
		
		assert os.path.isfile( os.path.join(path, f))
		fx =os.path.join(path, f)
		
		dt = {'sample_name':str,'v_gene':str,'v_gene_ties':str,'j_gene':str,'j_gene_ties':str, 'cdr3_amino_acid':str,'frequency':'Float64', 'templates':'Int64'}
		df = pd.read_csv( fx ,sep = "\t", usecols= select_these_columns, dtype =dt)

		df = fix_adaptive(dft = df)
		#print(df.head(2))
		rename_these_columns = {'sample_name':'subject',
								'v_gene':'v_reps',
								'j_gene':'j_reps', 
								'cdr3_amino_acid':'cdr3',
								'templates':'count',
								'frequency': 'freq'}
		df.rename(columns =rename_these_columns, inplace = True)
		#print(df.head(2))
		
		ind = df['cdr3'].apply(lambda x : _valid_cdr3(x))
		rows_pre = df.shape[0] # rows before checking if cdr3 are actually valid
		df = df[ind]
		sys.stdout.write(f" Adding {df.shape[0]} of {rows_pre} rows from {f}\n")
		lx.append(df[['subject','v_reps','j_reps','cdr3','count','freq']])
	
	bar.finish()
	return pd.concat(lx).reset_index(drop = True)


df = assemble_emerson_cmv_negative()
bar = IncrementalBar('Write .to_csv.    ', max = 2, suffix='%(percent)d%%')
bar.next()
#df.to_csv("emerson_cmv_negative.csv", index = False)
df[['v_reps', 'j_reps', 'cdr3', 'count', 'freq', 'subject']].to_csv("emerson_human_beta_t_cmvneg.tsv", sep = "\t", index = False)
bar.next();bar.finish()

