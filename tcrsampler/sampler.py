import os
import sys
import pandas as pd
import numpy as np
import warnings
import random 
import time
from progress.bar import IncrementalBar

class TCRsampler():
	"""
	A Python Class for sampling a known CDR3 TCR repertoire matching a 
	a known frequency of V-gene and J-gene usages. This class, can 
	either use default or custom reference sets.
	

	Attributes
	----------
	organism : str
		must be 'human' or 'mouse'
	db_path : str or None
		path to database (can be left blank for default db)
	db_file : str
		filename of database ('alphabeta_db.tsv')
	ref_path:
		path to database (can be left blank for default db)
	ref_file : str
		filenames of refernce CDR3s

	Methods 
	------_
	make_ref_dict : method
	
	sample_ref_dict : method 
	
	sample : method

	Example
	-------
	>>> from tcrsampler.sampler import TCRsampler
	>>> t = TCRsampler(organism = 'human',db_file = 'alphabeta_db.tsv', ref_file= 'new_nextgen_chains_human_B.tsv')
	>>> t.make_ref_dict()
	>>> t.sample(vj_usage = [('TRBV7-9*05', 'TRBJ1-3*01', 3), ('TRBV7-9*05', 'TRBJ1-4*01', 10)], multiple= 1)
	[['CASRPRLGVGAPIVSGNTIYF', 'CASRPRLGVGAPIVSGNTIYF', 'CASRPRLGVGAPIVSGNTIYF'],
	['CASNPGGGNEKLFF', 'CASRGSPNEKLFF', 'CASSVTSKNEKLFF', 'CASCCRDSPSGKLFF', 'CASTGDNEKLFF', 'CASTGDNEKLFF', 'CASSRDSLLGEKLFF', 'CASSVTSKNEKLFF', 'CAIRAGTPPEKLFF', 'CAIRAGTPPEKLFF']]
	
	Notes
	-----

	The db_file specifies all possible TCR V and J genes. 

	The ref_file is a reference of CDR3s

	For alpha-chain:

	db_file = "alphabeta_db.tsv", ref_file= 'new_nextgen_chains_human_A.tsv'

	For beta-chain:
	
	db_file = "alphabeta_db.tsv", ref_file= 'new_nextgen_chains_human_B.tsv'

	"""
	def __init__(self,
				organism ='human',
				db_path=None, 
				db_file=None, 
				ref_path=None, 
				ref_file=None):
		"""
		Parameters
		----------
		organism : str
			must be 'human' or 'mouse'
		db_path : str or None
			path to database (can be left blank for default db)
		db_file : str
			filename of database ('alphabeta_db.tsv')
		ref_path:
			path to database (can be left blank for default db)
		ref_file : str
			filenames of refernce CDR3s
		"""
		valid_organisms = ['human','mouse']
		if organism not in valid_organisms:
			raise KeyError(f'You must specify an organism in {valid_organisms}')
		self.organism = organism
		self.db_path = db_path
		self.ref_path = ref_path
		
		if self.db_path is None:
			self.db_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"db")
		if self.ref_path is None:
			self.ref_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"db")
		
		self.db_file = db_file
		self.ref_file = ref_file

		self.db = os.path.join(self.db_path, self.db_file)
		assert os.path.isfile(self.db)
		
		self.ref = os.path.join(self.ref_path, self.ref_file)
		assert os.path.isfile(self.ref)

		self.db_df = pd.read_csv(self.db, sep= '\t')
		self.db_df = self.db_df[self.db_df['organism'] == self.organism]

		self.ref_df = pd.read_csv(self.ref, sep = '\t')


	def make_ref_dict(self, df = None, max_sample = 10):
		"""
		Parameters
		----------
		df : DataFrame
			reference DataFrame (with columns: v_reps, j_reps, and cdr3)
		max_sample : int 
			maximum number of CDR3s for each V,J gene usage combination CDR3 
		
		Assigns
		-------
		self.ref_dict keyed on (v_gene,j_gene) points to list of CDR3s
		"""
		d = dict()
		if df is None:
			df = self.ref_df

		df_grouped= df.groupby(['v_reps','j_reps'])
		
		bar = IncrementalBar('Making Reference', max = df_grouped.ngroups, suffix='%(percent)d%%')

		for group_name, df_group in df_grouped:
			bar.next()
			v,j = group_name
			for v in sorted(v.split(",")):
				for j in sorted(j.split(",")):
					if v not in self.db_df.id.to_list():
						warnings.warn(f"{v} was not recognized in the TCRsample.db_df")
					if j not in self.db_df.id.to_list():
						warnings.warn(f"{j} was not recognized in the TCRsample.db_df")
						
					d[(v,j)] = list()
						
					for row_index, row in df_group.reset_index(drop = True).iterrows():
						if row_index < max_sample:
							
							cdr3 = row['cdr3']
							d[(v,j)].append(cdr3)
						else:
							break
		bar.finish()
		sys.stdout.write("TCRsampler reference dictionary created.\n")
		sys.stdout.write("Sample a CDR3 by .sample_ref_dict('TRBV1*01', 'TRBJ1-1*01',1)\n")
		sys.stdout.write("Sample a CDR3 repertoire by .sample([['TRBV1*01', 'TRBJ1-1*01',2], ['TRBV2*01', 'TRBJ1-1*01',1]])\n")
		
		self.ref_dict = d

	def sample_ref_dict(self, v,j, n = 1, seed =1):
		"""
		Sample a list of length n for a given v,j pair 

		Paramaters
		----------
		v : str
			v-gene name
		j : str
			j-gene name
		n : int 
			number to sample
		seed : int 
			random seed 

		Returns
		-------
		smp : list of sampled cdr3s of length n

		Example
		-------
		>>> from tcrsampler.sampler import TCRsampler
		>>> t = TCRsampler(organism = 'human',db_file = 'alphabeta_db.tsv', ref_file= 'new_nextgen_chains_human_B.tsv')
		>>> t.make_ref_dict()
		>>> t.sample_ref_dict('TRBV7-9*05','TRBJ1-3*01', 1)
		['CASRPRLGVGAPIVSGNTIYF']
		>>> t.sample_ref_dict('TRBV7-9*05','TRBJ1-4*01', 3)
		['CASNPGGGNEKLFF', 'CASRGSPNEKLFF', 'CASSVTSKNEKLFF']

		"""
		try:
			random.seed(seed)
			l = self.ref_dict[(v,j)]
			smp = random.choices(l, k =n)
			return smp 
		except KeyError:
			warnings.warn(f"({v},{j}) was not recognized in the TCRsample.ref_dict")
			return None

	def sample(self, vj_usage:tuple, multiple:int = 1, flatten = False):
		"""
		Sample multiple list of cdr3s of

		To sample more deeply increase multiple from 1.
		
		Parameters
		----------
		vj_usage : list of tuples or lists
			e.g., [('TRBV7-9*05', 'TRBJ1-3*01', 3), ('TRBV7-9*05', 'TRBJ1-4*01', 10)], multiple= 10)
		multiple : int
			Depth of sampling.
		flatten : bool 
			If true return a flat list.

		Returns
		-------

		r : list of lists 

		Example
		-------
		>>> t.sample(vj_usage = [('TRBV7-9*05', 'TRBJ1-3*01', 3), ('TRBV7-9*05', 'TRBJ1-4*01', 10)], multiple= 1)
		[['CASRPRLGVGAPIVSGNTIYF', 'CASRPRLGVGAPIVSGNTIYF', 'CASRPRLGVGAPIVSGNTIYF'],
		['CASNPGGGNEKLFF', 'CASRGSPNEKLFF', 'CASSVTSKNEKLFF', 'CASCCRDSPSGKLFF', 'CASTGDNEKLFF', 'CASTGDNEKLFF', 'CASSRDSLLGEKLFF', 'CASSVTSKNEKLFF', 'CAIRAGTPPEKLFF', 'CAIRAGTPPEKLFF']]

		"""

		r = [self.sample_ref_dict(v,j,multiple*n) for v,j,n in vj_usage]
		if flatten:
			return list(np.concatenate(r))
		else:
			return r

	def _valid_cdr3(self, cdr3):
		"""
		Examples
		--------
		>>> _valid_cdr3("AAAA")
		True
		>>> _valid_cdr3("AA.A")
		False
		"""
		amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
		valid = np.all([aa in amino_acids for aa in cdr3])
		return valid


	def clean_mixcr(self, filename, approximate = True):
		df = pd.read_csv(filename, sep = "\t")
		
		select_these_columns = ['bestVGene','bestJGene', 'aaSeqCDR3','cloneCount']
		rename_these_columns = {'bestVGene':'v_reps',
								'bestJGene':'j_reps', 
								'aaSeqCDR3':'cdr3',
								'cloneCount':'count'}

		df.rename(columns =rename_these_columns, inplace = True)
		df['v_reps'] = df['v_reps'].apply(lambda x: f"{x}*01")
		df['j_reps'] = df['j_reps'].apply(lambda x: f"{x}*01")
		
		ind = df['cdr3'].str.startswith('C') & df['cdr3'].str.endswith('F')
		df = df[ind]

		ind = df['cdr3'].apply(lambda x: self._valid_cdr3(cdr3= x))
		df = df[ind]

		df = df[rename_these_columns.values()]
		l = list()
		if approximate:
			# this divides each group by the group min to reduced overall size

			dfg = df.groupby(['v_reps','j_reps'])
			bar = IncrementalBar('Processing Mixcr', max = dfg .ngroups, suffix='%(percent)d%%')

			for group_name, df_group in dfg:
				bar.next()
				group_min= df_group['count'].min()
				#df_group['countaprox'] = df_group['count'].copy().apply(lambda x : x // group_min)
				x =  df_group['count'].apply(lambda x : x // group_min)
				#l.append(df_group.reset_index(drop = True))
				l.append(x)
		#	df = pd.concat(l).reset_index()
			df['count'] = np.concatenate(l)
			del dfg
			bar.finish()
		# expand
		ndf = df.reindex(df.index.repeat(repeats = df['count']))
		assert df['count'].sum() == ndf.shape[0]
		del df
		# shuffle

		ndf = ndf.sample(ndf.shape[0]).reset_index(drop = True)
		return ndf[['v_reps','j_reps','cdr3']]
