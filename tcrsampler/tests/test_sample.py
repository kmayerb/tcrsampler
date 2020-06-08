import pytest
from tcrsampler.sampler import TCRsampler
import pandas as pd
import numpy as np

t = TCRsampler(organism = 'human',db_file = 'alphabeta_db.tsv', ref_file= 'new_nextgen_chains_human_B.tsv')
t.make_ref_dict()

def test_TCRsampler_init_w_db_file():
	t = TCRsampler(organism = 'human', db_file = 'alphabeta_db.tsv', ref_file= 'new_nextgen_chains_human_B.tsv')

def test_TCRsampler_init_w_db_file_KeyError():
	with pytest.raises(KeyError):
			t = TCRsampler(organism = 'hamster', db_file = 'alphabeta_db.tsv', ref_file= 'new_nextgen_chains_human_B.tsv')

def test_TCRsample():


	r = t.sample(vj_usage = [('TRBV7-9*05', 'TRBJ1-3*01', 3), ('TRBV7-9*05', 'TRBJ1-4*01', 10)], multiple= 1)
	assert isinstance(r, list)
	assert r[0] == ['CASRPRLGVGAPIVSGNTIYF', 'CASRPRLGVGAPIVSGNTIYF', 'CASRPRLGVGAPIVSGNTIYF']
	assert r[1] == ['CASNPGGGNEKLFF', 'CASRGSPNEKLFF', 'CASSVTSKNEKLFF', 'CASCCRDSPSGKLFF', 'CASTGDNEKLFF', 'CASTGDNEKLFF', 'CASSRDSLLGEKLFF', 'CASSVTSKNEKLFF', 'CAIRAGTPPEKLFF', 'CAIRAGTPPEKLFF']

def test_TCRsample_invalid_v():

	with pytest.warns(None) as w:
		r = t.sample(vj_usage = [('TRGV-9*05', 'TRBJ1-3*01', 3), ('TRBV7-9*05', 'TRBJ1-4*01', 10)], multiple= 1)
	assert len(w) >=1


def test_TCRsample_invalid_j():

	with pytest.warns(None) as w:
		r = t.sample(vj_usage = [('TRBV7-9*05', 'TRBX1-3*01', 3), ('TRBV7-9*05', 'TRBJ1-4*01', 10)], multiple= 1)
	assert len(w) >=1

def test_TCRsample_from_DataFrame():
	df = pd.DataFrame( { "v_gene":['TRBV7-9*05','TRBV7-9*05','TRBV7-9*05'],
						"j_gene":['TRBJ1-3*01','TRBJ1-4*01','TRBJ1-4*01']})

	# efficient way to generate a set of synethetic seqs with same gene frequency as input
	
	gene_usage = df.groupby(['v_gene','j_gene']).size()
	
	# v_gene      j_gene
	# TRBV7-9*05  TRBJ1-3*01    1
 	#          	  TRBJ1-4*01    2

	syn_seq = t.sample( gene_usage.reset_index().to_dict('split')['data'], multiple = 1)
	[['CASRPRLGVGAPIVSGNTIYF'], ['CASNPGGGNEKLFF', 'CASRGSPNEKLFF']]
	assert list(np.concatenate(syn_seq)) == ['CASRPRLGVGAPIVSGNTIYF', 'CASNPGGGNEKLFF', 'CASRGSPNEKLFF']

	# Or add a simulated cdr3 as a column to a dataframe
	df['sampledcdr3'] = df.apply(lambda x : t.sample_ref_dict(v=x['v_gene'],j=x['j_gene'],n=1)[0], axis = 1)
	df['sampledcdr3'].to_list() == ['CASRPRLGVGAPIVSGNTIYF', 'CASNPGGGNEKLFF', 'CASRGSPNEKLFF']

def test_TCRsample_from_DataFrame_flatten():
	df = pd.DataFrame( { "v_gene":['TRBV7-9*05','TRBV7-9*05','TRBV7-9*05'],
						"j_gene":['TRBJ1-3*01','TRBJ1-4*01','TRBJ1-4*01']})

	# efficient way to generate a set of synethetic seqs with same gene frequency as input
	
	gene_usage = df.groupby(['v_gene','j_gene']).size()
	syn_seq = t.sample( gene_usage.reset_index().to_dict('split')['data'], multiple = 1, flatten = True)
	assert syn_seq == ['CASRPRLGVGAPIVSGNTIYF', 'CASNPGGGNEKLFF', 'CASRGSPNEKLFF']