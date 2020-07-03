import pytest 
import os
import pandas as pd
from tcrsampler.sampler import TCRsampler
import numpy as np

def test_TCRsampler_init():
	t = TCRsampler()

def test_TCRsampler_init_default():
	t = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
	assert t.default_bkgd == 'britanova_human_beta_t_cb.tsv.sampler.tsv'
	assert isinstance(t.ref_df, pd.DataFrame)
	assert isinstance(t.ref_dict, dict)
	assert 'TRBV2*01' in t.v_freq.keys()

def test_TCRsampler_clean_mixcr():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	assert isinstance(t.ref_df, pd.DataFrame)

def test_TCRsampler_build():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	t.build_background()
	assert isinstance(t.ref_dict, dict)
	assert isinstance(t.ref_dict.popitem()[1], pd.DataFrame)

def test_TCRsampler_build_vj_components():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	t.build_background()
	assert np.isclose(np.sum([k for _,k in t.vj_freq.items()]), 1.0)
	assert np.isclose(np.sum([k for _,k in t.j_freq.items()]), 1.0)
	assert np.isclose(np.sum([k for _,k in t.v_freq.items()]), 1.0)
	assert np.isclose(np.sum([k for _,k in t.vj_occur_freq.items()]), 1.0)
	assert np.isclose(np.sum([k for _,k in t.v_occur_freq.items()]), 1.0)
	assert np.isclose(np.sum([k for _,k in t.j_occur_freq.items()]), 1.0)
	
def test_TCRsampler_build_singleton():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	t.build_background(stratify_by_subject = True, make_singleton = True)
	r = t.sample_background('TRBV9*01','TRBJ2-7*01', n =10 )


def test_TCRsampler_build_stratified():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	t.build_background(stratify_by_subject = True)
	r = t.sample_background('TRBV9*01','TRBJ2-7*01', n =10 )


def test_prob_sampler_sample_background():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	t.build_background()
	r = t.sample_background('TRBV9*01','TRBJ2-7*01', n =10 )
	assert r == ['CASSRTGSLADEQYF',
				 'CASSATGVVSAQYF',
				 'CASSAWGQVYEQYF',
				 'CASSVSGSPYEQYF',
				 'CASSAWGQVYEQYF',
				 'CASSAWGQVYEQYF',
				 'CASRWGEQYF',
				 'CASSGDDWEQYF',
				 'CASSATGTSGPYEQYF',
				 'CASSSRTSGSNSEQYF']

def test_prob_sampler_sample():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	t.build_background()
	r = t.sample([['TRBV9*01','TRBJ2-7*01', 2]])
	assert r == [['CASSRTGSLADEQYF', 'CASSATGVVSAQYF']]
	r = t.sample([['TRBV9*01','TRBJ2-7*01', 2]], flatten = True)
	assert r == ['CASSRTGSLADEQYF', 'CASSATGVVSAQYF']
	r = t.sample([['TRBV9*01','TRBJ2-7*01', 2],['TRBV7-7*01', 'TRBJ2-4*01', 4] ])
	assert r == [['CASSRTGSLADEQYF', 'CASSATGVVSAQYF'],
 					['CASSLGQAARGIQYF', 'CASSLGQAARGIQYF', 'CASSLGQAARGIQYF', 'CASSLGQAARGIQYF']]


def test_prob_sampler_sample_key_warn():
	t = TCRsampler()
	fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	t.clean_mixcr(filename = fn)
	t.build_background()
	with pytest.warns(None):
		r = t.sample([['TRBV999*01','TRBJ2-7*01', 2]])
	assert r == [[None]]

def test_v_j_freq_estimates():
	d = {'Unnamed: 0': {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}, 'v_reps': {0: 'TRBV24-1*01',  1: 'TRBV5-1*01',  2: 'TRBV7-2*01',  3: 'TRBV3-1*01',  4: 'TRBV7-3*01'}, 'j_reps': {0: 'TRBJ2-1*01',  1: 'TRBJ2-5*01',  2: 'TRBJ2-3*01',  3: 'TRBJ2-5*01',  4: 'TRBJ2-3*01'}, 'cdr3': {0: 'CATRQDNEQFF',  1: 'CASSLEETQYF',  2: 'CASSLADTQYF',  3: 'CASSQETQYF',  4: 'CASSLAGGTDTQYF'}, 'count': {0: 252, 1: 166, 2: 113, 3: 98, 4: 89}, 'freq': {0: 0.0003726818302818776,  1: 0.0002454967612174273,  2: 0.00016711526516608003,  3: 0.00014493182288739684,  4: 0.00013162175752018694}, 'subject': {0: 'A5-S11.txt',  1: 'A5-S11.txt',  2: 'A5-S11.txt',  3: 'A5-S11.txt',  4: 'A5-S11.txt'}}
	df = pd.DataFrame(d)
	t = TCRsampler()
	t.ref_df = df
	t.build_background()
	assert  t.v_occur_freq == {  'TRBV3-1*01': 0.2,
								 'TRBV5-1*01': 0.2,
								 'TRBV7-2*01': 0.2,
								 'TRBV7-3*01': 0.2,
								 'TRBV24-1*01': 0.2}
	assert t.j_occur_freq == {'TRBJ2-1*01': 0.2, 'TRBJ2-3*01': 0.4, 'TRBJ2-5*01': 0.4}

								
def test_classmethod_calls():
	import tcrsampler as ts
	ts.TCRsampler.currently_available_backgrounds()
	ts.TCRsampler.download_all_background_files(dry_run = True)
	ts.TCRsampler.download_background_file(download_file = 'wiraninha_sampler.zip', dry_run = True, download_from = "dropbox")
	ts.TCRsampler.get_download_address()

def test_currently_available_backgrounds():
	import tcrsampler as ts
	r = ts.TCRsampler.currently_available_backgrounds()
	assert isinstance(r, list)

def test_get_download_address():
	import tcrsampler as ts
	r = ts.TCRsampler.get_download_address()
	assert isinstance(r, list)

	select_files = ["wiraninha_sampler.zip", 
	                "ruggiero_mouse_sampler.zip", 
	                "ruggiero_human_sampler.zip",
	                "emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip",
	                "britanova_human_beta_t_cb.tsv.sampler.tsv.zip"]
	for sf in select_files:
		assert sf in r





