import pytest 
import os
import pandas as pd
from tcrsampler.sampler import TCRsampler


def test_TCRsampler_init():
	t = TCRsampler()

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


