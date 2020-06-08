import pytest
from tcrsampler.sampler import TCRsampler
import pandas as pd
import os
import numpy as np

def test_cleaner():
	t = TCRsampler(organism = 'human',db_file = 'alphabeta_db.tsv', ref_file= 'new_nextgen_chains_human_B.tsv')
	filename = os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
	reference = t.clean_mixcr(filename)
	t.make_ref_dict(df = reference, max_sample = 3 )
	assert isinstance(t.ref_dict[('TRBV1*01', 'TRBJ1-1*01')], list)
	assert isinstance(t.sample([['TRBV1*01', 'TRBJ1-1*01',1]]), list)
