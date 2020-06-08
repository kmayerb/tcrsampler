import os
import sys
import pandas as pd
import numpy as np
import warnings
import random 
import time
from progress.bar import IncrementalBar

class TCRsampler():
  def __init__(self):
    self.ref_df = None
    self.ref_dict = None

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

  def clean_mixcr(self, filename):

    bar = IncrementalBar('Clean Mixcr     ', max = 1, suffix='%(percent)d%%')
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
    df = df[['v_reps','j_reps','cdr3', 'count']].copy().reset_index(drop = True)
    bar.next()
    bar.finish()
    self.ref_df = df

  def clean_reps_format():
    pass


  def build_background(self, df = None, max_rows = 100):
    """
    Parameters
    ----------
    d: dict 

    Assigns
    -------
    self.ref_dict : dict
      dictionary keyed on (v,j) tuples pointing to pd.DataFrame
    """
    if df is None:
      df = self.ref_df
    dfg =df.groupby(['v_reps','j_reps'])
    d = dict()
    bar = IncrementalBar('Build Background', max = dfg.ngroups, suffix='%(percent)d%%')
    for i,group in dfg:
      bar.next()
      if group.shape[0] > 0:
        if group.shape[0] > max_rows:
          n = max_rows
        else:
          n = group.shape[0]
        d[i] =group.sort_values(['count'], ascending = False).reset_index(drop = True).iloc[0:n,].copy()

    bar.finish()
    self.ref_dict = d

  def sample_background(self,v,j,n=1, d= None, depth = 1, seed =1):
    """
    Parameters
    ----------
    d: dict 
      defaults to self.ref_dict

    Returns
    -------
  

    """
    if d is None:
      d = self.ref_dict
    
    assert isinstance(v, str)
    assert isinstance(j, str)
    assert isinstance(d, dict)
    assert isinstance(depth, int)
    assert isinstance(seed, int)

    try:
      subdf = d[(v,j)]
  
      selection_probability = \
        subdf['count'] / np.sum(subdf['count'])

      np.random.seed(seed) 
      
      probabalistic_selection_index = \
        np.random.choice( range(subdf.shape[0]),
        size = n * depth,
        p=selection_probability)

      return subdf.iloc[probabalistic_selection_index,]['cdr3'].to_list()

    except KeyError:
      warnings.warn(f"({v},{j} gene usage not available")
      return [None]


  def sample(self, v_j_usage, depth = 1, seed = 1, flatten = False):
    """
    Parameters
    ----------
    v_j_usage : list of lists or list of tuples (v-gene. j-gene, n)
      e.g., ['TRBV9*01','TRBJ2-7*01', 2],['TRBV7-7*01', 'TRBJ2-4*01', 4] 
    depth : int
      mulitple of the number of times to sample for a given frequncy. 
    seed : int
      random number initialization
    flatten : bool
      if true return a single list, if false a list of lists 

    Example
    -------
    >>> t.sample([['TRBV9*01','TRBJ2-7*01', 2],['TRBV7-7*01', 'TRBJ2-4*01', 4] ])
    """
    assert isinstance(v_j_usage, list)
    assert isinstance(depth, int)
    assert isinstance(seed, int)

    result = list()
    for v,j,n in v_j_usage:
      r = self.sample_background(v = v, j = j, n = n, depth = depth, seed = seed)
      result.append(r)

    if flatten:
      result = list(np.concatenate(result))

    return result





