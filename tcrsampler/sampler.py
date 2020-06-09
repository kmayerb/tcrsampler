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
  Class for sampling CDR3s based on specified v-gene and j-gene usage.

  Attributes
  ----------
  ref_df : pd.DataFrame
    Dataframe of reference CDR$, but contain columns ['v_reps','j_reps','cdr3', 'count','freq']

  ref_dict : dict
    dictionary keyed on tuples that point to dataframe of CDR3s

  Notes
  -----
  The goal of tcrsampler is to allow representative CDR3s to be sampled from a background set. By default, the likelihood that a CDR3 is drawn from the background is proportional to the frequency of that
  CDR3 clone in the original sample. That is, more abundant clones will be sampled more frequently than less abundant clones. 
  """
  def __init__(self):
    self.ref_df = None
    self.ref_dict = None

  def _valid_cdr3(self, cdr3):
    """
    Return True only if CDR3 is comprised of valid upper case amino acid letters.

    Parameters
    ----------
    cdr3 : str
      string repressenting amino acid sequence

    Returns
    -------
    valid : bool

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
    """
    Parameters
    ----------
    filenname : str
      name of mixcr standard output flat .tsv

    Assigns
    -------
    self.ref_df : pd.DataFrame  
      DataFarme with columns ['v_reps','j_reps','cdr3', 'count']
    """

    bar = IncrementalBar('Clean Mixcr     ', max = 1, suffix='%(percent)d%%')
    df = pd.read_csv(filename, sep = "\t")
    
    if 'subject' in df.columns:
      select_these_columns = ['bestVGene','bestJGene', 'aaSeqCDR3','cloneCount', 'cloneFraction','subject']
    else:
      select_these_columns = ['bestVGene','bestJGene', 'aaSeqCDR3','cloneCount', 'cloneFraction']
    
    rename_these_columns = {'subject':'subject',
                            'bestVGene':'v_reps',
                            'bestJGene':'j_reps', 
                            'aaSeqCDR3':'cdr3',
                            'cloneCount':'count',
                            'cloneFraction':'freq'}

    df.rename(columns =rename_these_columns, inplace = True)
    df['v_reps'] = df['v_reps'].apply(lambda x: f"{x}*01")
    df['j_reps'] = df['j_reps'].apply(lambda x: f"{x}*01")
    
    ind = df['cdr3'].str.startswith('C') & df['cdr3'].str.endswith('F')
    df = df[ind]

    ind = df['cdr3'].apply(lambda x: self._valid_cdr3(cdr3= x))
    df = df[ind]
    df = df[rename_these_columns.values()]
    if 'subject' in df.columns:
      df = df[['v_reps','j_reps','cdr3', 'count','freq','subject']].copy().reset_index(drop = True)
    else:
      df = df[['v_reps','j_reps','cdr3', 'count','freq']].copy().reset_index(drop = True)
    bar.next()
    bar.finish()
    self.ref_df = df

  def build_background(self, df = None, max_rows = 100, stratify_by_subject = False, use_frequency= True, make_singleton = False):
    """
    Parameters
    ----------
    df : pd.DataFrame
      DataFrame with ['v_reps','j_reps','cdr3', 'count','freq'] columns
    max_rows : int
      Maximum clones per v,j pair (per subject)  
    stratify_by_subject : bool
      If True, max_rows will apply to v,j,subject. If False, max_rows applies to v,j
    use_frequency : bool
      If True, uses frequency for ranking rows. If False, uses raw counts. 
    make_singleton : bool
      If True, background is still sorted by frequency or counts, but final fequency and counts values are overridden
      and set to 1. 
    Assigns
    -------
    self.ref_dict : dict
      dictionary keyed on (v,j) tuples pointing to pd.DataFrame
    """
    if df is None:
      df = self.ref_df

    if use_frequency:
      col = 'freq'
    else:
      col = 'count'

    dfg =df.groupby(['v_reps','j_reps'])
    
    bar = IncrementalBar('Build Background', max = dfg.ngroups, suffix='%(percent)d%%')
    
    d = dict()
    
    for i,group in dfg:
    
      bar.next()
      # CASE:  stratify_by_subject = False
      if not stratify_by_subject:
        if group.shape[0] > 0:
          if group.shape[0] > max_rows:
            n = max_rows
          else:
            n = group.shape[0]
          d[i] =group.sort_values([col], ascending = False).reset_index(drop = True).iloc[0:n,].copy()
      else: # CASE:  stratify_by_subject = True, another nested groupby, and max_rows applies to subject
        dfgg = group.groupby(['subject'])
        for ii , ggroup in dfgg:
          if ggroup.shape[0] > 0:
            if ggroup.shape[0] > max_rows:
              n = max_rows
            else:
              n = ggroup.shape[0]
            if i not in d.keys():
              d[i] = ggroup.sort_values([col], ascending = False).reset_index(drop = True).iloc[0:n,].copy()
            else: 
              dii = ggroup.sort_values([col], ascending = False).reset_index(drop = True).iloc[0:n,].copy()
              d[i] = pd.concat([d[i],dii])
    bar.finish()
    
    if make_singleton:
      for k,v in d.items():
        v['count'] = 1
        v['freq']  = 1

    self.ref_dict = d

  def sample_background(self,v,j,n=1, d= None, depth = 1, seed =1, use_frequency= True ):
    """
    Parameters
    ----------
    v : str
      v-gene name e.g., 'TRBV10-1*01'
    j : str
      j-gene name e.g., 'TRBJ1-1*01'
    n : int
      number of cdr3 samples to draw for given v,j
    d : dict 
      Dictionary for sampling, generated in by .build_background
    depth : 

    seed : int
      random number generating seed
    use_frequency : bool
      If True, uses frequency for sampling proportionaly. If False, uses raw counts. 

    Returns
    -------
    r: list

    Example 
    -------
    >>> sample_background(v ='TRBV10-1*01', j ='TRBJ1-1*01',n=1, d= None, depth = 1, seed =1, use_frequency= True )
    """
    if d is None:
      d = self.ref_dict


    if use_frequency:
      col = 'freq'
    else:
      col = 'count'
    
    assert isinstance(v, str)
    assert isinstance(j, str)
    assert isinstance(d, dict)
    assert isinstance(depth, int)
    assert isinstance(seed, int)

    try:
      subdf = d[(v,j)]
  
      selection_probability = \
        subdf[ col ] / np.sum(subdf[ col ])

      np.random.seed(seed) 
      
      probabalistic_selection_index = \
        np.random.choice( range(subdf.shape[0]),
        size = n * depth,
        p=selection_probability)

      r = subdf.iloc[probabalistic_selection_index,]['cdr3'].to_list()
      return r

    except KeyError:
      warnings.warn(f"({v},{j} gene usage not available")
      r = [None]
      return r


  def sample(self, v_j_usage, depth = 1, seed = 1, flatten = False, use_frequency= True):
    """
    Sample a reference dictionary based on v and j gene usage 

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
    
    Returns 
    -------
    result : list
      list of lists if flatten is False, list if flatten is True 

    Example
    -------
    >>> t.sample([['TRBV9*01','TRBJ2-7*01', 2],['TRBV7-7*01', 'TRBJ2-4*01', 4] ])
    """
    assert isinstance(v_j_usage, list)
    assert isinstance(depth, int)
    assert isinstance(seed, int)

    result = list()
    for v,j,n in v_j_usage:
      r = self.sample_background(v = v, j = j, n = n, depth = depth, seed = seed, use_frequency= use_frequency)
      result.append(r)

    if flatten:
      result = list(np.concatenate(result))

    return result

