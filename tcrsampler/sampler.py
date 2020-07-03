import os
import sys
import pandas as pd
import numpy as np
import warnings
import random 
import time
from progress.bar import IncrementalBar
from tcrsampler.setup_db import install_all_next_gen, install_nextgen_data_to_db, select_files

__all__ = ['TCRsampler']

class TCRsampler():
  """ 
  Class for sampling CDR3s based on specified v-gene and j-gene usage.

  Attributes
  ----------
  default_background: str or None
    string name of background file, expected in db/ directory of package source if setup_db step are run.
  ref_df : pd.DataFrame
    Dataframe of reference CDR$, but contain columns ['v_reps','j_reps','cdr3', 'count','freq']
  ref_dict : dict
    dictionary keyed on tuples that point to dataframe of CDR3s

  vj_freq : dict
    dictionary keyed on  V,J gene name tuples pointing to frequency of V,J-gene pairings by FREQUENCY METHOD
  v_freq : dict
    dictionary keyed on  V gene name pointing to frequency of V-gene pairings by FREQUENCY METHOD
  j_freq : dict
    dictionary keyed on  J gene name pointing to frequency of J-gene pairings by FREQUENCY METHOD
  vj_occur_freq : dict
    dictionary keyed on tuples pointing to frequency of v,j pairings by OCCURRENCE METHOD
  v_occur_freq : dict
    dictionary keyed on  V gene name pointing to frequency of V-gene pairings by OCCURRENCE METHOD 
  j_occur_freq : dict 
    dictionary keyed on  J gene name pointing to frequency of J-gene pairings by OCCURRENCE METHOD 

  Notes
  -----
  The goal of tcrsampler is to allow representative CDR3s to be sampled from a background set. By default, the likelihood that a CDR3 is drawn from the background is proportional to the frequency of that
  CDR3 clone in the original sample. That is, more abundant clones will be sampled more frequently than less abundant clones. 
  
  FREQUENCY METHOD: V_J OCCURRENCE PROBABILITIES BY THE SEQUENCE FREQUENCY METHOD
  By this method, all seqs in the input are considered and probability of a V,J pairing
  is assumed to be proportional to frequency of V,J-derived sequences in the sample. This method has 
  the benefit of using all the data, but it is influenced by clonal expansions.

  OCCURRENCE METHOD: V_J OCCURRENCE PROBABILITIES BY THE UNIQUE N CLONES METHOD
  By this method, we take the top N clones from each subject and 
  V,J pairing probability is assumed to be proportional to frequncy of unique 
  clones. This method has a disadvantage that we only consider the number of 
  clones per sample equal to that in the least diverse sample; 
  however this method is not biased by clonal expansion.

  N is currently determined by taking the number of clones in the least diverse subject 
  (i.e., the subject with the fewest clones). In this case the ~102,000 clones are considered 
  from each sample. 

  Notes:
  Default data from Britanova OV, Shugay M, Merzlyak EM, Staroverov DB, Putintseva EV, Turchaninova MA, Mamedov IZ, Pogorelyy MV, Bolotin DA, Izraelson M, et al. Dynamics of individual T cell repertoires: from cord blood to centenarians. J Immunol. 2016;196:5005–5013. doi: 10.4049/jimmunol.1600005.

  """
  def __init__(self, default_background = None):
    self.default_bkgd = default_background
    self.ref_df = None
    self.ref_dict = None

    if default_background is not None:
      path_to_db = os.path.join(os.path.dirname(os.path.realpath(__file__)),'db')
      path_to_db_bkgd = os.path.join(path_to_db, self.default_bkgd)
      if not os.path.isfile(path_to_db_bkgd):
        raise OSError(f'{path_to_db_bkgd} default file not found. Download a default background using python -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"')
      else:
        bar = IncrementalBar(f'Loading {self.default_bkgd}', max = 2, suffix='%(percent)d%%')
        bar.next()
        if path_to_db_bkgd.endswith(".csv"):
          sep = ","
        elif path_to_db_bkgd.endswith(".tsv"):
          sep = "\t"
        self.ref_df= pd.read_csv(path_to_db_bkgd, sep = sep)
        bar.next()
        bar.finish()
        self.build_background()
  
  @classmethod
  def currently_available_backgrounds(self):
    """ 
    List available default background already in the /db/ folder.

    These file are immediately available for use as TCRsampler(default_background = 'file')

    Returns
    -------
    available_tsv_and_csv : list
    """
    path_to_db = os.path.join(os.path.dirname(os.path.realpath(__file__)),'db')
    available_tsv_and_csv = [f for f in os.listdir(path_to_db) if f.endswith('sv') and f not in ['alphabeta_db.tsv', 'gammadelta_db.tsv']]
    for filename in sorted(available_tsv_and_csv):
      sys.stdout.write(f"{filename}\n")
    available_tsv_and_csv =  sorted(available_tsv_and_csv)
    return available_tsv_and_csv

  @classmethod
  def download_all_background_files(self, dry_run = False):
    """
    Downloads, unzips, and installs all the priority default background files to package source tcrsampler/db/ folder. 

    Parameters
    ----------
    dry_run : bool 
      If True, only prints curl strings 

    Notes
    -----
    This function will only work on OSX/Linus systems.
    """
    install_all_next_gen(dry_run = dry_run)
  
  @classmethod
  def download_background_file(self, download_file, dry_run = False, download_from = "dropbox"):
    """
    Downloads, unzips, and installs a specific background files. 

    Parameters
    ----------
    download_file : str
      string of file to download (e.g. "britanova_human_beta_t_cb.tsv.sampler.tsv.zip", "emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip" , 
      "wiraninha_sampler.zip", "ruggiero_mouse_sampler.zip" , "ruggiero_human_sampler.zip" )
    dry_run : bool 
      If True, only prints curl strings 

    Notes
    -----
    Downloads, unzips, and installs file directly to package source tcrsampler/db/ folder. 

    These .zip file contain multiple backgound files. 

    This function will only work on OSX/Linus systems
    """
    install_nextgen_data_to_db(download_file = download_file, download_from = download_from, dry_run = dry_run)

  @classmethod
  def get_download_address(self):
    """
    Provides list of availble .zip files that can be downloaded individually. 
    """
    for file in select_files:
      sys.stdout.write(f"\t{file}\n") 
      sys.stdout.write(f"Download file {file} with TCRsampler.download_background_file(download_file = '{file}')\n")
      sys.stdout.write(f"\nDownload all background files with TCRsampler.download_all_background_files()\n")
    return select_files


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

  def clean_mixcr(self, filename = None, df = None):
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
    if df is None:
      df = pd.read_csv(filename, sep = "\t")
    
    bar = IncrementalBar('Clean Mixcr     ', max = 1, suffix='%(percent)d%%')
    
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

  def build_background( self, 
                        df = None, 
                        max_rows = 100, 
                        stratify_by_subject = False, 
                        use_frequency= True, 
                        make_singleton = False):
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

    #bar = IncrementalBar('Build Background', max = dfg.ngroups, suffix='%(percent)d%%')
    
    bar = IncrementalBar('V-J Frequency Dict ', max = 3, suffix='%(percent)d%%')
    # V_J OCCURRENCE PROBABILITIES BY THE SEQUENCE FREQUENCY METHOD
    vj = dfg['freq'].sum().reset_index()

    # Divide frequency by sum of total frequncy across all subjects, assumes most subject frequency columns sum roughly to 1
    bar.next()
    vj = vj.assign(freq = np.divide(vj.freq,vj.freq.sum()))

    # Assign a tuple key_dictionary to < .vj_freq >
    self.vj_freq = {(v,j):f for v,j,f in vj.to_dict('split')['data']}
    bar.next()
    # Calculate the marginal V-gene frequencies
    dfg = df.groupby(['v_reps'])
    v_only = dfg['freq'].sum().reset_index()
    v_only = v_only.assign(freq = np.divide(v_only.freq, v_only.freq.sum()))
    self.v_freq = {v:f for v,f in v_only.to_dict('split')['data']}
          # TEST THAT assert np.isclose(np.sum([k for _,k in t.v_freq.items()]), 1.0)
    # Calculate the marginal J-gene frequencies
    dfg = df.groupby(['j_reps'])
    j_only = dfg['freq'].sum().reset_index()
    j_only = j_only.assign(freq = np.divide(j_only.freq, j_only.freq.sum()))
    self.j_freq = {j:f for j,f in j_only.to_dict('split')['data']}
          
    bar.next();bar.finish()

    # V_J OCCURRENCE PROBABILITIES BY THE UNIQUE N CLONES METHOD (see NOTES)
    bar = IncrementalBar('V-J Occurrence Dict', max = 3, suffix='%(percent)d%%')
    df_gb_sub = df.\
      sort_values(['freq','subject'], ascending = False).\
      groupby(['subject'])
      # Count # of clone per subject.
    nclone_per_subject = [group.shape[0] for _,group in df_gb_sub]
      # Set N to be the number of clones in the subject with the fewest clones
    N = np.min(nclone_per_subject )
      # Take only the top N from each subject, and count how often these occured
      # (note thate we ignore the frequency which reflects clonal expansion)
    bar.next()
    vj_occur = df_gb_sub.head(N).groupby(['v_reps', 'j_reps']).\
      size().\
      reset_index().\
      rename(columns = {0:'occur'})
    vj_occur = vj_occur.assign(occur = np.divide(vj_occur.occur, vj_occur.occur.sum()))
    assert np.isclose(vj_occur.occur.sum(), 1.0, rtol = 0.000001)
    self.vj_occur_freq = {(v,j):f for v,j,f in vj_occur.to_dict('split')['data']}

    bar.next()
    # Calculate the marginal V-gene frequencies
    v_occur = df_gb_sub.head(N).groupby(['v_reps']).\
      size().\
      reset_index().\
      rename(columns = {0:'occur'})
    v_occur = v_occur.assign(occur = np.divide(v_occur.occur, v_occur.occur.sum()))
    assert np.isclose(vj_occur.occur.sum(), 1.0, rtol = 0.000001)
    self.v_occur_freq = {v:f for v,f in v_occur.to_dict('split')['data']}

    # Calculate the marginal J-gene frequencies
    j_occur = df_gb_sub.head(N).groupby(['j_reps']).\
      size().\
      reset_index().\
      rename(columns = {0:'occur'})
    j_occur = j_occur.assign(occur = np.divide(j_occur.occur, j_occur.occur.sum()))
    assert np.isclose(j_occur.occur.sum(), 1.0, rtol = 0.000001)
    self.j_occur_freq = {j:f for j,f in j_occur.to_dict('split')['data']}

    bar.next();bar.finish()


    dfg =df.groupby(['v_reps','j_reps'])
    bar = IncrementalBar('Build Background   ', max = dfg.ngroups, suffix='%(percent)d%%')
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

