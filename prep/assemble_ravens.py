"""
assemble_ravens.py (formerly 03_prep_ravens.py)

Septemeber 1, 2020 
kmayerb


Produces
--------

"ravens_human_delta_t.sampler.tsv"

"ravens_human_gamma_t.sampler.tsv"


Process
-------

This file describes the steps used to make a human gamma and delta background set. 

1. The following files were downloaded to from SRA using the script (https://github.com/kmayerb/kutils/blob/master/sra_to_fastq.py)
2. Then the .fasta file were run through mixcr using (/Volumes/Samsung_T5/kmayerbl/gd/SRA/01_batch_mixcr.py)  w/in milaboratory/mixcr:3-imgt
3. And the /Volumes/Samsung_T5/kmayerbl/gd/SRA/03_prep_ravens.py


Files Include:

gs = ['SRR5130201.1.clns.txt', 'SRR5130216.1.clns.txt','SRR5130220.1.clns.txt', 'SRR5130234.1.clns.txt', 'SRR5130282.1.clns.txt', 'SRR5130297.1.clns.txt','SRR5130299.1.clns.txt', 'SRR5130301.1.clns.txt']
gs = ['SRR5130231.1.clns.txt','SRR5130249.1.clns.txt','SRR5130255.1.clns.txt','SRR5130260.1.clns.txt','SRR5130264.1.clns.txt','SRR5130267.1.clns.txt','SRR5130268.1.clns.txt','SRR5130285.1.clns.txt']

Notes
-----

Ravens https://www.nature.com/articles/ni.3686

Published: 20 February 2017
Human γδ T cells are quickly reconstituted after stem-cell transplantation and show adaptive clonal expansion in response to viral infection

Sarina Ravens, Christian Schultze-Florey, Solaiman Raha, Inga Sandrock, Melanie Drenker, Linda Oberdörfer, Annika Reinhardt, Inga Ravens, Maleen Beck, Robert Geffers, Constantin von Kaisenberg, Michael Heuser, Felicitas Thol, Arnold Ganser, Reinhold Förster, Christian Koenecke & Immo Prinz 
Nature Immunology volume 18, pages393–401(2017)Cite this article

"""



import pandas as pd
import math
from tcrdist.adpt_funcs import _valid_cdr3

df = pd.read_csv("/Volumes/Samsung_T5/kmayerbl/gd/SRA/output/all_files.clns.tsv", sep = "\t")
df = df[['allVHitsWithScore',  'allJHitsWithScore', 'aaSeqCDR3', 'nSeqCDR3', 'cloneCount', 'cloneFraction','source']]

cache  = []
for i,row in df.iterrows():
	nsplit = len(row['allJHitsWithScore'].split(",")) * len(row['allJHitsWithScore'].split(","))
	for j in row['allJHitsWithScore'].split(","):
		for v in row['allVHitsWithScore'].split(","):
			new_row = row.copy()
			new_row['allJHitsWithScore'] = j
			new_row['allVHitsWithScore'] = v
			new_row['cloneFraction'] = new_row['cloneFraction'] / nsplit 
			new_row['cloneCount'] = math.ceil(new_row['cloneCount'] / nsplit)
			cache.append(new_row.copy())


dfnew = pd.concat(cache, axis = 1)
dfnew = dfnew.transpose()

# patche one sets all to *01 allele and remove score
patch = lambda s : s.split('(')[0].replace("00", "01")

swaps = {'TRAV14DV4*01':  'TRAV14/DV4*01',
 'TRAV14DV4*02':  'TRAV14/DV4*02',
 'TRAV14DV4*03':  'TRAV14/DV4*03',
 'TRAV14DV4*04':  'TRAV14/DV4*04',
 'TRAV23DV6*01':  'TRAV23/DV6*01',
 'TRAV23DV6*02':  'TRAV23/DV6*02',
 'TRAV23DV6*03':  'TRAV23/DV6*03',
 'TRAV23DV6*04':  'TRAV23/DV6*04',
 'TRAV29DV5*01':  'TRAV29/DV5*01',
 'TRAV29DV5*02':  'TRAV29/DV5*02',
 'TRAV36DV7*01':  'TRAV36/DV7*01',
 'TRAV36DV7*02':  'TRAV36/DV7*02',
 'TRAV36DV7*03':  'TRAV36/DV7*03',
 'TRAV36DV7*04':  'TRAV36/DV7*04',
 'TRAV38-2DV8*01':'TRAV38-2/DV8*01'}

def swap(x, swaps):
	if swaps.get(x) is not None:
		return swaps.get(x)
	else:
		return x

# patch2 swaps possible misnamed genes
patch2 = lambda x : swap(x,swaps=swaps)

dfnew['v'] = df.allVHitsWithScore.apply(patch).apply(patch2)
dfnew['j'] = df.allJHitsWithScore.apply(patch).apply(patch2)


# rename columns
new_cols = {'v':'v_reps', 'j':'j_reps', 'aaSeqCDR3':'cdr3', 'cloneCount':'count','cloneFraction':'freq','source':'subject'}
dfnew = dfnew[['v','j', 'aaSeqCDR3', 'cloneCount','cloneFraction','source']].rename(columns = new_cols).reset_index(drop = True)


dfnew = dfnew[dfnew.cdr3.apply(lambda x : _valid_cdr3(x))]
dfnew = dfnew[dfnew.cdr3.apply(lambda x : len(x) > 5)]

dfg = dfnew[(dfnew.v_reps.str.startswith('TRG') | dfnew.j_reps.str.startswith('TRG') ) ]
dfd = dfnew[(dfnew.v_reps.str.startswith('TRD') | dfnew.j_reps.str.startswith('TRD') ) ]


# Only include correct files
gs = ['SRR5130201.1.clns.txt', 'SRR5130216.1.clns.txt','SRR5130220.1.clns.txt', 'SRR5130234.1.clns.txt', 'SRR5130282.1.clns.txt', 'SRR5130297.1.clns.txt','SRR5130299.1.clns.txt', 'SRR5130301.1.clns.txt']
ds = ['SRR5130231.1.clns.txt','SRR5130249.1.clns.txt','SRR5130255.1.clns.txt','SRR5130260.1.clns.txt','SRR5130264.1.clns.txt','SRR5130267.1.clns.txt','SRR5130268.1.clns.txt','SRR5130285.1.clns.txt']
dfg = dfg[dfg.subject.apply(lambda x: x in gs)]
dfd = dfd[dfd.subject.apply(lambda x: x in ds)]

# Correct fequency to account for the removal and duplications
dfg = dfg.merge(dfg[['freq','subject']].groupby(['subject']).sum().reset_index(), how ='left', on = 'subject')
dfg['freq'] = dfg['freq_x']/ dfg['freq_y']
print(dfg[['freq','subject']].groupby(['subject']).sum())

dfd = dfd.merge(dfd[['freq','subject']].groupby(['subject']).sum().reset_index(), how ='left', on = 'subject')
dfd['freq'] = dfd['freq_x']/ dfd['freq_y']
print(dfd[['freq','subject']].groupby(['subject']).sum())

# Test that these will work with TCRsampler

from tcrsampler.sampler import TCRsampler

from tcrdist import repertoire_db
ref = repertoire_db.RefGeneSet(db_file = 'alphabeta_gammadelta_db.tsv')
ref.generate_all_genes()
ref.all_genes
ref.all_genes['human'].keys()


tsd = TCRsampler()
tsd.ref_df = dfd
tsd.build_background()
# find potential missing:
print([x for x in tsd.v_freq.keys()])
print([x for x in tsd.v_freq.keys() if x not in ref.all_genes['human'].keys()])
assert len([x for x in tsd.v_freq.keys() if x not in ref.all_genes['human'].keys()]) == 0
print([x for x in tsd.j_freq.keys()])
print([x for x in tsd.j_freq.keys() if x not in ref.all_genes['human'].keys()])
assert len([x for x in tsd.j_freq.keys() if x not in ref.all_genes['human'].keys()]) == 0


tsg = TCRsampler()
tsg.ref_df = dfg
tsg.build_background()
print([x for x in tsg.v_freq.keys()])
print([x for x in tsg.v_freq.keys() if x not in ref.all_genes['human'].keys()])
assert len([x for x in tsg.v_freq.keys() if x not in ref.all_genes['human'].keys()]) == 0
print([x for x in tsg.j_freq.keys()])
print([x for x in tsg.j_freq.keys() if x not in ref.all_genes['human'].keys()])
assert len([x for x in tsg.j_freq.keys() if x not in ref.all_genes['human'].keys()]) == 0

dfd[['v_reps','j_reps','cdr3','count','freq','subject']].to_csv("ravens_human_delta_t.sampler.tsv",index = False)
dfg[['v_reps','j_reps','cdr3','count','freq','subject']].to_csv("ravens_human_gamma_t.sampler.tsv",index = False)

