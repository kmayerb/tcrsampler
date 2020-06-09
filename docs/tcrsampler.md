# tcrsampler

Sample TCRs according to VDJ gene usage frequency

[![Build Status](https://travis-ci.com/kmayerb/tcrsampler.svg?branch=master)](https://travis-ci.com/kmayerb/tcrsampler)
[![Coverage Status](https://coveralls.io/repos/github/kmayerb/tcrsampler/badge.svg?branch=master)](https://coveralls.io/github/kmayerb/tcrsampler?branch=master)
[![PyPI version](https://badge.fury.io/py/tcrsampler.svg)](https://badge.fury.io/py/tcrsampler)

The current version of tcrsamper is a pre-release and is subject to change.


## Install 

```
pip install tcrsampler
```

## Get Background Reference Data

The files are:

* **britanova_chord_blood.csv** (80 MB) - It contains beta-chain repertoires from 8 umbilical cord samples [download](https://www.dropbox.com/s/rkbce72njcei4y8/britanova_chord_blood.csv?dl=1)

* **emerson_cmv_negative.csv** (120 MB) - It contains beta-chain repertoires from 7 CMV negative individuals aged 18-21 [download](https://www.dropbox.com/s/yrozbowxtqumjfl/emerson_cmv_negative.csv?dl=1).

If you use these files in your research, cite the appropriate primary studies:

* Britanova OV, Shugay M, Merzlyak EM, Staroverov DB, Putintseva EV, Turchaninova MA, Mamedov IZ, Pogorelyy MV, Bolotin DA, Izraelson M, et al. Dynamics of individual T cell repertoires: from cord blood to centenarians. J Immunol. 2016;196:5005–5013. doi: 10.4049/jimmunol.1600005.

* Emerson RO, DeWitt WS, Vignali M, Gravley J, Hu JK, Osborne EJ, Desmarais C, Klinger M, Carlson CS, et al. Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nat Genet. 2017;49:659–665. doi: 10.1038/ng.3822. 

Running the optional command will install reference data sets to the `tcrsampler/db/` folder within your Python environment (e.g., `python3.7/site-packages/tcrampler/db`)

```
python -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"
```

## Build the Reference Sampler

Building a tcrsampler is straightforward.

```python
import os
import pandas as pd
from tcrsampler.sampler import TCRsampler
fn = 'britanova_chord_blood.csv' 
t = TCRsampler()
t.ref_df = pd.read_csv(fn)
t.build_background()
```

After `.build_backgroud` the TCRsampler instnace is ready to draw CDR3s from the background. CDRs are returned to match 
a given v and j pairing. By default, the likelihood of drawing a given CDR3 is proportional to its frequency in 
the samples used to build the background.

```python
t.sample_background(v ='TRBV10-1*01', j ='TRBJ1-1*01',n=1, depth = 1, seed =1, use_frequency= True ) 
```

Sample an entire repertoire using `.sample`. This example is a trivial repertoire with only 3 representatives, but it illustrates
the I/O of the `.sample` method.

```python
t.sample([['TRBV10-1*01', 'TRBJ1-1*01', 1],['TRBV9*01','TRBJ2-7*01', 2]]) 
```
```
[['CASSDSMNTEAFF'], ['CASSVGTNSYEQYF', 'CASSAGTNYEQYF']]
```


### Build directly from a MiXCR file

MiXCR is free for academic and non-profit use (see [License](https://mixcr.readthedocs.io/en/master/license.html#license)). MiXCR is a universal framework that processes big immunome data from raw sequences to quantitated clonotypes. 

tcrsampler can build a reference sampler directly from a MiXCR output. 

This method enables the definition of custom background, which could be subject or HLA-type specific.

```python
import os
import pandas as pd
from tcrsampler.sampler import TCRsampler
t = TCRsampler()
fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
t.clean_mixcr(filename = fn)
t.build_background()
```


### Build Options


The first step is to build a background dictionary of CDR3s from repertoire data.

During the build, clones are sorted in descending order from most frequent to least for each v,j gene pairing.

`.sample` parameters control how many CDRs are included in the build.

* `max_rows` (int) specifies the maximum clones per v-gene, j-gene pair that will be included in the build.
* `stratify_by_subject` (bool) Set to True by default, max_rows will apply for each subject. 
For example, if max_rows = 100, the 100 most frequent clones for each v,j gene usage pair will be included from each subject.
If set to False, the number of max_rows will apply to the most frequent clone from any subject. 
* `use_frequency` (bool) Set to True by default, uses frequency for ranking rows not counts.
* `make_singleton` (bool), Set to False by default, however, if True background is still sorted by frequency, but the final frequency and counts values are overridden and set to 1. As a consequence, every included CDR3 has an equal likelihood of being sampled at the sampling phase.


## The Sampler

The goal of tcrsampler is to allow representative CDR3s to be sampled from a background set. By default, the likelihood that a CDR3 is drawn from the background is proportional to the frequency of that
CDR3 clone in the original sample. That is, more abundant clones will be sampled more frequently than less abundant clones. 

* `depth` how many CDR3 to draw per input v,j combination. 
* `seed` random seed effects the stochastic nature of the sampling.
* `use_frequency` set to True by default, uses frequency for determining the probability of picking a CDR3. If False, use raw counts.
* `flatten` Set to False by default, however, to true a single list is returned.


Sampling:

```python
t.sample([['TRBV10-1*01', 'TRBJ1-1*01', 1],['TRBV9*01','TRBJ2-7*01', 2]], depth = 1)
```
```
[['CASSDSMNTEAFF'], ['CASSVGTNSYEQYF', 'CASSAGTNYEQYF']]
```

```python
t.sample([['TRBV10-1*01', 'TRBJ1-1*01', 1],['TRBV9*01','TRBJ2-7*01', 2]], depth = 2)
```

```
[['CASSDSMNTEAFF', 'CASSAGRTEAFF'],
 ['CASSVGTNSYEQYF', 'CASSAGTNYEQYF', 'CASSVGTGGYEQYF', 'CASSGTTYEQYF']]
```

```python
t.sample([['TRBV10-1*01', 'TRBJ1-1*01', 1],['TRBV9*01','TRBJ2-7*01', 2]], depth = 3)
```

```
[['CASSDSMNTEAFF', 'CASSAGRTEAFF', 'CASSEIPGSVMNTEAFF'],
 ['CASSVGTNSYEQYF','CASSAGTNYEQYF', 'CASSVGTGGYEQYF', 'CASSGTTYEQYF', 'CASSDSYEQYF', 'CASSVDSSYEQYF']]
```

Sampling a reperoire repressented as a Pandas DataFrame is also possible. Let's consider the following example:

```python
df = pd.DataFrame( { "v_gene":['TRBV9*01','TRBV9*01','TRBV9*01'],
					"j_gene":['TRBJ2-6*01','TRBJ2-6*01','TRBJ2-5*01']})
gene_usage = df.groupby(['v_gene','j_gene']).size()
sampled_rep = t.sample( gene_usage.reset_index().to_dict('split')['data'], flatten = True)
```

```
['CASSVGLKETQYF', 'CASSVQHSGANVLTF', 'CASSLDPQGTGANVLTF']
```

Alternatively, a sampled CDR3 can be added to a column of an existing DataFrame.

```python
import random
df = pd.DataFrame( { "v_gene":['TRBV9*01','TRBV9*01','TRBV9*01'],
					 "j_gene":['TRBJ2-6*01','TRBJ2-6*01','TRBJ2-5*01']})
df['sampledcdr3'] = df.apply(lambda x : \
				t.sample_background(	v = x['v_gene'],\
							j = x['j_gene'],\
							n = 1,\
							seed = random.randint(1,1000))[0],\
							axis = 1)
```

```
     v_gene      j_gene          sampledcdr3
0  TRBV9*01  TRBJ2-6*01  CASSVAQAGTHSGANVLTF
1  TRBV9*01  TRBJ2-6*01   CASSLRPGISSGANVLTF
2  TRBV9*01  TRBJ2-5*01      CASSTSIGAQETQYF
```


