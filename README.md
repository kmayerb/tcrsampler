# tcrsampler

Sample TCRs according to VDJ gene usage frequency

[![Build Status](https://travis-ci.com/kmayerb/tcrsampler.svg?branch=master)](https://travis-ci.com/kmayerb/tcrsampler)

## Install 

```
pip install tcrsampler
```

Install default reference files:

```
python -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"
```

### Example

```python
from tcrsampler.sampler import TCRsampler
t = TCRsampler(organism = 'human',db_file = 'alphabeta_db.tsv', ref_file= 'new_nextgen_chains_human_B.tsv')
t.make_ref_dict(max_sample = 100)
t.sample_ref_dict('TRBV1*01', 'TRBJ1-1*01',1)
t.sample([['TRBV1*01', 'TRBJ1-1*01',1], ['TRBV2*01', 'TRBJ1-1*01',1]])
```