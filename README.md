# tcrsampler

Sample TCRs according to VDJ gene usage frequency

[![Build Status](https://travis-ci.com/kmayerb/tcrsampler.svg?branch=master)](https://travis-ci.com/kmayerb/tcrsampler)
[![Coverage Status](https://coveralls.io/repos/github/kmayerb/tcrsampler/badge.svg?branch=master)](https://coveralls.io/github/kmayerb/tcrsampler?branch=master)
[![PyPI version](https://badge.fury.io/py/tcrsampler.svg)](https://badge.fury.io/py/tcrsampler)

## Install 

```
pip install tcrsampler
```

### Example

```python
import os
import pandas as pd
from tcrsampler.sampler import TCRsampler
t = TCRsampler()
fn= os.path.join('tcrsampler' ,'tests', 'pmbc_mixcr_example_data.txt')
t.clean_mixcr(filename = fn)
t.build_background(stratify_by_subject = True, use_frequency = True)
t.sample([['TRBV9*01','TRBJ2-7*01', 2],['TRBV7-7*01', 'TRBJ2-4*01', 4] ], depth = 1, use_frequency = True)
```

```
[['CASSRTGSLADEQYF', 'CASSATGVVSAQYF'],
 ['CASSLGQAARGIQYF', 'CASSLGQAARGIQYF', 'CASSLGQAARGIQYF', 'CASSLGQAARGIQYF']]
```

## Download Default Background

```python
python -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"
```
