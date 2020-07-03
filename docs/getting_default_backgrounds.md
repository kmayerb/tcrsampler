
# Default Background Files

tcrsampler provides background data files for human and mouse alpha and beta TCR chains. 

Once installed these files reside in the `/tcrsampler/db/` folder of the package source. 


## Download All Backgrounds at Time of Install

While the number of background files is relatively small. The easist way to get all background data files is with a single commandline command.

```
python -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"
```

## Download Backgrounds As Needed

However, in the future, we envision that some users may wish to extert more control over which background files are installed or choose to install a 
background file ants needed in an ieractive session. 

```python
In [1]: from tcrsampler.sampler import TCRsampler

In [2]: TCRsampler.currently_available_backgrounds()
Out[2]: []

In [3]: TCRsampler.get_download_address()
```

```
	wiraninha_sampler.zip
Download file wiraninha_sampler.zip with TCRsampler.download_background_file(download_file = 'wiraninha_sampler.zip')

Download all background files with TCRsampler.download_all_background_files()
	ruggiero_mouse_sampler.zip
Download file ruggiero_mouse_sampler.zip with TCRsampler.download_background_file(download_file = 'ruggiero_mouse_sampler.zip')

Download all background files with TCRsampler.download_all_background_files()
	ruggiero_human_sampler.zip
Download file ruggiero_human_sampler.zip with TCRsampler.download_background_file(download_file = 'ruggiero_human_sampler.zip')

Download all background files with TCRsampler.download_all_background_files()
	emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip
Download file emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip with TCRsampler.download_background_file(download_file = 'emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip')

Download all background files with TCRsampler.download_all_background_files()
	britanova_human_beta_t_cb.tsv.sampler.tsv.zip
Download file britanova_human_beta_t_cb.tsv.sampler.tsv.zip with TCRsampler.download_background_file(download_file = 'britanova_human_beta_t_cb.tsv.sampler.tsv.zip')

Download all background files with TCRsampler.download_all_background_files()
```

```python
In [4]: TCRsampler.download_background_file(download_file = 'wiraninha_sampler.zip')
```

```
RUNNING: curl -o /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wiraninha_sampler.zip https://www.dropbox.com/s/ily0td3tn1uc7bi/wiraninha_sampler.zip?dl=1 -L
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
100 6882k  100 6882k    0     0  3702k      0  0:00:01  0:00:01 --:--:-- 13.0M
Archive:  /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wiraninha_sampler.zip
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_g8a.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_s_4.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_s_8.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_s_r.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_t_1.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_t_4.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_t_8.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_t_p.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_alpha_t_r.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_g8a.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_s_4.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_s_8.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_s_r.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_t_1.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_t_4.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_t_8.tsv.sampler.tsv
  inflating: /Users/kmayerbl/TCRDIST/tcrsampler/tcrsampler/db/wirasinha_mouse_beta_t_p.tsv.sampler.tsv
```

Now the backgrounds are available.

```In [5]: TCRsampler.currently_available_backgrounds()
Out[5]:
['wirasinha_mouse_alpha_g8a.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_s_4.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_s_8.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_s_r.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_t_1.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_t_4.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_t_8.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_t_p.tsv.sampler.tsv',
 'wirasinha_mouse_alpha_t_r.tsv.sampler.tsv',
 'wirasinha_mouse_beta_g8a.tsv.sampler.tsv',
 'wirasinha_mouse_beta_s_4.tsv.sampler.tsv',
 'wirasinha_mouse_beta_s_8.tsv.sampler.tsv',
 'wirasinha_mouse_beta_s_r.tsv.sampler.tsv',
 'wirasinha_mouse_beta_t_1.tsv.sampler.tsv',
 'wirasinha_mouse_beta_t_4.tsv.sampler.tsv',
 'wirasinha_mouse_beta_t_8.tsv.sampler.tsv',
 'wirasinha_mouse_beta_t_p.tsv.sampler.tsv',
 'wirasinha_mouse_beta_t_r.tsv.sampler.tsv']
 ```

For instance, 

```pyhton
In [6]: t = TCRsampler(default_background = 'wirasinha_mouse_alpha_t_4.tsv.sampler.tsv')
```
```
Loading wirasinha_mouse_alpha_t_4.tsv.sampler.tsv |████████████████████████████████| 100%
V-J Frequency Dict  |████████████████████████████████| 100%
V-J Occurrence Dict |████████████████████████████████| 100%
Build Background    |████████████████████████████████| 100%
```

Or, for instance,

```python 
In [7]: TCRsampler.download_background_file(download_file = 'britanova_human_beta_t_cb.tsv.sampler.tsv.zip')
In [8]: t = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
```

```
Loading britanova_human_beta_t_cb.tsv.sampler.tsv |████████████████████████████████| 100%
V-J Frequency Dict  |████████████████████████████████| 100%
V-J Occurrence Dict |████████████████████████████████| 100%
Build Background    |████████████████████████████████| 100%
```

