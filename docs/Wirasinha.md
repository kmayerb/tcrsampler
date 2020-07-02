



## Data Source

Wirasinha RC, Singh M, Archer SK, et al. αβ T-cell receptors with a central CDR3 cysteine are enriched in CD8αα intraepithelial lymphocytes and their thymic precursors.* 
Immunol Cell Biol. 2018;96(6):553-561. doi:10.1111/imcb.12047

[Data](https://figshare.com/articles/Wirasinha_migec_txt_zip/5766591)

[Description](https://researchdata.edu.au/wirasinhamigectxtzip/1307755)

The file contains mouse T cell receptor (TCR) sequences collected by multiplex PCR amplification of 
cDNA molecules followed by Illumina sequencing. Sequences were aligned to the mouse genome using MIGEC software 
(see doi: 10.1038/nmeth.2960 for details). 

Except for the header row, each row contains information about a unique TCR nucleotide sequence. 

* Columns 1-11 contain output from MIGEC software. 
* Columns 12-16 contain information about sample origin, detailed as follows: 
  * column 12 "sample_id" is an identifier for the sample of origin; 
  * column 13 "tcr_b" specifies the TCR beta chain status of donor mouse ("poly" = C57BL/6 inbred mouse strain); 
  * column 14 "organ" specifies organ of origin (thymus, "gut" = small intestine and spleen); 
  * column 15 "subset" specifies the T cell subset of origin 
    * "t_p" = pre-selection thymocytes; 
    * "t_1" = wave 1 thymocyte
    * "t_r" = thymic T-reg
    * "t_4" = CD4+ naive thymocytes
    * "t_8" = CD8+ naive thymocytes; 
    * "g8a" = small intestinal CD8aa+ intraepithelial lymphocytes; 
    * "s_r" = splenic T-reg; 
    * "s_4" = splenic CD4+ naive T cells and 
    * "s_8" = splenic CD8+ naive T cells; 
  * column 16 is identical to column 15 except that column 16 shows whether the sample of origin was a CD25+ or CD25– fraction of a thymic or splenic T-reg population.


"poly" = C57BL/6 inbred mouse strain: C57BL/6J mice are resistant to audiogenic seizures, have a relatively low bone density, and develop age related hearing loss. They are also susceptible to diet-induced obesity, type 2 diabetes, and atherosclerosis. Macrophages from this strain are resistant to the effects of anthrax lethal toxin.

The YAe-62.8 (YAe62) TCR, from a single peptide IAb mouse, has this ability to recognize both MHCI and various alleles of MHCII. We recently reported the structure of the YAe62 TCR bound to the MHCII molecule, IAb, plus a known peptide, p3K (Dai et al., 2008). However the nature of YAe62’s cross reactivity with class I was unknown.

Anti-3K peptide T cell receptor (B3K506)
The vector of anti-3K peptide T cell receptor (TCR) is constructed for the engineering of T cell to target 3K peptide. The T cells are genetically modified through transduction with a lentiviral vector expressing 3K peptide-specific T cell receptor.


``` python

import pandas as pd
from tcrsampler.sampler import TCRsampler

fn = '/Volumes/Samsung_T5/kmayerbl/tcr_data/wirasinha/Wirasinha.migec.txt'
df = pd.read_csv(fn, sep= '\t')

def wirasinha_subset(df, subset = 't_p', tcr_b = 'poly', chain = 'b'):
 ind1 = df['subset'] == subset
 ind2 = df['tcr_b'] == tcr_b
 ind3 = df['chain'] == chain
 sdf = df[(ind1&ind2&ind3)]
 return sdf

to_mixcr_headers = \
				{'count'        :'cloneCount',
				'freq'          :'cloneFraction', 
				'cdr3nt'        :'nSeqCDR3',
				'cdr3aa'        :'aaSeqCDR3',
				'v'             :'bestVGene',
				'j'             :'bestJGene',
				'd'             :'bestDGene',
				'VEnd'          :'VEnd',
				'DStart'        :'DStart' ,
				'DEnd'          :'DEnd',
				'JStart'        :'JStart',
				'sample_id'     :'subject',
				'tcr_b'         :'tcr_b',
				'chain'         :'chain',
				'mouse'         :'mouse',
				'organ'         :'organ',
				'subset'        :'subset',
				'subset_primary':'subest'}

import re

def _strip_allele(s):
 try:
  gr = re.match(string = s , pattern ='(.*)\*')
  return gr[1]
 except:
  return s
  
def strip(df):
	df['bestVGene'] = df['bestVGene'].apply(lambda x: _strip_allele(x))
	df['bestJGene'] = df['bestJGene'].apply(lambda x: _strip_allele(x))
	return df
 
wirasinha_mouse_beta_tp_c57bl6  = wirasinha_subset(df = df, subset = 't_p', tcr_b = 'poly', chain = 'b').copy().rename(columns = to_mixcr_headers)
wirasinha_mouse_beta_tp_c57bl6['bestVGene'] = wirasinha_mouse_beta_tp_c57bl6['bestVGene'].apply(lambda x: _strip_allele(x))
wirasinha_mouse_beta_tp_c57bl6['bestJGene'] = wirasinha_mouse_beta_tp_c57bl6['bestJGene'].apply(lambda x: _strip_allele(x))

wirasinha_mouse_alpha_tp_c57bl6 = wirasinha_subset(df = df, subset = 't_p', tcr_b = 'poly', chain = 'a').copy().rename(columns = to_mixcr_headers)
wirasinha_mouse_alpha_tp_c57bl6['bestVGene'] = wirasinha_mouse_alpha_tp_c57bl6['bestVGene'].apply(lambda x: _strip_allele(x))
wirasinha_mouse_alpha_tp_c57bl6['bestJGene'] = wirasinha_mouse_alpha_tp_c57bl6['bestJGene'].apply(lambda x: _strip_allele(x))


wirasinha_mouse_beta_t4_c57bl6  = wirasinha_subset(df = df, subset = 't_4', tcr_b = 'poly', chain = 'b').copy().rename(columns = to_mixcr_headers)
wirasinha_mouse_beta_t4_c57bl6['bestVGene'] = wirasinha_mouse_beta_t4_c57bl6['bestVGene'].apply(lambda x: _strip_allele(x))
wirasinha_mouse_beta_t4_c57bl6['bestJGene'] = wirasinha_mouse_beta_t4_c57bl6['bestJGene'].apply(lambda x: _strip_allele(x))

wirasinha_mouse_alpha_t4_c57bl6 = wirasinha_subset(df = df, subset = 't_4', tcr_b = 'poly', chain = 'a').copy().rename(columns = to_mixcr_headers)
wirasinha_mouse_alpha_t4_c57bl6['bestVGene'] = wirasinha_mouse_alpha_t4_c57bl6['bestVGene'].apply(lambda x: _strip_allele(x))
wirasinha_mouse_alpha_t4_c57bl6['bestJGene'] = wirasinha_mouse_alpha_t4_c57bl6['bestJGene'].apply(lambda x: _strip_allele(x))

wirasinha_mouse_beta_t8_c57bl6  = wirasinha_subset(df = df, subset = 't_8', tcr_b = 'poly', chain = 'b').copy().rename(columns = to_mixcr_headers)
wirasinha_mouse_beta_t8_c57bl6['bestVGene'] = wirasinha_mouse_beta_t8_c57bl6['bestVGene'].apply(lambda x: _strip_allele(x))
wirasinha_mouse_beta_t8_c57bl6['bestJGene'] = wirasinha_mouse_beta_t8_c57bl6['bestJGene'].apply(lambda x: _strip_allele(x))

wirasinha_mouse_alpha_t8_c57bl6 = wirasinha_subset(df = df, subset = 't_8', tcr_b = 'poly', chain = 'a').copy().rename(columns = to_mixcr_headers)
wirasinha_mouse_alpha_t8_c57bl6['bestVGene'] = wirasinha_mouse_alpha_t8_c57bl6['bestVGene'].apply(lambda x: _strip_allele(x))
wirasinha_mouse_alpha_t8_c57bl6['bestJGene'] = wirasinha_mouse_alpha_t8_c57bl6['bestJGene'].apply(lambda x: _strip_allele(x))

backgrounds = \
	{'wirasinha_mouse_beta_tp_c57bl6' :wirasinha_mouse_beta_tp_c57bl6,
	'wirasinha_mouse_alpha_tp_c57bl6' :wirasinha_mouse_alpha_tp_c57bl6, 
	'wirasinha_mouse_beta_t4_c57bl6'  :wirasinha_mouse_beta_t4_c57bl6, 
	'wirasinha_mouse_alpha_t4_c57bl6' :wirasinha_mouse_alpha_t4_c57bl6, 
	'wirasinha_mouse_beta_t8_c57bl6'  :wirasinha_mouse_beta_t8_c57bl6,
	'wirasinha_mouse_alpha_t8_c57bl6' :wirasinha_mouse_alpha_t8_c57bl6}

for k,bkgd in backgrounds.items():
 print(f"Testing {k} for import into TCRsampler")
 t = TCRsampler()
 t.clean_mixcr(df = bkgd)
 t.build_background()

for k,bkgd in backgrounds.items():
 bkgd.to_csv(f'{k}.tsv', sep = '\t')
 
```

