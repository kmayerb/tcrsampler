"""
July 2, 2020

This script cleans wirasinha/Wirasinha.migec.txt and produces:

 15M Jul  2 11:46 wirasinha_mouse_beta_t_p.tsv
5.5M Jul  2 11:46 wirasinha_mouse_beta_t_p.tsv.sampler.tsv
8.2M Jul  2 11:46 wirasinha_mouse_beta_t_1.tsv
2.9M Jul  2 11:46 wirasinha_mouse_beta_t_1.tsv.sampler.tsv
 11M Jul  2 11:46 wirasinha_mouse_beta_t_r.tsv
3.7M Jul  2 11:46 wirasinha_mouse_beta_t_r.tsv.sampler.tsv
 17M Jul  2 11:46 wirasinha_mouse_beta_t_4.tsv
6.1M Jul  2 11:46 wirasinha_mouse_beta_t_4.tsv.sampler.tsv
 11M Jul  2 11:46 wirasinha_mouse_beta_t_8.tsv
3.8M Jul  2 11:46 wirasinha_mouse_beta_t_8.tsv.sampler.tsv
 10M Jul  2 11:46 wirasinha_mouse_beta_s_r.tsv
3.4M Jul  2 11:46 wirasinha_mouse_beta_s_r.tsv.sampler.tsv
6.1M Jul  2 11:46 wirasinha_mouse_beta_g8a.tsv
2.2M Jul  2 11:46 wirasinha_mouse_beta_g8a.tsv.sampler.tsv
 20M Jul  2 11:46 wirasinha_mouse_beta_s_4.tsv
6.9M Jul  2 11:46 wirasinha_mouse_beta_s_4.tsv.sampler.tsv
 16M Jul  2 11:46 wirasinha_mouse_beta_s_8.tsv
5.6M Jul  2 11:46 wirasinha_mouse_beta_s_8.tsv.sampler.tsv
2.5M Jul  2 11:46 wirasinha_mouse_alpha_t_p.tsv
877K Jul  2 11:46 wirasinha_mouse_alpha_t_p.tsv.sampler.tsv
4.4M Jul  2 11:46 wirasinha_mouse_alpha_t_1.tsv
1.5M Jul  2 11:47 wirasinha_mouse_alpha_t_1.tsv.sampler.tsv
5.5M Jul  2 11:47 wirasinha_mouse_alpha_t_r.tsv
1.9M Jul  2 11:47 wirasinha_mouse_alpha_t_r.tsv.sampler.tsv
8.3M Jul  2 11:47 wirasinha_mouse_alpha_t_4.tsv
2.8M Jul  2 11:47 wirasinha_mouse_alpha_t_4.tsv.sampler.tsv
5.5M Jul  2 11:47 wirasinha_mouse_alpha_t_8.tsv
1.9M Jul  2 11:47 wirasinha_mouse_alpha_t_8.tsv.sampler.tsv
5.7M Jul  2 11:47 wirasinha_mouse_alpha_s_r.tsv
1.9M Jul  2 11:47 wirasinha_mouse_alpha_s_r.tsv.sampler.tsv
2.0M Jul  2 11:47 wirasinha_mouse_alpha_g8a.tsv
719K Jul  2 11:47 wirasinha_mouse_alpha_g8a.tsv.sampler.tsv
 15M Jul  2 11:47 wirasinha_mouse_alpha_s_4.tsv
5.3M Jul  2 11:47 wirasinha_mouse_alpha_s_4.tsv.sampler.tsv
8.1M Jul  2 11:47 wirasinha_mouse_alpha_s_8.tsv
2.8M Jul  2 11:47 wirasinha_mouse_alpha_s_8.tsv.sampler.tsv

 Crucially these inputs have the following essential input columns and can be recognized by TCRsampler.clean_mixcr()

    {'subject':'subject',
	'bestVGene':'v_reps',
	'bestJGene':'j_reps', 
	'aaSeqCDR3':'cdr3',
	'cloneCount':'count',
	'cloneFraction':'freq'}

Notes
------
The script splits the (650 MB) file into smaller files based on chain and subset

   subset tcr_b chain chain_name                                  subset_definition                       filename                      strain
0     t_p  poly     b       beta                           pre-selection thymocytes   wirasinha_mouse_beta_t_p.tsv  C57BL6 inbred mouse strain
1     t_1  poly     b       beta                                  wave 1 thymocytes   wirasinha_mouse_beta_t_1.tsv  C57BL6 inbred mouse strain
2     t_r  poly     b       beta                                       thymic T-reg   wirasinha_mouse_beta_t_r.tsv  C57BL6 inbred mouse strain
3     t_4  poly     b       beta                              CD4+ naive thymocytes   wirasinha_mouse_beta_t_4.tsv  C57BL6 inbred mouse strain
4     t_8  poly     b       beta                              CD8+ naive thymocytes   wirasinha_mouse_beta_t_8.tsv  C57BL6 inbred mouse strain
5     s_r  poly     b       beta                                      splenic T-reg   wirasinha_mouse_beta_s_r.tsv  C57BL6 inbred mouse strain
6     g8a  poly     b       beta  small intestinal CD8aa+ intraepithelial lympho...   wirasinha_mouse_beta_g8a.tsv  C57BL6 inbred mouse strain
7     s_4  poly     b       beta                         splenic CD4+ naive T cells   wirasinha_mouse_beta_s_4.tsv  C57BL6 inbred mouse strain
8     s_8  poly     b       beta                         splenic CD8+ naive T cells   wirasinha_mouse_beta_s_8.tsv  C57BL6 inbred mouse strain
9     t_p  poly     a      alpha                           pre-selection thymocytes  wirasinha_mouse_alpha_t_p.tsv  C57BL6 inbred mouse strain
10    t_1  poly     a      alpha                                  wave 1 thymocytes  wirasinha_mouse_alpha_t_1.tsv  C57BL6 inbred mouse strain
11    t_r  poly     a      alpha                                       thymic T-reg  wirasinha_mouse_alpha_t_r.tsv  C57BL6 inbred mouse strain
12    t_4  poly     a      alpha                              CD4+ naive thymocytes  wirasinha_mouse_alpha_t_4.tsv  C57BL6 inbred mouse strain
13    t_8  poly     a      alpha                              CD8+ naive thymocytes  wirasinha_mouse_alpha_t_8.tsv  C57BL6 inbred mouse strain
14    s_r  poly     a      alpha                                      splenic T-reg  wirasinha_mouse_alpha_s_r.tsv  C57BL6 inbred mouse strain
15    g8a  poly     a      alpha  small intestinal CD8aa+ intraepithelial lympho...  wirasinha_mouse_alpha_g8a.tsv  C57BL6 inbred mouse strain
16    s_4  poly     a      alpha                         splenic CD4+ naive T cells  wirasinha_mouse_alpha_s_4.tsv  C57BL6 inbred mouse strain
17    s_8  poly     a      alpha                         splenic CD8+ naive T cells  wirasinha_mouse_alpha_s_8.tsv  C57BL6 inbred mouse strain


	# 1. We only used data from Jackson C57BL6 inbred mouse strain
	# 2. The original data also includes transgenic strains, which were not included, since in all cases subsets 
	# were defined such that tcr_b == 'poly'
	# 3. Currently, TCRsampler can't deal with ambiguity. 
	# 4. Therefore we consider only V and J-genes, striping allele information
	# 5. Moreover, in cases where migec proposes two genes, we pick the first, consistent with bestVGene assignment


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
"""


import pandas as pd
import sys
import re
from tcrsampler.sampler import TCRsampler

wirasinha_definitions = {	't_p' : 'pre-selection thymocytes',
							't_1' : 'wave 1 thymocytes',
							't_r' : 'thymic T-reg',
							't_4' : 'CD4+ naive thymocytes',
							't_8' : 'CD8+ naive thymocytes',
							'g8a' : 'small intestinal CD8aa+ intraepithelial lymphocytes',
							's_r' : 'splenic T-reg',
							's_4' : 'splenic CD4+ naive T cells',
							's_8' : 'splenic CD8+ naive T cells'}

wirasinha_to_mixcr_headers = {	'count'        :'cloneCount',
								'freq'          :'cloneFraction', 
								'cdr3nt'        :'nSeqCDR3',
								'cdr3aa'        :'aaSeqCDR3',
								'bestv'         :'bestVGene',
								'bestj'         :'bestJGene',
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

def subset_wirasinha(df, subset = 't_p', tcr_b = 'poly', chain = 'b'):
	"""
	Subsets the Wirasinha.migec.txt dataframe
	"""
	ind1 = df['subset'] == subset
	ind2 = df['tcr_b'] == tcr_b
	ind3 = df['chain'] == chain
	sdf = df[(ind1&ind2&ind3)]
	return sdf.copy()


def _strip_allele(gene_string):
	"""
	Example
	-------
	>>> _strip_allele('TRBV12-2*01'
	>>> 'TRBV12-2'

	>>> df[['v','j']] = df[['v','j']].apply(lambda x: x.apply(_strip_allele))
	"""
	try:
		regex_groups = re.match(string = gene_string , pattern ='(.*)\*')
		return regex_groups[1]
	except AttributeError:
		return gene_string

def _pick_best(gene_string):
	"""
	Example
	-------
	>>> _pick_best('TRAV12-2*01,TRAV12D-1*01')
    TRAV12-2*01
    >>>_pick_best('TRAV12-2*01'
    'TRAV12-2*01'
	"""
	return gene_string.split(",")[0]

def _add_allele(gene_string, allele = "*01"):
	"""
	Example
	-------
	>>> _add_allele('TRBV12-2')
	'TRBV12-2*01'
	"""
	return f"{gene_string}{allele}"

if __name__ == "__main__":

	# Iterature over out control dataframe, 
	# 1. subsetting wirasinha, 
	# 2. stripping allle
	# 3. Changing headers to match mixcr outputs
	# 4. write file to .tsv
	dfb = pd.DataFrame({'subset' : ["t_p", "t_1", "t_r", "t_4", "t_8","s_r", "g8a", "s_4", "s_8"],
				 'tcr_b'   : ['poly','poly','poly','poly','poly','poly','poly','poly','poly'],
				 'chain'   : ['b','b','b','b','b','b','b','b','b']})
	dfa = pd.DataFrame({'subset' : ["t_p", "t_1", "t_r", "t_4", "t_8","s_r", "g8a", "s_4", "s_8"],
				  'tcr_b'   : ['poly','poly','poly','poly','poly','poly','poly','poly','poly'],
				  'chain'   : ['a','a','a','a','a','a','a','a','a']})
	df = pd.concat([dfb,dfa]).reset_index(drop = True)

	df['chain_name'] = df['chain'].apply(lambda x : {'b':'beta','a':'alpha'}[x])
	df['subset_definition'] = df['subset'].apply(lambda x : wirasinha_definitions [x])
	df['filename'] = 'wirasinha_mouse_' + df['chain_name'] + '_' +  df['subset'] + '.tsv'
	df['strain'] = 'C57BL6 inbred mouse strain'
	print(df)

	wirasinha = pd.read_csv('/Volumes/Samsung_T5/kmayerbl/tcr_data/wirasinha/Wirasinha.migec.txt', sep= '\t')
	for i,row in df.iterrows():
		sdf = subset_wirasinha(df = wirasinha, subset = row['subset'], tcr_b = row['tcr_b'], chain = row['chain'])
		sdf[['bestv','bestj']] = sdf[['v','j']].apply(lambda x: x.apply(_pick_best))
		sdf[['bestv','bestj']] = sdf[['bestv','bestj']].apply(lambda x: x.apply(_strip_allele))
		sdf = sdf.rename(columns = wirasinha_to_mixcr_headers)
		sys.stdout.write(f"Writing {row['filename']}\n")
		sdf.to_csv(row['filename'], sep = "\t")

		sys.stdout.write(f"Testing {row['filename']} for import into TCRsampler\t")
		t = TCRsampler()
		t.clean_mixcr(filename = row['filename'])
		t.build_background()
		print("\n")
		print(t.ref_df.head(3))
		name = f"{row['filename']}.sampler.tsv"
		sys.stdout.write(f"Writing {name} \t")
		t.ref_df.to_csv(name, sep = "\t", index = False)
