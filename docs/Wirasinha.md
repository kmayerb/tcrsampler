



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
python assemble_wirasinha.py
```


