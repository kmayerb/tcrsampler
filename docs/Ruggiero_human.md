# Data Source


## Bioproject PRJNA287162: TCR ligation anchored-magnetically captured PCR (TCR-LA-MC PCR) for TCR α-and β-chain diversity dissection

High-resolution analysis of the human T-cell receptor repertoire

Eliana Ruggiero, Jan P. Nicolay, Raffaele Fronza, Anne Arens, Anna Paruzynski, Ali Nowrouzi, Gökçe Ürenden, Christina Lulay, Sven Schneider, Sergij Goerdt, Hanno Glimm, Peter H. Krammer, Manfred Schmidt & Christof von Kalle 
Nature Communications volume 6, Article number: 8081 (2015) Cite this article

## SRA Study: SRP059581

NCBI : [SRP059581](https://www.ncbi.nlm.nih.gov/sra/?term=SRP059581)

SRA Study : [SRP059581](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP059581)

SRA Browser: [mus musculus](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=1&WebEnv=NCID_1_11090927_130.14.18.48_5555_1593620834_465491811_0MetA0_S_HStore&f=organism_s%3An%3Amus%2520musculus%3Ac&o=acc_s%3Aa)

## Getting only the Health Donor Controls 

```python
 July 2, 2020

# https://www.ncbi.nlm.nih.gov/sra/?term=SRP059581
import os
import pandas as pd 
df = pd.read_csv('SRP059581_sra_result.csv')

ind1 = df['Experiment Title'].apply(lambda x : x.find("Alpha") != -1) 
ind2 = df['Experiment Title'].apply(lambda x : x.find("TCRLAMC PCR Healthy donor") != -1)
ind3 = df['Experiment Title'].apply(lambda x : x.find("PBMC") != -1)
df[(ind1&ind2&ind3)][['Experiment Accession','Experiment Title']]
#     Experiment Accession                                   Experiment Title Sample Accession
# 87            SRX1074237  TCRLAMC PCR Healthy donor 6 PBMC cDNA Alpha chain        SRS973088
# 89            SRX1074235  TCRLAMC PCR Healthy donor 5 PBMC cDNA Alpha chain        SRS973090
# 101           SRX1074222  TCRLAMC PCR Healthy donor 4 PBMC cDNA Alpha chain        SRS973077
# 103           SRX1074220  TCRLAMC PCR Healthy donor 3 PBMC cDNA Alpha chain        SRS973003
# 105           SRX1074218  TCRLAMC PCR Healthy donor 2 PBMC cDNA Alpha chain        SRS973004
# 107           SRX1074216  TCRLAMC PCR Healthy donor 1 PBMC cDNA Alpha chain        SRS973005
selections = df[(ind1&ind2&ind3)][['Experiment Accession','Experiment Title']]
runs = pd.read_csv('SRP059581_SraRunInfo.csv')
selections_url = selections.merge(runs, how = 'left', left_on = ['Experiment Accession'], right_on = ['Experiment'])[['Experiment Accession','Experiment Title','Run','download_path']]
for i,row in selections_url.iterrows():
	cmd = f"wget {row['download_path']}"
	print(cmd)
	os.system(cmd)

# https://www.ncbi.nlm.nih.gov/sra/?term=SRP059581
import os
import pandas as pd 
df = pd.read_csv('SRP059581_sra_result.csv')
#   Experiment Accession                                  Experiment Title         Run                                      download_path
# 0           SRX1074238  TCRLAMC PCR Healthy donor 6 PBMC cDNA Beta chain  SRR2079472  https://sra-downloadb.be-md.ncbi.nlm.nih.gov/s...
# 1           SRX1074236  TCRLAMC PCR Healthy donor 5 PBMC cDNA Beta chain  SRR2079469  https://sra-downloadb.be-md.ncbi.nlm.nih.gov/s...
# 2           SRX1074223  TCRLAMC PCR Healthy donor 4 PBMC cDNA Beta chain  SRR2079457  https://sra-downloadb.be-md.ncbi.nlm.nih.gov/s...
# 3           SRX1074221  TCRLAMC PCR Healthy donor 3 PBMC cDNA Beta chain  SRR2079455  https://sra-downloadb.be-md.ncbi.nlm.nih.gov/s...
# 4           SRX1074219  TCRLAMC PCR Healthy donor 2 PBMC cDNA Beta chain  SRR2079453  https://sra-downloadb.be-md.ncbi.nlm.nih.gov/s...
# 5           SRX1074217  TCRLAMC PCR Healthy donor 1 PBMC cDNA Beta chain  SRR2079451  https://sra-downloadb.be-md.ncbi.nlm.nih.gov/s.
ind1 = df['Experiment Title'].apply(lambda x : x.find("Beta") != -1) 
ind2 = df['Experiment Title'].apply(lambda x : x.find("TCRLAMC PCR Healthy donor") != -1)
ind3 = df['Experiment Title'].apply(lambda x : x.find("PBMC") != -1)
ind4 = df['Experiment Title'].apply(lambda x : x.find("PBMC cDNA Beta chain") != -1)
selections = df[(ind1&ind2&ind3&ind4)][['Experiment Accession','Experiment Title']]
runs = pd.read_csv('SRP059581_SraRunInfo.csv')
selections_url = selections.merge(runs, how = 'left', left_on = ['Experiment Accession'], right_on = ['Experiment'])[['Experiment Accession','Experiment Title','Run','download_path']]
print(selections_url)
for i,row in selections_url.iterrows():
	cmd = f"wget {row['download_path']}"
	print(cmd)
	os.system(cmd)
```


## Convert SRA to fastq

```
import os 

path_to_fastq_dump = '/Users/kmayerbl/sratoolkit.2.9.6-1-mac64/bin/fastq-dump'
path_to_inspect = '/Volumes/Samsung_T5/kmayerbl/tcr_data/ruggiero/SRA/healthy_human_pmbc_alpha/'
files = [f for f in os.listdir(path_to_inspect) if f.endswith('.1')]
for f in files:
	cmd = f"{path_to_fastq_dump} {os.path.join(path_to_inspect,f)} --outdir {path_to_inspect}"
	print(cmd)
	os.system(cmd)

path_to_inspect = '/Volumes/Samsung_T5/kmayerbl/tcr_data/ruggiero/SRA/healthy_human_pmbc_beta/'
files = [f for f in os.listdir(path_to_inspect) if f.endswith('.1')]
for f in files:
	cmd = f"{path_to_fastq_dump} {os.path.join(path_to_inspect,f)} --outdir {path_to_inspect}"
	print(cmd)
	os.system(cmd)
```


## Use nextflow to run mixcr on these files

```groovy 
// Simple Nextflow Routine to run mixcr on a batch of fastq files in a particular folder. 
params.species = 'hsa'
params.input_folder = 'input'
params.output_folder = 'output2'

files = Channel.fromPath(params.input_folder + '/*.fastq')

process mix {

	container 'quay.io/kmayerb/dmixcr:0.0.1'

	publishDir params.output_folder, mode: 'copy', overwrite: false

	errorStrategy 'finish'
		
	input:
	file(fastq) from files

	output:
	file("${fastq}.clns.best.txt")

	script:
	"""
	mixcr align ${fastq} ${fastq}.vdjca --species ${params.species} 
	mixcr assemble ${fastq}.vdjca ${fastq}.clns
	mixcr exportClones -cloneId -count -fraction -vGene -jGene -vHit -jHit -vHits -jHits -aaFeature CDR3 -nFeature CDR3 ${fastq}.clns ${fastq}.clns.best.txt
	"""
}
```

