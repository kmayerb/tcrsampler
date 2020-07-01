

## Bioproject PRJNA287162: TCR ligation anchored-magnetically captured PCR (TCR-LA-MC PCR) for TCR α-and β-chain diversity dissection

High-resolution analysis of the human T-cell receptor repertoire

Eliana Ruggiero, Jan P. Nicolay, Raffaele Fronza, Anne Arens, Anna Paruzynski, Ali Nowrouzi, Gökçe Ürenden, Christina Lulay, Sven Schneider, Sergij Goerdt, Hanno Glimm, Peter H. Krammer, Manfred Schmidt & Christof von Kalle 
Nature Communications volume 6, Article number: 8081 (2015) Cite this article

## SRA Study: SRP059581

NCBI : [SRP059581](https://www.ncbi.nlm.nih.gov/sra/?term=SRP059581)

SRA Study : [SRP059581](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP059581)

SRA Browser: [mus musculus](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=1&WebEnv=NCID_1_11090927_130.14.18.48_5555_1593620834_465491811_0MetA0_S_HStore&f=organism_s%3An%3Amus%2520musculus%3Ac&o=acc_s%3Aa)

### SRA toolkit conversion to fastq

```
cd sratoolkit.2.9.6-1-mac64
./fastq-dump /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079521/SRR2079521.1 --outdir /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079521/
./fastq-dump /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079522/SRR2079522.1 -outdir /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079522/
```

### mixcr clonotyping 

#### Beta Chain

```
cd /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079522/
docker run -v /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079522/:/work -it milaboratory/mixcr:3-imgt
mixcr align SRR2079522.1.fastq SRR2079522.1.vdjca --species mmu
mixcr assemble SRR2079522.1.vdjca SRR2079522.1.clns
mixcr exportClones SRR2079522.1.clns SRR2079522.1.clns.txt
mixcr exportAlignments SRR2079522.1.vdjca  SRR2079522.1.vdjca.txt
```

### Alpha Chain
```
cd /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079521/
docker run -v /Volumes/Samsung_T5/kmayerbl/tcr_data/SRP059581/SRR2079521/:/work -it milaboratory/mixcr:3-imgt
mixcr align SRR2079521.1.fastq SRR2079521.1.vdjca --species mmu
mixcr assemble SRR2079521.1.vdjca SRR2079521.1.clns
mixcr exportClones SRR2079521.1.clns SRR2079521.1.clns.txt
mixcr exportAlignments SRR2079521.1.vdjca  SRR2079521.1.vdjca.txt
```
