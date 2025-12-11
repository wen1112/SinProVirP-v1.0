## SinProVirP: Signature-Protein based Virome Profiling tool
![Version](https://img.shields.io/badge/version-1.0-blue)
![last commit](https://img.shields.io/badge/last_commit-2025.03.05-blue)
![GitHub](https://img.shields.io/github/license/wen1112/SinProVirP-v1.0)
![Badge](https://hitscounter.dev/api/hit?url=https%3A%2F%2Fgithub.com%2Fwen1112%2FSinProVirP-v1.0&label=Visitors&icon=github&color=%23198754&message=&style=flat&tz=UTC)
![GitHub all releases](https://img.shields.io/github/downloads/wen1112/SinProVirP-v1.0/total)
![GitHub stars](https://img.shields.io/github/stars/wen1112/SinProVirP-v1.0?style=social)






### A signature protein-based tool developed for genus-level profiling of the human gut virome.

<img width="868" alt="image" src="https://github.com/user-attachments/assets/25e4a88a-2b57-4f6f-919e-e0b416024fbd">


## Citation
[SinProVirP](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(25)00286-3)
```
@article{yang2025signature,
  title={A signature-protein-based approach for accurate and efficient profiling of the human gut virome},
  author={Yang, Fangming and Xiong, Liwen and Li, Min and Feng, Xuyang and Ren, Huahui and Shi, Zhun and Zhong, Huanzi and Li, Junhua},
  journal={Cell Reports Methods},
  year={2025},
  publisher={Elsevier},
  doi={10.1016/j.crmeth.2025.101250},
  url={https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(25)00286-3}
}
```

## Installation  

- **Step One - Clone the repository**  
Clone the repository to the desired location on your system using the following command:
```
git clone https://github.com/wen1112/SinProVirP-v1.0.git
```

- **Step Two - Create a new conda environment**  
Create a new conda environment via the following command, replacing /full/path/to/cloned/repo with the appropriate path to your cloned repository:
```
conda env create -n SinProVirP_env --file /full/path/to/cloned/repo/SinProVirP_env.yaml
```

Activate the environment by executing:
```
conda activate SinProVirP_env
```

If you have trouble creating the environment using the above commands, you can alternatively follow the instructions to create a enviroment.
```
conda create -n SinProVirP_env python=3.11
conda activate SinProVirP_env

conda install bioconda::snakemake
conda install samtools=1.6
conda install bedtools=2.31.0
conda install diamond=2.0.14
conda install biopython
```

If you do not already have conda installed, please install using the instructions provided [here](https://developers.google.com/earth-engine/guides/python_install-conda/).  
If you have some problem with snakemake install, please install using the instructions provided [here](https://snakemake.readthedocs.io/en/v7.25.0/getting_started/installation.html)

## Test Your Installation
### Input file

+ If the input file is Paired-end, then it will 'cat R1.fq R2.fq > R1R2.fq' as the final input file because of diamond required SE reads input.
+ The input file each col should split by '\t'.
+ Paired-end reads
  |  sample_id   |     fq1      |     fq2     |
  | :----------: |  :--------:  | :--------:  |
  |     s1       |  s1.R1.fq.gz | s1.R2.fq.gz |
  |     s2       |  s2.R1.fq.gz | s2.R2.fq.gz |
  |     s3       |  s3.R1.fq.gz | s3.R2.fq.gz |
  |     s4       |  s4.R1.fq.gz | s4.R2.fq.gz |

+ Also you can 'cat R1.fq R2.fq > R1R2.fq' before you input file.
  |  sample_id   |       fq       |
  | :----------: |  :----------:  |
  |     s1       |  s1.R1R2.fq.gz |
  |     s2       |  s2.R1R2.fq.gz |
  |     s3       |  s3.R1R2.fq.gz |
  |     s4       |  s4.R1R2.fq.gz |
 
+ Single-end reads
  |  sample_id   |     fq1      |    fq2     |
  | :----------: |  :--------:  | :--------: |
  |     s1       |  s1.R1.fq.gz |            |
  |     s2       |  s2.R1.fq.gz |            |
  |     s3       |  s3.R1.fq.gz |            |
  |     s4       |  s4.R1.fq.gz |            |

  
### Quick start

+ We need three file: **'work.sh', 'config.yaml' and 'sample_file.txt'**
+ Please remeber to **change the path in it**.  

```
cd /full/path/to/cloned/repo/test_data
bash work.sh
```

+ work.sh:
```
conda activate SinProVirP_env
snakemake -s /full/path/to/cloned/repo/Snakefile.py \
--configfile /full/path/to/cloned/repo/test_data/config.yaml \
--jobs 99 --cores 8 
```
+ config.yaml
```
pipeline_directory: /full/path/to/cloned/repo

# Sample file specifies sample names and names of files containing sample reads
# Format: Tab-delimited, two columns
# sample_name  reads_R1R2file
# See example (samp_file.txt) in the testing folder
sample_file: /full/path/to/sample/file

# In which directory should results be outputted?
outdir: /full/path/to/output/directory
```

+ sample_file.txt - split each col by **"\t"**.  
```
sample1 /full/path/to/cloned/repo/test_data/test1_R1R2.fastq
sample2 /full/path/to/cloned/repo/test_data/test2_R1R2.fastq
```


+ The script output log may just like this:
```
Building DAG of jobs...
Relative file path './test2_R1R2.fastq' starts with './'. This is redundant and strongly discouraged. It can also lead to inconsistent results of the file-matching approach used by Snakemake. You can simply omit the './' for relative file paths.
Using shell: /bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       add_taxonomy
        1       all
        1       calculate_relative_abudance
        1       diamond_blastx
        1       filter_mapped_reads
        1       get_marker_coverage
        1       index_bam_file
        1       reads_per_marker_protein
        1       sam_to_bam
        1       sort_bam_file
        10

2024-05-30 22:03:25.334570
The sample ./test_data/sample1/sample1.sorted.bam.idxstats.addTaxonomy.final.abundance.txt of SinProVirP pipline is done!
```

### Output file

```
./
├── sample1.bam
├── sample1.sam
├── sample1.sorted.bai
├── sample1.sorted.bam
├── sample1.sorted.bam.coverage
├── sample1.sorted.bam.filter.coverage.idxstats
├── sample1.sorted.bam.idxstats
├── sample1.sorted.bam.idxstats.abundance
└── sample1.sorted.bam.idxstats.addTaxonomy.final.abundance.txt

0 directories, 9 files
```

+ *.final.abundance.txt

  |    VCID   | Relative_abundance | VC_lineage | VC_lifestyle | VC_host_lineage |
  | :-------: |  :--------------:  | :--------------: | :--------------: |  :--------------: |
  | VC_3803_0 | 0.0123774546010874 | Viruses;Uroviricota;Caudoviricetes;Caudovirales;Siphoviridae;; |  Temperate | Bacteria;Firmicutes;Clostridia;Eubacteriales;Oscillospiraceae;Ruminococcus; |
  | VC_2046_0 | 0.0362526365705308 | Viruses;Uroviricota;Caudoviricetes;Caudovirales;Myoviridae;; | Temperate | Bacteria;Firmicutes;Clostridia;Eubacteriales;Oscillospiraceae;Faecalibacterium; |
  | VC_2673_0 | 0.0140144228362927 | Viruses;Uroviricota;Caudoviricetes;;;; | Unknown | Bacteria;Firmicutes;Clostridia;Eubacteriales;Oscillospiraceae;Oscillibacter; |
  | VC_6076_0 | 0.0063759751235582 | Viruses;Uroviricota;Caudoviricetes;Caudovirales;Myoviridae;; | Virulent | Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides; |


## About parameters
+ If you want to adjust the parameters about SinProVirP pipeline used, you could change the 'config.yaml' file. 
```
# The value setting
threads: 8 # see usage instructions - can increase if you have more threads available
IDENTITY: 90 # diamond balst used identity
COVERAGE: 80 # diamond balst used coverage
MARKER_COVERAGE: 0.7 # SinProVirP profiling used marker coverage
MARKER_RATIO: 0.5 # SinProVirP profiling used marker ratio
```

## DCS Cloud avaliable
* https://cloud.stomics.tech/library/#/tool/detail/workflow/WF09202503109bppPb/1.0.0?zone=st
* DCS Help: https://www.stomics.tech/helpcenter
