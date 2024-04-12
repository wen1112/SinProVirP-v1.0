# VCP
## Accurate and Rapid Metagenomic Virome Profiling using a Marker Protein-based Approach
### A novel marker-based virome profiling approach based on “viral clusters (VCs)-marker Protein Clusters (mPCs)”, which show higher accuracy, comprehensiveness and efficiency, is considerably support large-scale human virome studies to decipher the relationship between gut viruses and human health and disease.

<img width="440" alt="image" src="https://github.com/wen1112/VCP_v1.0/assets/71487119/f95bb9a5-3985-414e-a8c6-7cdbb9d05039">

## Citation
xxxxxxxxxxxxxxxxxxx  

## Installation  

- **Step One - Clone the repository**  
Clone the repository to the desired location on your system using the following command:
```
git clone https://github.com/wen1112/VCP_v1.0.git
```

- **Step Two - Create a new conda environment**  
Create a new conda environment via the following command, replacing /full/path/to/cloned/repo with the appropriate path to your cloned repository:
```
conda env create -n VCP_env --file /full/path/to/cloned/repo/vcp_env.yaml
```

Activate the environment by executing:
```
conda activate VCP_env
```

If you have trouble creating the environment using the above commands, you can alternatively follow the instructions to create a enviroment.
```
conda create -n VCP_env
conda activate VCP_env

conda install -c conda-forge python=3.9
conda install -c conda-forge biopython=1.81
conda install -c conda-forge -c bioconda snakemake=7.25.0
conda install -c conda-forge samtools=1.6
conda install -c conda-forge bedtools=2.31.0
conda install -c conda-forge diamond=2.0.14
```

If you do not already have conda installed, please install using the instructions provided [here](https://developers.google.com/earth-engine/guides/python_install-conda/).  
If you have some problem with snakemake install, please install using the instructions provided [here](https://snakemake.readthedocs.io/en/v7.25.0/getting_started/installation.html)

## Test Your Installation
### Iuput file

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

+ Please remeber to change the path in 'config.yaml' & 'sample_file.txt'.
```
cd /full/path/to/cloned/repo/test_data
bash work.sh
```

work.sh:
```
snakemake -s /full/path/to/cloned/repo/Snakefile.py \
--configfile /full/path/to/cloned/repo/test_data/config.yaml \
--jobs 99 --cores 1 --max-threads 16
```


### Output file
+ *.abundance

