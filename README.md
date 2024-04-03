# VCP
## Accurate and Rapid Metagenomic Virome Profiling using a Marker Protein-based Approach
### A novel marker-based virome profiling approach based on “viral clusters (VCs)-marker Protein Clusters (mPCs)”, which show higher accuracy, comprehensiveness and efficiency, is considerably support large-scale human virome studies to decipher the relationship between gut viruses and human health and disease.

<img width="440" alt="image" src="https://github.com/wen1112/VCP_v1.0/assets/71487119/f95bb9a5-3985-414e-a8c6-7cdbb9d05039">

## Citation
xxxxxxxxxxxxxxxxxxx  

## Quick Start
### Installation
- **Step One - Clone the repository**  
Clone the repository to the desired location on your system using the following command:
```
git clone
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
conda create -n mvp
conda activate mvp

conda install python (testing on 3.9)
conda install biopython
conda install snakemake
conda install samtools
conda install bedtools
conda install diamond
```

If you do not already have conda installed, please install using the instructions provided [here](https://developers.google.com/earth-engine/guides/python_install-conda/).


### Test Your Installation



