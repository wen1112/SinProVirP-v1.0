conda activate VCP_env

snakemake -s /hwfsxx1/ST_HN/P17Z10200N0246/User/xiongliwen/QuantifyPhageProject/snakemake/00.VCP_smk_V1.0/Snakefile.py \
--configfile /hwfsxx1/ST_HN/P17Z10200N0246/User/xiongliwen/QuantifyPhageProject/snakemake/00.VCP_smk_V1.0/test_data/config.yaml \
--jobs 99 --cores 1 --max-threads 16
