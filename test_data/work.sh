conda activate SinProVirP_env

snakemake -s /full/path/to/cloned/repo/Snakefile.py \
--configfile /full/path/to/cloned/repo/test_data/config.yaml \
--jobs 99 --cores 8
