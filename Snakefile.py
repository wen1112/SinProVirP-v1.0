'''
2024-05-30
VCP v1.0
This is a snakemake rule file script for the VCP pipline.
'''

# Step 1 - Import modules
import snakemake
import sys
import os
from os.path import join
import logging
import subprocess
from datetime import datetime


#logger = logging.getLogger(__name__)
#logging.basicConfig(filename='vcp.log', level=logging.INFO)

# Step 2 - Get information from the config file & read the input file list
configfile: "config.yaml"

def get_sample_reads(sample_file,output_dir):
    sample_reads = {}
    with open (sample_file, 'r') as sf:
        for line in sf:
            line = line.rstrip('\n').split('\t')
            if (len(line) == 2): # single end
                reads = line[1]
                sample = line[0]
                os.makedirs(f"{output_dir}/{sample}", exist_ok=True)
                sample_reads[sample] = reads
                #logger.info("Your input file is single end, please make sure the sample file is - cat R1.fq R2.fq > R1R2.fq")
                if not os.path.exists(reads):
                    sys.exit('No such file or directory: ' + reads)
            elif (len(line) == 3): # paired end specified
                sample = line[0]
                os.makedirs(f"{output_dir}/{sample}", exist_ok=True)
                if not (os.path.exists(line[1]) and os.path.exists(line[2])):
                    sys.exit('No such file or directory, please check your input: ' + line[1] + ' or ' + line[2])

                subprocess.run(f"cat {line[1]} {line[2]} > {output_dir}/{sample}/{sample}_R1R2.fq", shell=True)
                reads = f"{output_dir}/{sample}/{sample}_R1R2.fq"
                sample_reads[sample] = reads
                #logger.info("Because your input data is paired end, so we cat R1.fq and R2.fq to R1R2.fq.")
            else :
                sys.exit('There is no required sample file, please check your input.')

    return (sample_reads)               

# get output directory
outdir = config['outdir']

# call the function, get the sample reads
sample_reads = get_sample_reads(config['sample_file'],outdir)
sample_names = sample_reads.keys()


# Step 3 - Define the rules
rule all:
    input:
        expand(f"{outdir}/{{samp}}/{{samp}}.sorted.bam.idxstats.addTaxonomy.final.abundance.txt", samp=sample_names)

rule diamond_blastx:
    input:
        reads = lambda wildcards: sample_reads[wildcards.samp]
    output:
        sam = f"{outdir}/{{samp}}/{{samp}}.sam"
    params:
        pipe_directory = config['pipeline_directory'],
        ID = config['ID'],
        COVER = config['COVER']
    threads: config['threads']
    shell:"""
        diamond blastx --db {params.pipe_directory}/database/PhageMarkerProtein.V6.Diamond.dmnd \
        -q {input.reads} -o {output.sam} --max-target-seqs 1 --outfmt 101 --evalue 1e-6 --unal 0 \
        --id {params.ID} --query-cover {params.COVER}
        """

rule sam_to_bam:
    input:
        sam = f"{outdir}/{{samp}}/{{samp}}.sam"
    output:
        bam = f"{outdir}/{{samp}}/{{samp}}.bam"
    params: pipe_directory = config['pipeline_directory']
    shell: """
        samtools view -bSh -T {params.pipe_directory}/database/all.324056.cancidate.pc.subdist.n324056.representPR.rmPhageHomo.rmMarkerHomo.addRecovered.rmBacHMM.Shinkage.above3.sameVC.faa {input.sam} > {output.bam}
        """

rule sort_bam_file:
        input:
            bam = f"{outdir}/{{samp}}/{{samp}}.bam"
        output:
            sort_bam = f"{outdir}/{{samp}}/{{samp}}.sorted.bam"
        threads: config['threads']
        shell:
            "samtools sort -@ {threads} {input.bam} -o {output.sort_bam}"

rule index_bam_file:
        input:
            sort_bam = f"{outdir}/{{samp}}/{{samp}}.sorted.bam"
        output:
            bai = f"{outdir}/{{samp}}/{{samp}}.sorted.bai"
        threads: config['threads']
        shell:
            "samtools index -@ {threads} {input.sort_bam} {output.bai}"

rule reads_per_marker_protein:
        input:
            sort_bam = f"{outdir}/{{samp}}/{{samp}}.sorted.bam",
            bai = f"{outdir}/{{samp}}/{{samp}}.sorted.bai"
        output:
            idxstats = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.idxstats"
        shell:
            "samtools idxstats {input.sort_bam} > {output.idxstats}"

rule get_marker_coverage:
        input:
            sort_bam = f"{outdir}/{{samp}}/{{samp}}.sorted.bam",
            bai = f"{outdir}/{{samp}}/{{samp}}.sorted.bai"
        output:
            coverage = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.coverage"
        shell:
            "bedtools genomecov -ibam {input.sort_bam} > {output.coverage}"

rule filter_mapped_reads:
        input:
            coverage = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.coverage",
            idxstats = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.idxstats"
        output:
            filter_idxstats = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.filter.coverage.idxstats"
        params:
            MARKER_COVERAGE = config['MARKER_COVERAGE'],
            pipe_directory = config['pipeline_directory']
        shell: """
        perl {params.pipe_directory}/script/filter.coverage.pl \
        {input.coverage} {input.idxstats} {output.filter_idxstats} {params.MARKER_COVERAGE}
        """

rule calculate_relative_abudance:
        input:
            filter_idxstats = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.filter.coverage.idxstats"
        output:
            abundance = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.idxstats.abundance"
        params:
            pipe_directory = config['pipeline_directory'],
            MARKER_RATIO = config['MARKER_RATIO']
        shell: """
            python {params.pipe_directory}/script/reads_abundance.V6.py \
            {input.filter_idxstats} {params.pipe_directory}/database/all.324056.cancidate.pc.subdist.n324056.representPR.rmPhageHomo.rmMarkerHomo.addRecovered.rmBacHMM.Shinkage.above3.sameVC \
            {params.MARKER_RATIO} {output.abundance}
            """

rule add_taxonomy:
        input:
            abundance = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.idxstats.abundance",
        output:
            final = f"{outdir}/{{samp}}/{{samp}}.sorted.bam.idxstats.addTaxonomy.final.abundance.txt"
        params:
            pipe_directory = config['pipeline_directory']
        shell: """
            python {params.pipe_directory}/script/add_taxonomy.py \
            {params.pipe_directory}/database/20230821_VC_tax_lifestyle_host.txt \
            {input.abundance} {output.final}
            """


#logger.info("The VCP pipline is over!")
