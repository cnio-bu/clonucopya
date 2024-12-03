import pandas as pd
from snakemake.io import expand, touch
from glob import glob 
import os

###### Config file and sample sheets #####

configfile: "config/config_multi_vcf.yaml"


##### Target rules #####

# Get sample names from mutation files
SAMPLES = glob_wildcards(os.path.join(config["samples"]["mutations"], 
                        "{sample}.mut.vcf")).sample

# Define target rule
rule all:
    input:
        expand("results/mutation_prep/{sample}_prep.mut.tsv", sample=SAMPLES)


# FORMATTING MUTATION VCF

include: "rules/01_mutation_prep.smk"
