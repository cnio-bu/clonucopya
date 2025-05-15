from glob import glob
import os


rule mutation_prep:
    input:
        mutations = lambda wildcards: samples.loc[wildcards.sample, "mutations"]
    output:
        "results/mutation_prep/{project}/{sample}_prep.mut.tsv"
    log:
        "logs/mutation_prep/{project}/{sample}.log"
    benchmark:
        "logs/mutation_prep/{project}/{sample}.bmk"
    conda:
        "../envs/mutation_prep.yaml"
    shell:
        """
        python scripts/mutations_formatting.py --input_vcf {input} --output_file {output} > {log} 2>&1
    """
