from glob import glob
import os


rule mutation_prep:
    input:
        mutations = lambda wildcards: samples.loc[wildcards.sample, "mutations"]
    output:
        f"results/{experiment}/mutation_prep/{{sample}}_prep.mut.tsv"
    log:
        f"logs/mutation_prep/{{sample}}.log"
    benchmark:
        f"logs/mutation_prep/{{sample}}.bmk"
    conda:
        "../envs/mutation_prep.yaml"
    shell:
        """
        python scripts/mutations_formatting.py --input_vcf {input} --output_file {output} > {log} 2>&1
    """
