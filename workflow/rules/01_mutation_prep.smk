from glob import glob
import os


rule mutation_prep:
    input:
        mutations = os.path.join(config["samples"]["mutations"],"{sample}.mut.vcf")
    output:
        "results/mutation_prep/{sample}_prep.mut.tsv"
    log:
        "logs/mutation_prep/{sample}.log"
    benchmark:
        "logs/mutation_prep/{sample}.bmk"
    conda:
        "../envs/mutation_prep.yaml"
    shell:
        """
        python scripts/mutations_formatting.py --input_vcf {input} --output_file {output} > {log} 2>&1
    """
