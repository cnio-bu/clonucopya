from glob import glob
import os


rule cnv_prep:
    input:
#        cnvs = os.path.join(config["samples"]["cnvs"],"{sample}.cnv.tsv")
        cnvs = lambda wildcards: samples.loc[wildcards.sample, "cnvs"]
    output:
#        "results/cnv_prep/{sample}_prep.cnv.tsv"
        f"results/{experiment}/cnv_prep/{{sample}}_prep.cnv.tsv"
    log:
        f"logs/cnv_prep/{experiment}/{{sample}}.log"
    benchmark:
        f"logs/cnv_prep/{experiment}/{{sample}}.bmk"
    conda:
        "../envs/cnv_prep.yaml"
    params:
        sex = lambda wildcards: samples.loc[wildcards.sample, "sex"]
    shell:
        """
        python scripts/cnv_formatting.py --input_file {input} --output_file {output} --sex {params.sex} > {log} 2>&1
    """

