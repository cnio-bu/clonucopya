from glob import glob
import os


rule cnv_prep:
    input:
        cnvs = os.path.join(config["samples"]["cnvs"],"{sample}.cnv.tsv")
    output:
        "results/cnv_prep/{sample}_prep.cnv.tsv"
    log:
        "logs/cnv_prep/{sample}.log"
    benchmark:
        "logs/cnv_prep/{sample}.bmk"
    conda:
        "../envs/cnv_prep.yaml"
    params:
        sex = config["samples"]["sex"]
    shell:
        """
        python scripts/cnv_formatting.py --input_file {input} --output_file {output} --sex {params.sex} > {log} 2>&1
    """

