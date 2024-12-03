rule pre_pyclone-vi:
    input:
        SNVs = "????/annotated/snv_{CASES}.vcf",
        CNVs = "????/annotated/cnv_{CASES}.tsv",
    output:
        "results/pre_pvi/{PATIENT}/input_{CASES}.tsv",
    log:
        "logs/pre_pvi/{PATIENT}_{CASES}.log",
    benchmark:
        "logs/pre_pvi/{PATIENT}_{CASES}.smk",
    conda:
        "envs/commons.yaml",
    threads:
         config["resources"]["pyclone-vi"]["threads"],
    resources:
        mem = config["resources"]["pyclone-vi"]["mem"],
        walltime=config["resources"]["pyclone-vi"]["walltime"],
    params:
        sample = "{CASES}" # Arreglarlo en config,revisar con el Snakefile original
    shell:
        """
        python scripts/cross_snv_cnv.py {params.sample} {input.SNVs} {input.CNVs} {output} > {log} 2>&1
    """
