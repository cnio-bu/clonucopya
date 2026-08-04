def get_input_liftover(wildcards):
    study_to_pvi = globals().get("study_to_pvi", None)

    if study_to_pvi:
        pvi_path = study_to_pvi.get(wildcards.study, None)
        if pvi_path:
            return "results/{study}/pvi-start_prep/pvi-start_prep.tsv".format(
                study=wildcards.study
            )
    return f"results/{wildcards.study}/pyclone-vi_prep/combined_intersect_pvi.tsv"



rule mutation_liftover:
    input:
        get_input_liftover
    output:
        "results/{study}/mutation_liftover/pvi_checked.tsv"
    params:
        genome_reference = config["genome_reference"],
        liftOver = "resources/liftOver/hg19ToHg38.over.chain.gz",
    log:
        "logs/{study}/mutation_liftover/pvi_liftover.log"
    benchmark:
        "logs/{study}/mutation_liftover/pvi_liftover.bmk"
    conda:
        "../envs/mutation_liftover.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/pvi_liftover.py \
              --pvi_input {input} \
              --genome_reference {params.genome_reference} \
              --liftOver {params.liftOver} \
              --output_file {output} > {log} 2>&1
    """
