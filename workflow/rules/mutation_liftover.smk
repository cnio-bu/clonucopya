rule mutation_liftover:
    input:
        "results/{study}/apply_alias/pvi_alias_swap.tsv"
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