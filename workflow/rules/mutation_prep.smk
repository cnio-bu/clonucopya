rule mutation_prep:
    input:
        mutations = lambda wildcards: samples.loc[wildcards.sample, "mutations"]
    output:
        "results/{study}/mutation_prep/{sample}_prep.mut.tsv"
    params:
        snv_filter = config["just_snv"]
    log:
        "logs/{study}/mutation_prep/{sample}.log"
    benchmark:
        "logs/{study}/mutation_prep/{sample}.bmk"
    conda:
        "../envs/mutation_prep.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/mutations_formatting.py --input_vcf {input} --just_snv {params.snv_filter} --output_file {output} > {log} 2>&1
    """
