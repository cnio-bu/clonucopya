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



rule check_mutations:
    input:
        lambda wildcards: expand(
            "results/{study}/mutation_prep/{sample}_prep.mut.tsv",
            sample=samples_df[samples_df['study'] == wildcards.study]['sample_id'],
            study=wildcards.study,
            allow_missing=True
        )
    output:
        "results/{study}/mutation_prep/bam_checked/{sample}_check.mut.tsv"
    params:
        mut_path = directory("results/{study}/mutation_prep"),
        bam_path = config["bam_files"],
        outdir = directory("results/{study}/mutation_prep/bam_checked")
    log:
        "logs/{study}/mutation_prep/{sample}_check_mutations.log"
    benchmark:
        "logs/{study}/mutation_prep/{sample}_check_mutations.bmk"
    conda:
        "../envs/mutation_prep.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/bam_check.py --mut_path {params.mut_path} \\
                                    --bam_path {params.bam_path} \\
                                    --output_dir {params.outdir} \\
                                    > {log} 2>&1
        """
