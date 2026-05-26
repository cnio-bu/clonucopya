rule pvi_intesersect:
    wildcard_constraints:
        sample = "|".join(samples["sample_id"].tolist())
    input:
        mutations = "results/{study}/mutation_prep/bam_checked/{sample}_check.mut.tsv",
        cnvs = lambda wildcards: samples.loc[wildcards.sample, "cnvs"]
    output:
        "results/{study}/pyclone-vi_prep/{sample}_intersect_pvi.tsv"
    params:
        sample_id = lambda wildcards: wildcards.sample
    log:
        "logs/{study}/pyclone-vi_prep/{sample}.log"
    benchmark:
        "logs/{study}/pyclone-vi_prep/{sample}.bmk"
    conda:
        "../envs/intersect_mutations_cnv.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/intersect_mutations_cnv.py \
              --sample_id {params.sample_id} \
              --mutations {input.mutations} \
              --cnvs {input.cnvs} \
              --output_file {output} > {log} 2>&1
    """






rule format_pvi_intersect:
    input:
        lambda wildcards: expand(
            "results/{study}/pyclone-vi_prep/{sample}_intersect_pvi.tsv",
            sample=samples_df[samples_df['study'] == wildcards.study]['sample_id'],
            study=wildcards.study
        )
    output:
        pvi="results/{study}/pyclone-vi_prep/combined_intersect_pvi.tsv",
        phyclone_prep="results/{study}/pyclone-vi_prep/pvi_input_phyclone_formatted.tsv"
    params:
        samplesheet=config["samplesheet"]
    log:
        "logs/{study}/pyclone-vi_prep/concat.log"
    benchmark:
        "logs/{study}/pyclone-vi_prep/concat.bmk"
    conda:
        "../envs/intersect_mutations_cnv.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/format_pvi_intersect.py \
            --intersect_list {input} \
            --samplesheet {params.samplesheet} \
            --pvi_prep {output.pvi} \
            --phyclone {output.phyclone_prep} \
        > {log} 2>&1
        """

