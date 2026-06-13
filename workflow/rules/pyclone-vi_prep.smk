def get_check_mutations_input(wildcards):
    if config.get("bam_check", False):
        return "results/{study}/mutation_prep/bam_checked/{sample}_check.mut.tsv".format(
            study=wildcards.study,
            sample=wildcards.sample
        )
    else:
        return "results/{study}/mutation_prep/{sample}_prep.mut.tsv".format(
            study=wildcards.study,
            sample=wildcards.sample
        )


rule pvi_intesersect:
    wildcard_constraints:
        sample = "|".join(samples["sample_id"].tolist())
    input:
        mutations = get_check_mutations_input,
        cnas = lambda wildcards: samples.at[(wildcards.study, wildcards.sample), "cnas"]
    output:
        "results/{study}/pyclone-vi_prep/{sample}_intersect_pvi.tsv"
    params:
        sample_id = lambda wildcards: wildcards.sample
    log:
        "logs/{study}/pyclone-vi_prep/{sample}.log"
    benchmark:
        "logs/{study}/pyclone-vi_prep/{sample}.bmk"
    conda:
        "../envs/intersect_mutations_cna.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/intersect_mutations_cna.py \
              --sample_id {params.sample_id} \
              --mutations {input.mutations} \
              --cnas {input.cnas} \
              --output_file {output} > {log} 2>&1
    """






rule format_pvi_intersect:
    input:
        lambda wildcards: expand(
            "results/{study}/pyclone-vi_prep/{sample}_intersect_pvi.tsv",
            sample=samplesheet[samplesheet['study'] == wildcards.study]['sample_id'],
            study=wildcards.study
        )
    output:
        pvi="results/{study}/pyclone-vi_prep/combined_intersect_pvi.tsv",
    params:
        samplesheet=config["samplesheet"]
    log:
        "logs/{study}/pyclone-vi_prep/concat.log"
    benchmark:
        "logs/{study}/pyclone-vi_prep/concat.bmk"
    conda:
        "../envs/intersect_mutations_cna.yaml"
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
        > {log} 2>&1
        """
