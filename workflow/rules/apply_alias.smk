def get_pyclone_input(wildcards):
    study_to_pvi = globals().get("study_to_pvi", None)

    if study_to_pvi:
        pvi_path = study_to_pvi.get(wildcards.study, None)
        if pvi_path:
            return pvi_path.format(study=wildcards.study)

    return f"results/{wildcards.study}/pyclone-vi_prep/combined_intersect_pvi.tsv"


rule apply_alias:
    input:
        get_pyclone_input,
    output:
        "results/{study}/apply_alias/pvi_alias_swap.tsv",
    params:
        metadata = lambda wildcards: samplesheet.loc[
            samplesheet["study"] == wildcards.study, "metadata"
        ].iloc[0],
    log:
        "logs/{study}/apply_alias/alias_swap.log",
    benchmark:
        "logs/{study}/apply_alias/alias_swamp.bmk",
    conda:
        "../envs/apply_alias.yaml",
    threads:
        config["resources"]["default"]["threads"],
    resources:
        mem_mb  = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"],
    shell:
        """
        python scripts/apply_alias.py \
            --pvi         {input} \
            --metadata    {params.metadata} \
            --output_file {output} \
            > {log} 2>&1
        """