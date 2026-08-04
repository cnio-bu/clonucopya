def get_pyclone_input(wildcards):
    study_to_pvi = globals().get("study_to_pvi", None)

    if study_to_pvi:
        pvi_path = study_to_pvi.get(wildcards.study, None)
        if pvi_path:
            return pvi_path.format(study=wildcards.study)

    return f"results/{wildcards.study}/pyclone-vi_prep/combined_intersect_pvi.tsv"


rule pvi_start_prep:
    input:
        get_pyclone_input,
    output:
        "results/{study}/pvi-start_prep/pvi-start_prep.tsv",
    params:
        snv_filter = config["just_snv"],
        metadata = lambda wildcards: samplesheet.loc[
            samplesheet["study"] == wildcards.study, "metadata"
        ].iloc[0],
    log:
        "logs/{study}/pvi-start_prep/pvi-start_prep.log",
    benchmark:
        "logs/{study}/pvi-start_prep/pvi-start_prep.bmk",
    conda:
        "../envs/pvi-start_prep.yaml",
    threads:
        config["resources"]["default"]["threads"],
    resources:
        mem_mb  = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"],
    shell:
        """
        python scripts/pvi-start_prep.py \
            --pvi {input} \
            --just_snv {params.snv_filter} \
            --metadata    {params.metadata} \
            --output_file {output} \
            > {log} 2>&1
        """