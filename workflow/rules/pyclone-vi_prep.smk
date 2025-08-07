rule pvi_intesersect:
    input:
        mutations = "results/mutation_prep/{project}/{sample}_prep.mut.tsv",
        cnvs = lambda wildcards: samples.loc[wildcards.sample, "cnvs"]
    output:
        "results/pyclone-vi_prep/{project}/{sample}_intersect_pvi.tsv"
    params:
        sample_id = lambda wildcards: wildcards.sample
    log:
        "logs/pyclone-vi_prep/{project}/{sample}.log"
    benchmark:
        "logs/pyclone-vi_prep/{project}/{sample}.bmk"
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
            "results/pyclone-vi_prep/{project}/{sample}_intersect_pvi.tsv",
            sample=samples_df[samples_df['project'] == wildcards.project]['sample_id'],
            project=wildcards.project
        )
    output:
        pvi="results/pyclone-vi_prep/{project}/combined_intersect_pvi.tsv",
        phyclone_prep="results/pyclone-vi_prep/{project}/pvi_input_phyclone_formatted.tsv"
    params:
        samplesheet=config["samplesheet"]
    log:
        "logs/pyclone-vi_prep/{project}/concat.log"
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

