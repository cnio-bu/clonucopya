rule pvi_vep_prep:
    input:
        pvi_df="results/pyclone-vi/{project}/pvi_out.tsv",
    output:
        dir=directory("results/pvi_vep_prep/{project}")
    params:
        sample_id=lambda wildcards: wildcards.project
    log:
        "logs/pvi_vep_prep/{project}/pvi_vep_prep.log"
    benchmark:
        "logs/pvi_vep_prep/{project}/pvi_vep_prep.smk"
    conda:
        "../envs/pvi_vep_prep.yaml"
    threads: 
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime=240
    shell:
        """
        python scripts/vep_formatting.py --pvi_data {input.pvi_df} \
                                         --sample_id {params.sample_id} \
                                         --out_dir {output.dir} 2> {log}
        """

