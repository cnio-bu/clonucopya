rule pvi_vep_prep:
    input:
        pvi_df="results/pyclone-vi/{sample}_pvi_out.tsv",
    output:
        vep_input="results/pvi_vep_prep/{sample}_cluster_0.tsv"
#        dir=directory("results/pvi_vep_prep/{sample}")
    params:
        dir=directory("results/pvi_vep_prep"),
        sample_id=lambda wildcards: wildcards.sample
    log:
        "logs/pvi_vep_prep/{sample}.log"
    benchmark:
        "logs/pvi_vep_prep/{sample}.smk"
    conda:
        "../envs/pvi_vep_prep.yaml"
    threads: 
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/vep_formatting.py --pvi_data {input.pvi_df} \
                                         --sample_id {params.sample_id} \
                                         --out_dir {params.dir} 2> {log}
        """
