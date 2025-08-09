rule pvi_vep_prep:
    input:
        pvi_prep="results/pyclone-vi_prep/{project}/combined_intersect_pvi.tsv",
        pvi_results="results/pyclone-vi/{project}/pvi_out.tsv"
    output:
        dir=directory("results/pvi_vep_prep/{project}")
    params:
        study=lambda wildcards: wildcards.project
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
        python scripts/vep_formatting.py --pvi_prep {input.pvi_prep} \
                                         --pvi_data {input.pvi_results} \
                                         --study {params.study} \
                                         --out_dir {output.dir} 2> {log}
        """

