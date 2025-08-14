rule pvi_vep_prep:
    input:
        pvi_prep="results/{study}/pyclone-vi_prep/combined_intersect_pvi.tsv",
        pvi_results="results/{study}/pyclone-vi/pvi_out.tsv"
    output:
        dir=directory("results/{study}/pvi_vep_prep")
    params:
        study=lambda wildcards: wildcards.study
    log:
        "logs/{study}/pvi_vep_prep/pvi_vep_prep.log"
    benchmark:
        "logs/{study}/pvi_vep_prep/pvi_vep_prep.bmk"
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

