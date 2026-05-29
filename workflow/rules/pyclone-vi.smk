def get_pyclone_input(wildcards):
    study_to_pvi = globals().get("study_to_pvi", None)

    if study_to_pvi:
        pvi_path = study_to_pvi.get(wildcards.study, None)
        if pvi_path:
            return pvi_path.format(study=wildcards.study)

    return f"results/{wildcards.study}/pyclone-vi_prep/combined_intersect_pvi.tsv"



rule pyclone_vi:
    input:
        get_pyclone_input
    output:
        fit = "results/{study}/pyclone-vi/pvi_out.h5",
        result = "results/{study}/pyclone-vi/pvi_out.tsv"
    params:
        nclusters = config["params"]["pyclone-vi"]["num_clusters"],
        density = config["params"]["pyclone-vi"]["density"],
        ngrid = config["params"]["pyclone-vi"]["num_grid_points"],
        nrestarts = config["params"]["pyclone-vi"]["num_restarts"],
        seed = config["params"]["pyclone-vi"]["seed"]
    log:
        "logs/{study}/pyclone-vi/pvi.log"
    benchmark:
        "logs/{study}/pyclone-vi/pvi.bmk"
    conda:
        "../envs/pyclone-vi.yaml"
    threads: 
        config["resources"]["pyclone-vi"]["threads"]
    resources:
        mem_mb = config["resources"]["pyclone-vi"]["mem"],
        runtime = config["resources"]["pyclone-vi"]["walltime"]
    shell:
        """
        pyclone-vi fit -i {input} -o {output.fit} \
            -c {params.nclusters} \
            -d {params.density} \
            -r {params.nrestarts} \
            -g {params.ngrid} \
            --seed {params.seed} > {log} 2>&1 && \
        pyclone-vi write-results-file -i {output.fit} -o {output.result} >> {log} 2>&1
        """
