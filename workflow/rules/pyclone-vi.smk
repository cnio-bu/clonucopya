rule pyclone_vi:
    input:
        "results/{study}/mutation_liftover/pvi_checked.tsv"
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
