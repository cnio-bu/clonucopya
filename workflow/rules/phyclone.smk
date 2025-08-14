rule phyclone:
    input:
        pvi_input="results/{study}/pyclone-vi_prep/pvi_input_phyclone_formatted.tsv",
        pvi_output="results/{study}/pyclone-vi/pvi_out.tsv"
    output:
        clusters="results/{study}/phyclone/clusters.tsv",
        trace="results/{study}/phyclone/trace.pkl.gz",
        tree_nwk="results/{study}/phyclone/tree.nwk",
        tree_table="results/{study}/phyclone/tree_table.tsv"
    params:
        num_chains = config["params"]["phyclone"]["num_chains"],
        density = config["params"]["phyclone"]["density"],
        proposal = config["params"]["phyclone"]["proposal"],
        burnin = config["params"]["phyclone"]["burnin"],
        num_iters = config["params"]["phyclone"]["num_iters"],
        seed = config["params"]["phyclone"]["seed"],
        grid_size = config["params"]["phyclone"]["grid_size"]
    log:
        "logs/{study}/phyclone/phyclone.log"
    benchmark:
        "logs/{study}/phyclone/phyclone.bmk"
    conda:
        "../envs/phyclone.yaml"
    threads:
        config["resources"]["phyclone"]["threads"]
    resources:
        mem_mb = config["resources"]["phyclone"]["mem"],
        runtime = config["resources"]["phyclone"]["walltime"]
    shell:
        """
        python scripts/phyclone_cluster_formatting.py --input {input.pvi_output} --output {output.clusters} > {log}

        phyclone run -i {input.pvi_input} \
                     -c {output.clusters} \
                     -o {output.trace} \
                     --num-chains {params.num_chains} \
                     -d {params.density} \
                     --proposal {params.proposal} \
                     -n {params.num_iters} \
                     -b {params.burnin} \
                     --seed {params.seed} \
                     --grid-size {params.grid_size} \
                     --assign-loss-prob >> {log} 2>&1

         phyclone map -i {output.trace} \
                      -t {output.tree_nwk} \
                      -o {output.tree_table} >> {log} 2>&1
         """
