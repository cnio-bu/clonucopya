rule phyclone:
    input:
        pvi_input="results/pyclone-vi_prep/{project}/pvi_input_phyclone_formatted.tsv",
        pvi_output="results/pyclone-vi/{project}/pvi_out.tsv"
    output:
        clusters="results/phyclone/{project}/clusters.tsv",
        trace="results/phyclone/{project}/trace.pkl.gz",
        tree_nwk="results/phyclone/{project}/tree.nwk",
        tree_table="results/phyclone/{project}/tree.tsv"
    log:
        "logs/phyclone/{project}/phyclone.log"
    benchmark:
        "logs/phyclone/{project}/phyclone.smk"
    conda:
        "../envs/phyclone.yaml"
    threads:
        config["resources"]["phyclone"]["threads"]
    resources:
        mem_mb = config["resources"]["phyclone"]["mem"],
        runtime = config["resources"]["phyclone"]["walltime"]
    params: 
        num_chains = config["params"]["phyclone"]["num_chains"], 
        density = config["params"]["phyclone"]["density"], 
        proposal = config["params"]["phyclone"]["proposal"],
        burnin = config["params"]["phyclone"]["burnin"],
        num_iters = config["params"]["phyclone"]["num_iters"], 
        seed = config["params"]["phyclone"]["seed"],
        grid_size = config["params"]["phyclone"]["grid_size"]
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
