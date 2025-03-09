rule query_pandrugs:
    input:
        vep_dir="results/vep_annotation/{project}/annotations"
    output:
        pandrugs_dir=directory("results/query_pandrugs/{project}")
    params:
    log:
        "logs/query_pandrugs/{project}.log"
    benchmark:
        "logs/query_pandrugs/{project}.smk"
    conda: "../envs/query_pandrugs.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        for clone in {input.vep_dir}/*.vcf; do
            clone_id=$(basename $clone | grep -o 'cluster_[0-9]\+')
            mkdir -p {output.pandrugs_dir}/"$clone_id"  
            python scripts/query_pandrugs.py \
                       --vep_vcf $clone \
                       --out_dir {output.pandrugs_dir}/"$clone_id" > {log} 2>&1
        done
        """
