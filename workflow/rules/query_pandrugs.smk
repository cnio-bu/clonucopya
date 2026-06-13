rule query_pandrugs:
    input:
        vep_dir="results/{study}/vep_annotation/annotations"
    output:
        pandrugs_dir=directory("results/{study}/query_pandrugs")
    log:
        "logs/{study}/query_pandrugs/query_pandrugs.log"
    benchmark:
        "logs/{study}/query_pandrugs/query_pandrugs.bmk"
    conda: "../envs/query_pandrugs.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        for clone in {input.vep_dir}/*.vcf; do
            clone_id=$(basename "$clone" | grep -oE 'clone_-*[0-9]+')
            mkdir -p {output.pandrugs_dir}/"$clone_id"  
            python scripts/query_pandrugs.py \
                       --vep_vcf $clone \
                       --out_dir {output.pandrugs_dir}/"$clone_id" > {log} 2>&1
        done
        """

