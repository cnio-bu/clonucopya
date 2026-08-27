rule vep_annotation:
    input:
        vep_prep="results/{study}/mut_vep_prep",
        cache_dir="resources/vep/cache"
    output:
        stats=directory("results/{study}/vep_annotation/stats"),
        annotations=directory("results/{study}/vep_annotation/annotations"),
    params:
    log:
        "logs/{study}/vep_annotation/annotation.log"
    benchmark:
        "logs/{study}/vep_annotation/annotation.bmk"
    conda: "../envs/vep_annotation.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        r"""
        mkdir -p {output.annotations} {output.stats}
    
        for clone in {input.vep_prep}/*; do
    
            vcf_file=$(echo "$clone" | sed 's/.tsv$/.vcf/' | xargs basename -a)
            stat_file=$(echo "$clone" | sed 's/.tsv$/_summary.html/' | xargs basename -a)
    
            vep --cache \
                --cache_version 116 \
                --offline \
                --format ensembl \
                --vcf \
                --force_overwrite \
                --dir_cache {input.cache_dir} \
                --species homo_sapiens \
                --input_file "$clone" \
                --output_file "{output.annotations}/$vcf_file" \
                --stats_file "{output.stats}/$stat_file" \
                --assembly GRCh38 \
                --verbose
        done > {log} 2>&1
        """
