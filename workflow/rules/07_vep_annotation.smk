#rule vep_annotation:
#    input:
#        pvi_prep="results/pvi_vep_prep/{project}"
#    output:
#        dir=directory("results/vep_annotation/{project}")
#    params:
#        cache_dir="resources/vep/cache"
#    log:
#        "logs/vep_annotation/{project}.log"
#    benchmark:
#        "logs/vep_annotation/{project}.smk"
#    conda: "../envs/vep_annotation.yaml"
#    threads: config["resources"]["default"]["threads"]
#    resources:
#        mem_mb=config["resources"]["default"]["mem"],
#        runtime=config["resources"]["default"]["walltime"]
#    shell:
#        """
#        mkdir -p {output.dir}/{{annotations,stats}}
#
#        for clone in {input.pvi_prep}/*; do
#
#        vcf_file=$(echo $clone | sed 's/\.tsv$/.vcf/' | xargs basename -a)
#        stat_file=$(echo $clone | sed 's/\.tsv$/_summary.html/' | xargs basename -a)
#
#        vep --cache \
#        --cache_version 113 \
#        --offline \
#        --format ensembl \
#        --vcf \
#        --force_overwrite \
#        --dir_cache {params.cache_dir} \
#        --species homo_sapiens \
#        --input_file $clone \
#        --output_file {output.dir}/annotations/"$vcf_file" \
#        --stats_file {output.dir}/stats/"$stat_file" \
#        --assembly GRCh38 \
#        --verbose > {log}
#        done > {log} 2>&1
#        """


rule vep_annotation:
    input:
        pvi_prep="results/pvi_vep_prep/{project}"
    output:
        stats=directory("results/vep_annotation/{project}/stats"),
        annotations=directory("results/vep_annotation/{project}/annotations"),
#        dir=directory("results/vep_annotation/{project}")
    params:
        cache_dir="resources/vep/cache"
    log:
        "logs/vep_annotation/{project}.log"
    benchmark:
        "logs/vep_annotation/{project}.smk"
    conda: "../envs/vep_annotation.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        mkdir -p {output.annotations} {output.stats}

        for clone in {input.pvi_prep}/*; do

        vcf_file=$(echo $clone | sed 's/\.tsv$/.vcf/' | xargs basename -a)
        stat_file=$(echo $clone | sed 's/\.tsv$/_summary.html/' | xargs basename -a)

        vep --cache \
        --cache_version 113 \
        --offline \
        --format ensembl \
        --vcf \
        --force_overwrite \
        --dir_cache {params.cache_dir} \
        --species homo_sapiens \
        --input_file $clone \
        --output_file {output.annotations}/"$vcf_file" \
        --stats_file {output.stats}/"$stat_file" \
        --assembly GRCh38 \
        --verbose > {log}
        done > {log} 2>&1
        """
