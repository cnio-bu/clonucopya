rule vep_annotation:
    input:
        pvi_prep="results/pvi_vep_prep/{sample}"
    output:
        dir=directory("results/vep_annotation/{sample}")
    params:
        cache_dir="resources"
    log:
        "logs/vep_annotation/{sample}.log"
    benchmark:
        "logs/vep_annotation/{sample}.smk"
    conda: "../envs/vep_annotation.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        for clone in {input.pvi_prep}/*; do
        out_vcf=$(echo $clone | sed 's/\.tsv$/.vcf/')
        vcf_name=$(basename $out_vcf)
        mkdir -p {output.dir}
        vep --cache \
        --cache_version 113 \
        --offline \
        --format ensembl \
        --vcf \
        --force_overwrite \
        --dir_cache {params.cache_dir} \
        --species homo_sapiens \
        --input_file $clone \
        --output_file {output.dir}/$vcf_name \
        --no_stats \
        --assembly GRCh38 \
        --verbose > {log} 2>&1
        done
        """
