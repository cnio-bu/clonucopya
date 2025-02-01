rule vep_annotation:
    input:
        pvi_prep="results/pvi_vep_prep/{sample}_cluster_0.tsv"
    output:
#        vep_vcf="results/vep_annotation/{sample}/{sample}_{cluster}.vcf",
#        summary="results/vep_annotation/{sample}/{sample}_{cluster}.vcf_summary.html"
#results/vep_annotation/{sample}/{sample}_cluster_0.vcf
        vep_vcf="results/vep_annotation/{sample}_cluster_0.vcf",
        summary="results/vep_annotation/{sample}_cluster_0.vcf_summary.html"
    log:
        "logs/vep_annotation/{sample}.log"
    benchmark:
        "logs/vep_annotation/{sample}.smk"
    conda: "../envs/vep_annotation.yaml"
    threads: config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    params:
        cache_dir="resources"
    shell:
        """
        vep --cache \
        --cache_version 113 \
        --offline \
        --format ensembl \
        --vcf \
        --force_overwrite \
        --dir_cache {params.cache_dir} \
        --species homo_sapiens \
        --input_file {input.pvi_prep} \
        --output_file {output.vep_vcf} \
        --assembly GRCh38 \
        --verbose > {log} 2>&1

        """
