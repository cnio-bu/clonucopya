rule pvi_vcf_prep:
    input:
        pvi_df="results/pyclone-vi/{sample}_pvi_out.tsv",
        mutations=os.path.join(config["samples"]["mutations"],"{sample}.mut.vcf"),
    output:
        vep_input="results/pvi_vcf_prep/{sample}_cluster_0.tsv", 
    log:
        "logs/pvi_vcf_prep/{sample}.log"
    benchmark:
        "logs/pvi_vcf_prep/{sample}.smk"
    conda:
        "../envs/pvi_vcf_prep.yaml"
    threads: 
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem"],
        runtime = config["resources"]["default"]["walltime"]
    params:
        sample_id="{sample}",
        out_dir="results/pvi_vcf_prep"
    shell:
        """
        python scripts/vep_formatting.py --mutations {input.mutations} \
                                         --pvi_output {input.pvi_df} \
                                         --sample_id {params.sample_id} \
                                         --out_dir {params.out_dir} 2> {log}
        """
