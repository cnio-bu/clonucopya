rule pyclone_vi_prep:
    input:
        mutations = "results/mutation_prep/{sample}_prep.mut.tsv",
        cnvs = "results/cnv_prep/{sample}_prep.cnv.tsv"
    output:
        "results/pyclone-vi_prep/{sample}_intersect_pvi.tsv"
    log:
        "logs/pyclone-vi_prep/{sample}.log"
    benchmark:
        "logs/pyclone-vi_prep/{sample}.bmk"
    conda:
        "../envs/intersect_mutations_cnv.yaml"
    shell:
        """
        python scripts/intersect_mutations_cnv.py --mutations {input.mutations} --cnvs {input.cnvs} --output_file {output} > {log} 2>&1
    """
