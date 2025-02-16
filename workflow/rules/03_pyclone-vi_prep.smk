rule pyclone_vi_prep:
    input:
#        sample_id = lambda wildcards: samples.loc[wildcards.sample, "sample_id"],
        mutations = f"results/{experiment}/mutation_prep/{{sample}}_prep.mut.tsv",
        cnvs = f"results/{experiment}/cnv_prep/{{sample}}_prep.cnv.tsv"
    output:
        f"results/{experiment}/pyclone-vi_prep/{{sample}}_intersect_pvi.tsv"
    log:
        f"logs/{experiment}/pyclone-vi_prep/{{sample}}.log"
    benchmark:
        f"logs/{experiment}/pyclone-vi_prep/{{sample}}.bmk"
    conda:
        "../envs/intersect_mutations_cnv.yaml",
    params:
        sample_id = lambda wildcards: samples.loc[wildcards.sample, "sample_id"]
    shell:
        """
        python scripts/intersect_mutations_cnv.py --sample_id {params.sample_id} --mutations {input.mutations} --cnvs {input.cnvs} --output_file {output} > {log} 2>&1
    """
