def get_mutation_files(wildcards):
    samples = samples_df[samples_df['study'] == wildcards.study]['sample_id'].tolist()
    return [f"results/{wildcards.study}/mutation_prep/{sample}_prep.mut.tsv" for sample in samples]



rule report_panels:
    input:
        samplesheet = config["samplesheet"],
        intersect = "results/{study}/pyclone-vi_prep/combined_intersect_pvi.tsv",
        mut_files = get_mutation_files
    output:
        "results/{study}/report/components/report_panel.tsv"
    params:
        study = lambda wildcards: wildcards.study,
        mut_dir = lambda wildcards: f"results/{wildcards.study}/mutation_prep"
    log:
        "logs/{study}/report/report_panels.log"
    benchmark:
        "logs/{study}/report/report_panels.bmk"
    conda:
        "../envs/report_components.yaml"
    threads: 
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/build_panels.py \
             --study {params.study} \
             --samplesheet {input.samplesheet} \
             --mut_dir {params.mut_dir} \
             --intersect_combined {input.intersect} \
             --out_file {output} > {log} 2>&1
        """



rule plot_tree:
    input:
        nwk_file = "results/{study}/phyclone/tree.nwk"
    output:
        "results/{study}/report/components/clonal_tree.png"
    params:
        palette = "resources/clonucopya_palette.txt"
    log:
        "logs/{study}/report/plot_tree.log"
    benchmark:
        "logs/{study}/report/plot_tree.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/draw_clonal_tree.py \
            --nwk_file {input.nwk_file} \
            --palette {params.palette} \
            --out_file {output} > {log} 2>&1
        """



rule plot_sphere:
    input:
        tree_df = "results/{study}/phyclone/tree_table.tsv"
    output:
        directory("results/{study}/report/components/spheres_of_clones")
    params:
        palette = "resources/clonucopya_palette.txt"
    log:
        "logs/{study}/report/plot_sphere.log"
    benchmark:
        "logs/{study}/report/plot_sphere.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/draw_sphere_of_clones.py \
            --tree_df {input.tree_df} \
            --palette {params.palette} \
            --out_dir {output} > {log} 2>&1
        """



rule plot_vaf_heatmap:
    input:
        tree_df = "results/{study}/phyclone/tree_table.tsv",
        pvi_out = "results/{study}/pyclone-vi/pvi_out.tsv",
        mut_files = get_mutation_files,
        gene_alt = "results/{study}/report/components/gene_alterations.tsv"
    output:
        directory("results/{study}/report/components/vaf_heatmaps")
    params:
        mut_dir = lambda wildcards: f"results/{wildcards.study}/mutation_prep"
    log:
        "logs/{study}/report/plot_vaf_heatmap.log"
    benchmark:
        "logs/{study}/report/plot_vaf_heatmap.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/draw_vaf_heatmap.py \
            --tree_df {input.tree_df} \
            --pvi_out {input.pvi_out} \
            --mut_dir {params.mut_dir} \
            --gene_alterations {input.gene_alt} \
            --out_dir {output} > {log} 2>&1
        """



rule report_tables:
    input:
        tree_df = "results/{study}/phyclone/tree_table.tsv",
        pandrugs_dir = "results/{study}/query_pandrugs",
        mut_files = get_mutation_files
    output:
        "results/{study}/report/components/gene_alterations.tsv",
        "results/{study}/report/components/drug_prioritization.tsv"
    params:
        out_dir = directory("results/{study}/report/components"),
        mut_dir = lambda wildcards: f"results/{wildcards.study}/mutation_prep"
    log:
        "logs/{study}/report/report_tables.log"
    benchmark:
        "logs/{study}/report/report_tables.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/build_report_tables.py \
            --tree_df {input.tree_df} \
            --pandrugs_dir {input.pandrugs_dir} \
            --mut_dir {params.mut_dir} \
            --out_dir {params.out_dir} > {log} 2>&1
        """

