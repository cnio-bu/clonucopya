def get_mutation_files(wildcards):
    samples = samplesheet[samplesheet['study'] == wildcards.study]['sample_id'].tolist()
    return [f"results/{wildcards.study}/mutation_prep/{sample}_prep.mut.tsv" for sample in samples]



rule report_panels:
    input:
        samplesheet = config["samplesheet"],
        pvi_input = "results/{study}/mutation_liftover/pvi_checked.tsv",
        drug_prior = "results/{study}/report/components/drug_prioritization.tsv"
    output:
        "results/{study}/report/components/report_panel.tsv"
    params:
        study = lambda wildcards: wildcards.study,
        mut_dir = lambda wildcards: f"results/{wildcards.study}/mutation_prep",
        drug_filter = config["drug_filter"]
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
        python scripts/build_panels.py \\
             --study {params.study} \\
             --samplesheet {input.samplesheet} \\
             --mut_dir {params.mut_dir} \\
             --intersect_combined {input.pvi_input} \\
             --drug_prioritization {input.drug_prior} \\
             --drug_filter {params.drug_filter} \\
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
        python scripts/draw_clonal_tree.py \\
            --nwk_file {input.nwk_file} \\
            --palette {params.palette} \\
            --out_file {output} > {log} 2>&1
        """


rule plot_sphere:
    input:
        tree_df = "results/{study}/phyclone/tree_table.tsv"
    output:
        "results/{study}/report/components/spheres_of_clones/{sample}_sphere_of_clones.png"
    params:
        palette = "resources/clonucopya_palette.txt",
        out_dir = lambda wildcards: f"results/{wildcards.study}/report/components/spheres_of_clones"
    log:
        "logs/{study}/report/plot_sphere_{sample}.log"
    benchmark:
        "logs/{study}/report/plot_sphere_{sample}.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        mkdir -p {params.out_dir}
        python scripts/draw_sphere_of_clones.py \\
            --tree_df {input.tree_df} \\
            --palette {params.palette} \\
            --out_dir {params.out_dir} > {log} 2>&1
        """


rule plot_vaf_heatmap:
    input:
        tree_df = "results/{study}/phyclone/tree_table.tsv",
        pvi_input = "results/{study}/mutation_liftover/pvi_checked.tsv",
        gene_alt = "results/{study}/report/components/gene_alterations.tsv"
    output:
        "results/{study}/report/components/vaf_heatmaps/sampled/{sample}_sampled_vaf_heatmap.png"
    params:
        out_dir = lambda wildcards: f"results/{wildcards.study}/report/components/vaf_heatmaps"
    log:
        "logs/{study}/report/plot_vaf_heatmap_{sample}.log"
    benchmark:
        "logs/{study}/report/plot_vaf_heatmap_{sample}.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        mkdir -p {params.out_dir}/sampled
        python scripts/draw_vaf_heatmap.py \\
            --tree_df {input.tree_df} \\
            --pvi_input {input.pvi_input} \\
            --gene_alterations {input.gene_alt} \\
            --out_dir {params.out_dir} > {log} 2>&1
        """



rule report_tables:
    input:
        pvi_input = "results/{study}/mutation_liftover/pvi_checked.tsv",
        pandrugs_dir = "results/{study}/query_pandrugs"
    output:
        "results/{study}/report/components/gene_alterations.tsv",
        "results/{study}/report/components/drug_prioritization.tsv",
        "results/{study}/report/components/drug_summary.tsv"
    params:
        out_dir = directory("results/{study}/report/components")
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
        python scripts/build_report_tables.py \\
            --pvi_input {input.pvi_input} \\
            --pandrugs_dir {input.pandrugs_dir} \\
            --out_dir {params.out_dir} > {log} 2>&1
        """



rule mutation_contribution:
    input:
        tree_df = "results/{study}/phyclone/tree_table.tsv"
    output:
        "results/{study}/report/components/mutation_contribution/{sample}_mutation_contribution.png"
    params:
        palette = "resources/clonucopya_palette.txt",
        out_dir = lambda wildcards: f"results/{wildcards.study}/report/components/mutation_contribution"
    log:
        "logs/{study}/report/{sample}_plot_histogram.log"
    benchmark:
        "logs/{study}/report/{sample}_plot_histogram.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        mkdir -p {params.out_dir}
        python scripts/draw_clone_histogram.py \\
            --tree_df {input.tree_df} \\
            --palette {params.palette} \\
            --out_dir {params.out_dir} > {log} 2>&1
        """



rule study_composition:
    input:
        tree_df = "results/{study}/phyclone/tree_table.tsv"
    output:
        "results/{study}/report/components/{study}_clonal_composition.png"
    params:
        study = lambda wildcards: wildcards.study,
        palette = "resources/clonucopya_palette.txt",
        out_dir = lambda wildcards: f"results/{wildcards.study}/report/components"
    log:
        "logs/{study}/report/plot_study_composition.log"
    benchmark:
        "logs/{study}/report/plot_study_composition.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        mkdir -p {params.out_dir}
        python scripts/draw_study_clonal_composition.py \\
            --tree_df {input.tree_df} \\
            --study {params.study} \\
            --palette {params.palette} \\
            --out_dir {params.out_dir} > {log} 2>&1
        """


######################### TABLES FOR ALTERNATIVE WORKFLOW ###############################

rule report_panels_wf2:
    input:
        samplesheet = config["samplesheet"],
        phy_out = "results/{study}/phyclone/tree_table.tsv",
        pvi_input = "results/{study}/mutation_liftover/pvi_checked.tsv",
        drug_prior = "results/{study}/report/components/drug_prioritization.tsv"
    output:
        "results/{study}/report/components/report_panel_wf2.tsv"
    params:
        study = lambda wildcards: wildcards.study,
        drug_filter = config["drug_filter"]
    log:
        "logs/{study}/report/report_panels_wf2.log"
    benchmark:
        "logs/{study}/report/report_panels_wf2.bmk"
    conda:
        "../envs/report_components.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/build_panels_wf2.py --study {params.study} \\
                                           --samplesheet {input.samplesheet} \\
                                           --phy_out {input.phy_out} \\
                                           --pvi_input {input.pvi_input} \\
                                           --drug_prioritization {input.drug_prior} \\
                                           --drug_filter {params.drug_filter} \\
                                           --out_file {output} > {log} 2>&1
        """


