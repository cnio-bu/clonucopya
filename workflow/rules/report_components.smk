rule report_panels:
    input:
        samplesheet = config["samplesheet"],
        mut_dir = "results/mutation_prep/{project}",
        intersect = "results/pyclone-vi_prep/{project}/combined_intersect_pvi.tsv"
    output:
        "results/report/{project}/components/report_panel.tsv"
    params:
        project = lambda wildcards: wildcards.project
    log:
        "logs/report/{project}/report_panels.log"
    benchmark:
        "logs/report/{project}/report_panels.bmk"
    conda:
        "../envs/report_componets.yaml"
    threads: 
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/build_panels.py \
             --project {params.project} \
             --samplesheet {input.samplesheet} \
             --mut_dir {input.mut_dir} \
             --intersect_combined {input.intersect} \
             --out_file {output} > {log} 2>&1
        """



rule plot_tree:
    input:
        nwk_file = "results/phyclone/{project}/tree.nwk"
    output:
        "results/report/{project}/components/clonal_tree.png"
    params:
        palette = "resources/clonucopya_palette.txt"
    log:
        "logs/report/{project}/plot_tree.log"
    benchmark:
        "logs/report/{project}/plot_tree.bmk"
    conda:
        "../envs/report_componets.yaml"
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
        tree_df = "results/phyclone/{project}/tree.tsv"
    output:
        directory("results/report/{project}/components/spheres_of_clones")
    params:
        palette = "resources/clonucopya_palette.txt"
    log:
        "logs/report/{project}/plot_sphere.log"
    benchmark:
        "logs/report/{project}/plot_sphere.bmk"
    conda:
        "../envs/report_componets.yaml"
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
        tree_df = "results/phyclone/{project}/tree.tsv",
        pvi_out = "results/pyclone-vi/{project}/pvi_out.tsv",
        mut_dir = "results/mutation_prep/{project}"
    output:
        directory("results/report/{project}/components/vaf_heatmaps")
    log:
        "logs/report/{project}/plot_vaf_heatmap.log"
    benchmark:
        "logs/report/{project}/plot_vaf_heatmap.bmk"
    conda:
        "../envs/report_componets.yaml"
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
            --mut_dir {input.mut_dir} \
            --out_dir {output} > {log} 2>&1
        """



rule report_tables:
    input:
        tree_df = "results/phyclone/{project}/tree.tsv",
        pandrugs_dir = "results/query_pandrugs/{project}",
        mut_dir = "results/mutation_prep/{project}"
    output:
        "results/report/{project}/components/gene_alterations.tsv",
        "results/report/{project}/components/drug_priorization.tsv"
    params:
        out_dir = directory("results/report/{project}/components")
    log:
        "logs/report/{project}/report_tables.log"
    benchmark:
        "logs/report/{project}/report_tables.bmk"
    conda:
        "../envs/report_componets.yaml"
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
            --mut_dir {input.mut_dir} \
            --out_dir {params.out_dir} > {log} 2>&1
        """

