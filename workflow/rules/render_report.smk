def get_spheres_files(wildcards):
    samples = samples_df[samples_df['study'] == wildcards.study]['sample_id'].tolist()
    return [
        f"results/{wildcards.study}/report/components/spheres_of_clones/{sample}_sphere_of_clones.png"
        for sample in samples
    ]


def get_heatmap_files(wildcards):
    samples = samples_df[samples_df['study'] == wildcards.study]['sample_id'].tolist()
    return [
        f"results/{wildcards.study}/report/components/vaf_heatmaps/sampled/{sample}_sampled_vaf_heatmap.png"
        for sample in samples
    ]


rule render_report:
    input:
        report_panel        = "results/{study}/report/components/report_panel.tsv",
        clonal_tree         = "results/{study}/report/components/clonal_tree.png",
        spheres_of_clones   = get_spheres_files,
        vaf_heatmaps        = get_heatmap_files,
        gene_alterations    = "results/{study}/report/components/gene_alterations.tsv",
        drug_prioritization = "results/{study}/report/components/drug_prioritization.tsv",
        template            = "resources/templates/study_report.template",
        logo                = "resources/templates/img/clonucopya_logo.png",
    output:
        "results/{study}/report/{study}_report.pdf"
    params:
        study_path = "results/{study}"
    log:
        "logs/{study}/report/render_report.log"
    benchmark:
        "logs/{study}/report/render_report.bmk"
    conda:
        "../envs/render_report.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        runtime=config["resources"]["default"]["walltime"]
    shell:
        """
        python scripts/render_study_report.py \\
             --study {params.study_path} \\
             --output_pdf {output} \\
             --template {input.template} \\
             --logo {input.logo} \\
             &> {log}
        """
