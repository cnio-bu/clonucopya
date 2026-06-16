rule render_report_wf2:
    input:
        report_panel        = "results/{study}/report/components/report_panel_wf2.tsv",
        drug_summary        = "results/{study}/report/components/drug_summary.tsv",
        clonal_tree         = "results/{study}/report/components/clonal_tree.png",
        clonal_histogram    = "results/{study}/report/components/clonal_histogram.png",
        spheres_of_clones   = sphere_clones,
        vaf_heatmaps        = vaf_heatmaps,
        gene_alterations    = "results/{study}/report/components/gene_alterations.tsv",
        drug_prioritization = "results/{study}/report/components/drug_prioritization.tsv",
        template            = "resources/templates/pvi-start_study_report.template",
        logo                = "resources/templates/img/clonucopya_logo.png",
    output:
        "results/{study}/report/{study}_report_wf2.pdf"
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
