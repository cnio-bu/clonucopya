rule render_report:
    input:
        study_path = "results/{study}"
    output:
        "results/{study}/report/{study}_report.pdf"
    params:
        template = "resources/templates/study_report.template",
        logo = "resources/templates/img/clonucopya_logo.png"
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
        python scripts/render_study_report.py \
             --study {input.study_path} \
             --output_pdf {output} \
             --template {params.template} \
             --logo {params.logo} \
        """

