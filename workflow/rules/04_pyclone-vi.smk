rule pyclone-vi:
    input:
        "results/pre_pvi/{PATIENT}/input_{CASES}.tsv"
    output:
        fit = "results/pyclone-vi/{PATIENT}/output_{CASES}.h5",
        result = "results/pyclone-vi/{PATIENT}/output_{CASES}.tsv"
    log:
        "logs/pyclone-vi/{PATIENT}_{CASES}.log",
    benchmark:
        "logs/pyclone-vi/{PATIENT}_{CASES}.smk",
    conda:"envs/pyclone-vi.yaml"
    threads: config["resources"]["pyclone-vi"]["threads"]
    resources:
        mem = config["resources"]["pyclone-vi"]["mem"],
        walltime = config["resources"]["pyclone-vi"]["walltime"]
    shell:"""
        pyclone-vi fit -i {input} -o {output.fit} -c 40 -d beta-binomial -r 10 &&
        pyclone-vi write-results-file -i {output.fit} -o {output.result} > {log} 2>&1
    """
