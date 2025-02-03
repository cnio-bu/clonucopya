rule get_vep_cache:
    output:
        directory("resources/vep/cache")
    params:
        species="homo_sapiens",
        build="GRCh38",
        release="113",
        type="merged"
    log:
        "logs/get_vep_cache/cache.log"
    resources:
        mem_mb=config["resources"]["default"]["mem"],
        disk_mb=40000,
        runtime=config["resources"]["default"]["walltime"]
    cache: "omit-software" # True
    wrapper:
        "v5.5.0/bio/vep/cache"

