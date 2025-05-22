rule pvi_intesersect:
    input:
        mutations = "results/mutation_prep/{project}/{sample}_prep.mut.tsv",
        cnvs = lambda wildcards: samples.loc[wildcards.sample, "cnvs"]
    output:
        "results/pyclone-vi_prep/{project}/{sample}_intersect_pvi.tsv"
    log:
        "logs/pyclone-vi_prep/{project}/{sample}.log"
    benchmark:
        "logs/pyclone-vi_prep/{project}/{sample}.bmk"
    conda:
        "../envs/intersect_mutations_cnv.yaml",
    params:
        sample_id = lambda wildcards: wildcards.sample
    shell:
        """
        python scripts/intersect_mutations_cnv.py \
              --sample_id {params.sample_id} \
              --mutations {input.mutations} \
              --cnvs {input.cnvs} \
              --output_file {output} > {log} 2>&1
    """


rule concat_and_purity_pvi:
    input:
        lambda wildcards: expand("results/pyclone-vi_prep/{project}/{sample}_intersect_pvi.tsv", 
                                 sample=samples_df[samples_df['project'] == wildcards.project]['sample_id'],
                                 project=wildcards.project)
    output:
        pvi = "results/pyclone-vi_prep/{project}/combined_intersect_pvi.tsv",
        phyclone_prep = "results/pyclone-vi_prep/{project}/pvi_input_phyclone_formatted.tsv"
    params:
        samplesheet = config["samplesheet"]
    log:
        "logs/pyclone-vi_prep/{project}/concat.log"
    conda:
        "../envs/intersect_mutations_cnv.yaml"
    shell:
        """
        python -c "
import pandas as pd
import sys

# List of pvi prep samples
input_files = {input!r}
pvi_preps = [pd.read_csv(file, sep='\t') for file in input_files]

# Concatenate
combined_pvi = pd.concat(pvi_preps, ignore_index=True)

# Drop duplicates

combined_pvi_dedup = combined_pvi.drop_duplicates()

combined_pvi_dedup.to_csv('{output.phyclone_prep}', sep='\t', index=False)

samplesheet = pd.read_csv('{params.samplesheet}')

tumour_content_dict = dict(zip(samplesheet['sample_id'], samplesheet['tumour_content']))

combined_pvi_dedup['tumour_content'] = combined_pvi_dedup['sample_id'].map(tumour_content_dict)

combined_pvi_dedup.to_csv('{output.pvi}', sep='\t', index=False)
" > {log} 2>&1
        """




