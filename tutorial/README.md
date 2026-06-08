# Hands-on tutorial
In this tutorial, we present an step-by-step guide to run Clonucopya with sample data obtained from Pyclone-VI supplementary files, specifically from Dream Challenge dataset (ICGC-TCGA DREAM Somatic Mutation Calling – Tumor Heterogeneity (SMC-Het) Challenge). For this demonstration, we will perform an study we sample P3. To lauch the analysis, there are a paired SNV (data/P3.vcf) and CNA (data/P3.txt) variant calling obtained from Mutec and Battemberg, repectively. 

## Installation

```bash
# Clone Clonucopya repository
git clone https://github.com/cnio-bu/clonucopya.git

# Create a conda environment
conda create -n clonucopya

# Activate the environment
conda activate clonucopya

# Intall Required Packages
mamba install snakemake apptainer snakemake-executor-plugin-slurm pandas
```

## CNA preprocessing

The CNA file is straight from Battemberg output, so a light preprocessing is required to meet the CNA file format of Clonucopya workflow: 


| Chrom 	| Start     	| End       	| major_cn 	| minor_cn 	| normal_cn 	|
|-------	|-----------	|-----------	|----------	|----------	|-----------	|
| chr1  	| 109367944 	| 109371874 	| 1        	| 0        	| 2         	|


Clonucopya provides a few script to format CNA calling output for caller such as Battemberg, ascat3 or Facets. If your file does not meet none of those format, you can custom the script to suit your case. In this example, we will use the Battemberg preprocessing script `workflow/scripts/process_cna_battemberg.py` as follows:

```python
python workflow/scripts/process_cna_battemberg.py --input_file tutorial/test/P3.txt --output_file tutorial/test/P3_cna.tsv
```

## Settings

Once we have the sample files ready to execute clonucopya, config file (config/config_template.yaml) and samplesheet (config/samplesheet_template.csv) must be set up. 

> We recommend to set the absolute path when a file or directory file is asked to avoid confussions. However, for the sake of this tutorial, we will set relative paths to simplify the explanation.

### Config file

This time we will be running Clonucopya with default execution parameters. So we only need to rename the file (config_template.yaml to config.yaml) or a different name but you have to change at workflow/Snakefile (configfile: "../config/config.yaml"). 

The we have to set following parameters at config.yaml:
* samplesheet: "./samplesheet.csv"
* bam_check: False (True if we have bam files).
* bam_files: "" (we don't have indexed bam files of the study so there is no need to set the path the directory).
* just_snv: True (False if we want to keep Indels in the dataset).

### Samplesheet file

We need to fill the csv file which has the following format:
| study 	| sample_id     	| sex       	| mutations 	| cnas 	| tumour_content 	|
|-------	|-----------	|-----------	|----------	|----------	|-----------	|
| dream_P3  	| tumour 	| female 	| ../tutorial/test/P3.vcf        	| tutorial/test/P3_cna.tsv        	| 0.85         	|


## Execution

Files are formatted and settings are configured. Let's run Clonucopya full-set mode: 

```bash
# Go to workflow directory
cd workflow/

snakemake --software-deployment-method conda -j unlimited --cache
```

> First successful execution will last over 7-8 hours. VEP's reference needs to be cached. The reason to that way is because we choose to automatize caching the reference instead of ask the user for a manual downloading of the reference and place it on the correct place which can lead to incorrect executions of Clonucopya. Moreover, it reduces the dependency of packages such as docker or other related container management systems. Once the reference is cached, each sample will be processed in less than 1 hour on the following runs.


## Results

For this use case, the output will look like this:

```
dream_P3/
├── pyclone-vi/
│   ├── pvi_out.h5
│   └── pvi_out.tsv
├── phyclone/
│   ├── clusters.tsv
│   ├── trace.pkl.gz
│   ├── tree.nwk
│   └── tree_table.tsv
├── pvi_vep_prep
│   └── dream_P3_clone_*.tsv
├── vep_annotation/
│   ├── annotations/
│   │   └── dream_P3_clone_*.vcf
│   └── stats/ 
│       └── dream_P3_clone_*_summary.html
├── query_pandrugs/
│   └── clone_*/
│       ├── dream_P3_clone_*_computation.tsv
│       ├── {dream_P3_clone_*_vscore.vcf
│       ├── dream_P3_clone_*_gene-drug.json
│       └── dream_P3_clone_*_gene-drug.csv
└── report/
    ├── dream_P3_report.pdf
    └── components/
        ├── clonal_tree.png
        ├── clonal_histogram.png
        ├── drug_prioritization_wf2.tsv
        ├── drug_summary_wf2.tsv
        ├── gene_alterations_wf2.tsv
        ├── report_panel.tsv
        ├── spheres_of_clones/
        │   └── tumour_sphere_of_clones.png
        └── vaf_heatmaps/
            ├── complete/
            │   └── tumour_complete_vaf_heatmap.png
            └── sampled/
                └── tumour_sampled_vaf_heatmap.png
```


