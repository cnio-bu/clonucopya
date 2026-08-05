# Hands-on tutorial: Full-Set Execution Mode
In this tutorial, we present an step-by-step guide to run Clonucopya with sample data obtained from [Pyclone_VI](https://zenodo.org/record/4268826) supplementary files, specifically from Dream Challenge dataset (ICGC-TCGA DREAM Somatic Mutation Calling – Tumor Heterogeneity (SMC-Het) Challenge). For this demonstration, we will perform an study we sample P3. To launch the analysis, there are a paired SNV (data/P3.vcf) and CNA (data/P3.txt) variant calling obtained from Mutec and Battemberg, repectively. These files have the following features:

**P10.vcf**
- File type: Simulated VCF (Variant Call Format) file containing genomic variants generated for testing purposes.
- Size: 883 Kb.
- Reference genome: GRCh37/hg37.
- Samples included: One tumour sample.
- Variant content: 8192 SNVs.
- Available information: For each variant, the VCF provides genomic position, reference and alternative alleles, and sample-level metrics such as genotype, read depth, and variant allele - frequency, which are used in downstream analyses.

**P10.txt**
- File type: Simulated CNA calling (Copy Number Aberration) file containing genomic variants generated for testing purposes.
- Size: 37 Kb.
- Reference genome: GRCh37/hg37.
- Samples included: One tumour sample.
- Variant content: 75 CNAs.
- Available information: For each variant, the VCF provides chromosome, start position, end position, major copy number, ninor Copy Number, and other calling-realted information.


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
# Run this commant from clonucopya/ or adapt the paths
python workflow/scripts/process_cna_battemberg.py --input_file tutorial/test/P10.txt --output_file tutorial/test/P10_cna.tsv
```

## Settings

Once we have the sample files ready to execute Clonucopya, config file (config/config_template.yaml) and samplesheet (config/samplesheet_template.csv) must be set up. 

> We recommend to set the absolute path when a file or directory is asked to avoid confussions. However, for the sake of this tutorial, we will set relative paths to simplify the explanation.

### Config file

This time we will be running Clonucopya with default execution parameters. So we only need to rename the file (config_template.yaml to config.yaml) or a different name but you have to change at workflow/Snakefile (configfile: "../config/config.yaml"). 

The we have to set following parameters at config.yaml:
* samplesheet: "../config/samplesheet.csv"
* bam_check: False (True if we have bam files).
* bam_files: "" (we don't have indexed bam files of the study so there is no need to set the path the directory).
* just_snv: True (False if we want to keep Indels in the dataset).

### Samplesheet file

We need to fill the csv file which has the following format:
| study 	| sample_id     	| sex       	| mutations 	| cnas 	| tumour_content 	|
|-------	|-----------	|-----------	|----------	|----------	|-----------	|
| dream_P10  	| tumour 	| unknown 	| ../tutorial/test/P10.vcf        	| ../tutorial/test/P10_cna.tsv        	| 0.85         	|


## Execution

Files are formatted and settings are configured. Let's run Clonucopya full-set mode: 

```bash
# Go to workflow directory
cd workflow/

snakemake --software-deployment-method conda -j unlimited --cache
```

> First successful execution will last over 7-8 hours. VEP's reference needs to be cached. The reason to that way is because we choose to automatize caching the reference instead of ask the user for a manual downloading of the reference and place it on the correct place which can lead to incorrect executions of Clonucopya. Moreover, it reduces the dependency of packages such as docker or other related container management systems. Once the reference is cached, each sample will be processed in less than 1 hour on the following runs.


## Results
Once the workflow execution is completed, you can check the results of the study at clonucopya/workflow/results. For this use case, the output will look like this:

```
dream_P10/
├── mutation_prep/
│   └── tumour_prep.mut.tsv
│   mutation_liftover/
│   └── pvi_checked.tsv
├── pyclone-vi_prep/
│   ├── dream_P10_intersect_pvi.tsv
│   └── combined_intersect_pvi.tsv
├── pyclone-vi/
│   ├── pvi_out.h5
│   └── pvi_out.tsv
├── phyclone/
│   ├── clusters.tsv
│   ├── trace.pkl.gz
│   ├── tree.nwk
│   └── tree_table.tsv
├── mut_vep_prep
│   └── dream_P10_clone_*.tsv
├── vep_annotation/
│   ├── annotations/
│   │   └── dream_P10_clone_*.vcf
│   └── stats/ 
│       └── dream_P10_clone_*_summary.html
├── query_pandrugs/
│   └── clone_*/
│       ├── dream_P10_clone_*_computation.tsv
│       ├── {dream_P10_clone_*_vscore.vcf
│       ├── dream_P10_clone_*_gene-drug.json
│       └── dream_P10_clone_*_gene-drug.csv
└── report/
    ├── dream_P10_report.pdf
    └── components/
        ├── clonal_tree.png
        ├── drug_prioritization.tsv
        ├── drug_summary.tsv
        ├── gene_alterations.tsv
        ├── report_panel.tsv
        ├── mutation_contribution/
        │   └── tumour_mutation_contribution.png
        ├── spheres_of_clones/
        │   └── tumour_sphere_of_clones.png
        └── vaf_heatmaps/
            ├── complete/
            │   └── tumour_complete_vaf_heatmap.png
            └── sampled/
                └── tumour_sampled_vaf_heatmap.png
```

> [!IMPORTANT]
> In the directory report_results you have a file `dream_P10_report.pdf` to check how would look like the results of the execution.

# Hands-on tutorial: PVI-Start Execution Mode
In this tutorial, we present an step-by-step guide to run Clonucopya with PEACE cohort data data obtained from [Hessey, S. et al. Nature, 2026](https://doi.org/10.1038/s41586-026-10428-4) supplementary files, specifically from from a lung adenocarcinoma patient (woman, 71 years) (ID: CRUKP7127) in the PEACE/TRACERx cohort. The information regarding tumor content was obtained from [Al Bakir, M. et al. Nature, 2023](https://doi.org/10.1038/s41586-023-05729-x) supplementary files. To launch the analysis, there is a file containing all the information of the mutations required to start from Pyclone-VI. This file have the following features:

**CRUKP7127_peace.tsv**
- File type: Pyclone-VI input file containing genomic variants generated for testing purposes.
- Size: 120 Kb.
- Reference genome: GRCh37/hg37.
- Samples included:
  - CRUKP7127_SU_T1-R1 - Primary (CRUKP7127_SU_T1-R1--df37093bhhy7)
  - CRUKP7127_SU_T1-R2 - Primary (CRUKP7127_SU_T1-R2--ff616f1ahhy7)
  - CRUKP7127_SU_T1-R3 - Primary (CRUKP7127_SU_T1-R3--6f2c7730hhy7)
  - CRUKP7127_SU_T1-R4 - Primary (CRUKP7127_SU_T1-R4--9ab5e2c1hhy7)
  - CRUKP7127_BR_LN01 - Lymph node relapse (CRUKP7127_BR_LN01--2a2a4b4dhwcg)
- Variant content: 1,819 mutations.
- Available information: For each mutation, the file provides mutation id with chromosome, genomic position, reference and alternative alleles, sample id, reference counts, alternative counts, normal copy number, major copy number, minor copy number, and tumor content.


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

## Settings

Once we have the sample files ready to execute Clonucopya, config file (config/config_pvi-start_template.yaml),  alias metadata file (metadata_pvi-start_template.csv), and samplesheet (config/samplesheet_pvi-start_template.csv) must be set up. 

> We recommend to set the absolute path when a file or directory is asked to avoid confussions. However, for the sake of this tutorial, we will set relative paths to simplify the explanation.

### Config file

This time we will be running Clonucopya with default execution parameters. So we only need to rename the file (config_pvi-start_template.yaml to config_pvi-start.yaml) or a different name but you have to change at workflow/Snakefile (configfile: "config_pvi-start.yaml"). 

The we have to set following parameters at config.yaml:
* samplesheet: "../config/samplesheet_pvi-start.csv"

### Metadata alias file

To abbreviate and improve the readability, we can customize the sample names of the stud filling the csv file which has the following format:
| sample_id 	| alias     |
|-------	|-----------	|
| CRUKP7127_SU_T1-R1--df37093bhhy7  	| SU_T1-R1 	|
| CRUKP7127_SU_T1-R2--ff616f1ahhy7  	| SU_T1-R2 	|
| CRUKP7127_SU_T1-R3--6f2c7730hhy7  	| SU_T1-R3 	|
| CRUKP7127_SU_T1-R4--9ab5e2c1hhy7  	| SU_T1-R4 	|
| CRUKP7127_BR_LN01--2a2a4b4dhwcg  	| BR_LN01 	|

### Samplesheet file

We need to fill the csv file which has the following format:
| study 	| pyclone_vi     | metadata |
|-------	|-----------	|-----------	|
| CRUKP7127  	| ../tutorial/test/CRUKP7127_peace.tsv 	| ../config/metadata_pvi-start.csv |



## Execution

Files are formatted and settings are configured. Let's run Clonucopya full-set mode: 

```bash
# Go to workflow directory
cd workflow/

snakemake -s pvi-start --software-deployment-method conda -j unlimited --cache
```

> First successful execution will last over 7-8 hours. VEP's reference needs to be cached. The reason to that way is because we choose to automatize caching the reference instead of ask the user for a manual downloading of the reference and place it on the correct place which can lead to incorrect executions of Clonucopya. Moreover, it reduces the dependency of packages such as docker or other related container management systems. Once the reference is cached, each sample will be processed in less than 1 hour on the following runs.


## Results

Once the workflow execution is completed, you can check the results of the study at clonucopya/workflow/results. For this use case, the output will look like this:

```
CRUKP7127/
├── pvi-start_prep/
│   └── pvi-start_prep.tsv
├── pyclone-vi/
│   ├── pvi_out.h5
│   └── pvi_out.tsv
├── phyclone/
│   ├── clusters.tsv
│   ├── trace.pkl.gz
│   ├── tree.nwk
│   └── tree_table.tsv
├── mut_vep_prep
│   └── CRUKP7127_clone_*.tsv
├── vep_annotation/
│   ├── annotations/
│   │   └── CRUKP7127_clone_*.vcf
│   └── stats/ 
│       └── CRUKP7127_clone_*_summary.html
├── query_pandrugs/
│   └── clone_*/
│       ├── CRUKP7127_clone_*_computation.tsv
│       ├── CRUKP7127_clone_*_vscore.vcf
│       ├── CRUKP7127_clone_*_gene-drug.json
│       └── CRUKP7127_clone_*_gene-drug.csv
└── report/
    ├── CRUKP7127_report_wf2.pdf
    └── components/
        ├── clonal_tree.png
        ├── drug_prioritization.tsv
        ├── drug_summary.tsv
        ├── gene_alterations.tsv
        ├── report_panel_wf2.tsv
        ├── mutation_contribution/
        │   └──{sample_id}_mutation_contribution.png
        ├── spheres_of_clones/
        │   └── {sample_id}_sphere_of_clones.png
        └── vaf_heatmaps/
            ├── complete/
            │   └── {sample_id}_complete_vaf_heatmap.png
            └── sampled/
                └── {sample_id}_sampled_vaf_heatmap.png
```

> [!IMPORTANT]
> In the directory report_results you have a file `CRUKP7127_report_wf2.pdf` to check how would look like the results of the execution.
