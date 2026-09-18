# Clonucopya Frequently Asked Questions (FAQ)

## General

**What is Clonucopya?**  
Clonucopya is an open-source, automated Snakemake workflow that integrates clonal inference with in silico drug prioritization. It reconstructs the clonal and subclonal evolutionary architecture of tumors from bulk DNA sequencing data and identifies FDA/EMA-approved and experimental candidate drugs targeting the specific genomic alterations within each clonal population.

**What problem does Clonucopya solve?**  
Existing drug prioritization tools typically work on a bulk tumor consensus genomic profile, overlooking the distinct therapeutic vulnerabilities of individual tumor clones. Clonucopya bridges this gap by combining subclonal inference (PyClone-VI + PhyClone) with clone-resolved drug prioritization (PanDrugs2), enabling the design of treatment strategies that simultaneously target all tumor clones to minimize resistance and relapse.

**Who is Clonucopya designed for?**  
The workflow is designed primarily for bioinformaticians running the analysis, but its outputs, PDF reports, tables, and figures, are intended to be interpretable by translational and clinical researchers without a bioinformatics background.

**Where can I find the code and documentation?**  
Clonucopya is available as open-source software under a GPL-3.0 license at [https://github.com/cnio-bu/clonucopya](https://github.com/cnio-bu/clonucopya), including full documentation and step-by-step tutorials.

**Who developed Clonucopya?**  
Clonucopya was developed by the Bioinformatics Unit at the Spanish National Cancer Research Centre (CNIO), Madrid, Spain. Contact: falshahrour [at] cnio [dot] es.

## Input Files

**What input files are required to run Clonucopya?**  
The minimum required inputs are:

- **SNV calls** in VCF format  
- **CNV calls** in CSV format  
- A **configuration file** in YAML format specifying run parameters  

**Can I provide additional sample metadata?**  
Yes. Complementary sample information such as sex and tumor purity can be provided via a TSV-formatted sample sheet to improve the accuracy of clonal inference.

**Can I use BAM files?**  
BAM files are optional. When provided, Clonucopya uses them to estimate variant allele frequencies (VAFs) and perform quality control on the input data.

**Can I skip upstream processing and start from PyClone-VI formatted input files?**  
Yes. Clonucopya accepts PyClone-VI formatted input files directly, allowing users to bypass the SNV/CNV intersection and preprocessing steps. This is useful for datasets with precomputed clonal inference results (use the `pvi-start` execution mode).

**Can I include indels in the analysis?**  
Yes, the inclusion of indels prior to clonal inference is optional and can be enabled through the configuration file.

**What reference genome is used?**  
By default, Clonucopya uses the GRCh38 human genome with Ensembl annotation release 103 for variant functional annotation via Ensembl VEP.

## Installation & Dependencies

**What operating systems does Clonucopya support?**  
Clonucopya supports Linux/UNIX and macOS environments.

**What version of Python is required?**  
Clonucopya requires Python v3.12.8 or compatible.

**How are software dependencies managed?**  
All software dependencies are managed through Conda (v24.11.3), making installation straightforward and reproducible.

**What are the key software components used internally?**  
The main tools integrated in the workflow are:

- **Snakemake** v9.20.0 → workflow management  
- **PyClone-VI** v0.1.6 → clonal inference and clustering  
- **PhyClone** v0.7.0 → phylogenetic tree reconstruction  
- **Ensembl VEP** v116.1 → functional variant annotation  
- **PanDrugs2 API** v2.3.0 → drug prioritization  
- **snakemake-executor-plugin-slurm** v2.6.1 → Snakemake cluster execution  
- **pandas** v2.2.2 → parse dataframes  
- **pysam** v0.24.0-0 → parse aligned read files (BAM format)  
- **requests** v2.32.3 → query PanDrugs2 database  
- **numpy** v2.1.2 → parse dataframes' columns and rows  
- **jinja2** v3.1.6 → create report from templates  
- **weasyprint** v66.0 → report rendering  
- **matplotlib** v3.10.1  
- **ete3** v3.1.3  

## Running the Workflow

**How do I configure a Clonucopya run?**  
All parameters are set in a single YAML configuration file. This file controls input paths, analysis parameters, and optional features, enabling easy setup and scalability.

**Can I run multiple samples at the same time?**  
Yes. Clonucopya's modular design processes samples independently, enabling dynamic parallel execution of multiple analyses simultaneously.

**Does Clonucopya support HPC environments?**  
Yes. In HPC environments, Snakemake submits each rule as an independent job with user-defined resource allocation. If a node fails or a system interruption occurs, only the incomplete tasks are automatically resubmitted, avoiding unnecessary recomputation.

**Can I run only part of the workflow?**  
Yes. Clonucopya supports both end-to-end and stepwise execution, giving users flexibility to run specific modules independently or resume from intermediate results.

**How do I ensure reproducible results?**  
Users can specify a fixed random seed for PyClone-VI analyses in the configuration file to guarantee reproducibility across runs.

## Methods & Analysis

**How does Clonucopya infer tumor clonal architecture?**  
Clonucopya first intersects SNV and CNV data and formats the overlapping variants as PyClone-VI input. PyClone-VI performs Bayesian clonal inference, clustering somatic mutations into clonal populations. PhyClone then reconstructs a phylogenetic tree from these clusters, assigning each cluster to a defined clone within the tumor hierarchy.

**How does drug prioritization work?**  
Somatic variants are functionally annotated with Ensembl VEP, and the annotated clone-resolved variants are queried against the PanDrugs2 API. PanDrugs2 integrates multi-source genomic and pharmacological evidence to generate a ranked, evidence-based list of therapeutic candidates for each individual clone.

**How comprehensive is the PanDrugs2 drug database?**  
PanDrugs2 v2.3.0 integrates 23 primary data sources, covering more than 74,000 drug-gene associations involving 4,642 genes and 14,659 unique compounds, including both FDA/EMA-approved drugs, drugs in clinical trials, and experimental treatments.

**What types of sequencing data are supported?**  
Clonucopya supports data from bulk DNA-seq experiments including whole-genome sequencing (WGS) and whole-exome sequencing (WES) data.

**What types of study designs are supported?**  
Clonucopya can analyze single-sample, paired-wise, multi-region, and longitudinal sequencing datasets.

## Outputs

**What outputs does Clonucopya generate?**  
Clonucopya produces three types of outputs:

1. A **self-contained PDF report** with a comprehensive summary of the clonal and therapeutic landscape  
2. **Tables** (TXT format) with clone- and drug-resolved results  
3. **Figures** (PNG format) for use in reports and presentations  

**What does the PDF report contain?**  
For each sample, the report includes:

- Per-sample quality metrics and detected variants after preprocessing  
- A phylogenetic tree and proportional visualizations of clone distribution  
- Clone-resolved VAF heatmaps showing evolutionary dynamics of somatic alterations  

**What tabular outputs are available?**  
Three tables are generated:

- **Variant-per-clone table:** somatic alterations assigned to each clonal population  
- **Drug-per-clone table:** therapeutic options prioritized for each individual clone  
- **Drug summary table:** aggregated treatment recommendations across clone groups  

## Performance & Benchmarking

**How long does a typical Clonucopya run take?**  
On a simulated dataset of 10 samples from the DREAM Challenge (paired SNV and CNV calls, no BAM files), Clonucopya generates the final report in approximately 12 hours and 50 minutes using 8 CPU cores.

## Error Handling & Troubleshooting

**Does Clonucopya validate input files before starting?**  
Yes. Prior to execution, Clonucopya performs pre-run validation checks including: sample-sheet consistency, BAM and VCF file integrity, and compliance of CNV inputs with documented format specifications.

**What happens if a sample fails during the run?**  
Failed executions preserve all intermediate files and use controlled retries to recover from transient errors without restarting completed steps. In HPC environments, only incomplete tasks are resubmitted upon failure.

**How can I debug a failed run?**  
Clonucopya generates detailed log files for each analysis step, facilitating troubleshooting and debugging of individual tasks.

## Citation

If you use Clonucopya in your research, please cite:

Sánchez-Cid G, León-Ramos C, Elena Piñeiro-Yáñez, Gómez-López G, Al-Shahrour F. *Clonucopya: A computational workflow for drug prioritization based on clonal tumor heterogeneity.*
