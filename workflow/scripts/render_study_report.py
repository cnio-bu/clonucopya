import os
import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, Template
from datetime import datetime
import weasyprint


def safe_float_format(value, decimals=2):
    """Convert value to float"""
    try:
        if pd.isna(value) or value is None:
            return "N/A"
        float_val = float(value)
        return f"{float_val:.{decimals}f}"
    except (ValueError, TypeError):
        return str(value) if value is not None else "N/A"

def safe_int_format(value):
    """Convert a value to integer"""
    try:
        if pd.isna(value) or value is None:
            return "N/A"
        return str(int(float(value)))
    except (ValueError, TypeError):
        return str(value) if value is not None else "N/A"

def create_jinja_env(template_dir="."):
    """Create jinja2 environment with custom filters"""
    env = Environment(loader=FileSystemLoader(template_dir))
    env.filters['safe_float'] = safe_float_format
    env.filters['safe_int'] = safe_int_format
    return env

def parse_clones(s):
    if pd.isna(s) or s == "":
        return []
    return [int(x.strip()) for x in str(s).split(",") if x.strip() != ""]

def get_clonucopya_version():
    try:
        return subprocess.check_output(
            ["git", "describe", "--tag"],
            text=True, stderr=subprocess.DEVNULL
        ).strip()
    except subprocess.CalledProcessError:
        return "no information available"

def prepare_report_data(study_path, drug_filter):
    """
    Prepare data for the template injection with absolute paths
    """
    study_name = os.path.basename(study_path)
    components_path = Path(study_path) / "report" / "components"
    

    # Load report panel
    rp1 = components_path / "report_panel.tsv"
    rp2 = components_path / "report_panel_wf2.tsv"

    if rp1.is_file():
        rp_path = rp1
    elif rp2.is_file():
        rp_path = rp2
    else:
        raise FileNotFoundError("No se encontraron report_panel*.tsv")

    samples_df = pd.read_csv(rp_path, sep='\t')

    # Drug Summary
    drug_sum_path = components_path / "drug_summary.tsv"
    drug_sum_top = None
    if drug_sum_path.exists():
        try:
            drug_sum = pd.read_csv(drug_sum_path, sep='\t')
            if not drug_sum.empty:
                if drug_filter == 'clinical':
                     drug_sum = drug_sum[(drug_sum['Status'] != 'EXPERIMENTAL') & (drug_sum['Interaction_Type'] != 'PATHWAY_MEMBER')]

                # Create support columns to sort by number of target clones
                drug_sum["Target_Clones_list"] = drug_sum["clone_list"].apply(parse_clones)
                drug_sum["n_clones"] = drug_sum["Target_Clones_list"].apply(lambda xs: len(set(xs)))

                # Sort interaction_type
                interaction_order = ["DIRECT_TARGET", "BIOMARKER", "PATHWAY_MEMBER"]
                drug_sum["Interaction_Type"] = pd.Categorical(
                            drug_sum["Interaction_Type"],
                            categories=interaction_order,
                            ordered=True
                )
                
                # FILTER SENSITIVITY  
                sensitivity = (
                    drug_sum[drug_sum["max_dScore"] > 0]
                    .sort_values(["n_clones", "Status", 'Interaction_Type', "max_dScore"], ascending=[False, True, True, False])
                    .head(25)
                    .copy()
                )

                # FILTER RESISTANCE
                resistance =(
                    drug_sum[drug_sum["max_dScore"] < 0]
                    .sort_values(["n_clones", "Status", 'Interaction_Type', "max_dScore"], ascending=[False, True, True, False])
                    .head(25)
                    .copy()
                )
                
                drug_sum_top = pd.concat([sensitivity, resistance], axis=0)
                drug_sum_top.drop(columns = ["n_clones", "Target_Clones_list", "Interaction_Type", "clone_list"], axis=1, inplace=True)
                
        except:
            drug_sum_top = None
    
    # Gene alterations
    gene_alterations_path = components_path / "gene_alterations.tsv"
    gene_alt = None
    
    if gene_alterations_path.exists():
        try:
            gene_alt = pd.read_csv(gene_alterations_path, sep='\t')
            if gene_alt.empty:
                gene_alt = None
            else:
                gene_alt = gene_alt[gene_alt["Clone"] != -1]
                gene_alt = gene_alt[(gene_alt['Impact'] == 'MODERATE') 
                | (gene_alt['Impact'] == 'HIGH')]
                # Sort by clone and impact
                impact_order = pd.CategoricalDtype(categories=["HIGH", "MODERATE"], ordered=True)
                gene_alt["Impact"] = gene_alt["Impact"].str.upper().astype(impact_order)
                gene_alterations_sorted = gene_alt.sort_values(by=["Clone", "Impact"], ascending=[True, True])
                # Select top 10 alterations per clone
                gene_alterations_top = gene_alterations_sorted.groupby("Clone", sort=False).head(10)
        except:
            gene_alt = None
        

    # Drug prioritization
    drug_prioritization_path = components_path / "drug_prioritization.tsv"
    drug_hits = None
    if drug_prioritization_path.exists():
        try:
            drug_hits = pd.read_csv(drug_prioritization_path, sep='\t')
            if drug_hits.empty:
                drug_hits = None
            else:
                drug_hits = drug_hits[drug_hits["Clone"] != -1]
                drug_hits["dScore"] = pd.to_numeric(drug_hits["dScore"], errors="coerce")

                # Sort Status and Intereaction type column
                status_order = ["APPROVED", "CLINICAL_TRIALS", "EXPERIMENTAL"]
                drug_hits["Status"] = pd.Categorical(drug_hits["Status"], categories=status_order, ordered=True)

                interaction_type_order = ["DIRECT_TARGET", "BIOMARKER", "PATHWAY_MEMBER"]
                drug_hits["Interaction Type"] = pd.Categorical(drug_hits["Interaction Type"], categories=interaction_type_order, ordered=True)

                if drug_filter == 'clinical':
                    drug_hits = drug_hits[(drug_hits['Status'] != 'EXPERIMENTAL') & (drug_hits['Interaction Type'] != 'PATHWAY_MEMBER')]
                
                # Grouping key
                group_keys = ["Clone", "Mutation ID", "Gene Symbol", "VAF"]
                
                drugs_df_sorted = drug_hits.sort_values(
                    by=group_keys + ["Status", "dScore", "gScore","Interaction Type"],
                    ascending=[True, True, True, True, True, False, False, True]
                )
                
                # Top 3 drugs for "each mutation of the clone" (an specific mutation may be targeted by one or more drugs)
                drugs_df_compact = (
                    drugs_df_sorted
                    .groupby(group_keys, sort=False)
                    .head(3)
                    .reset_index(drop=True)
                )
                
                # Colapse repeated key columns for readability
                collapse_cols = ["Clone", "Mutation ID", "Gene Symbol", "VAF"]
                drugs_df_compact[collapse_cols] = drugs_df_compact[collapse_cols].astype("string")
                
                within_group_index = drugs_df_compact.groupby(collapse_cols).cumcount()
                is_dup = within_group_index > 0
                
                drugs_df_compact.loc[is_dup, collapse_cols] = ""
                
        except Exception as e:
            print(f"[ERROR] fail to process drug prioritization results: {repr(e)}")
            drugs_df_compact = None
    
    # Absolute paths to images
    clonal_tree_image = os.path.abspath(str(components_path / "clonal_tree.png"))
    study_clonal_composition = os.path.abspath(str(components_path / f"{study_name}_clonal_composition.png"))
    clonal_histogram_images = {}
    clonal_proportions_images = {}
    clone_alterations_images = {}

    for _, sample in samples_df.iterrows():
        sample_id = sample['sample_id']

        # Sphere of clones
        sphere_path = components_path / "spheres_of_clones" / f"{sample_id}_sphere_of_clones.png"
        if sphere_path.exists():
            clonal_proportions_images[sample_id] = os.path.abspath(str(sphere_path))
        else:
            print(f"[WARN] sphere of clones not found for {sample_id}: {sphere_path}")

            # Clonal histogram
        histogram_path = components_path / "mutation_contribution" / f"{sample_id}_mutation_contribution.png"
        if histogram_path.exists():
            clonal_histogram_images[sample_id] = os.path.abspath(str(histogram_path))
        else:
            print(f"[WARN] clonal histogram not found for {sample_id}: {histogram_path}")

        # VAF heatmap
        heatmap_path = components_path / "vaf_heatmaps" / "sampled" / f"{sample_id}_sampled_vaf_heatmap.png"
        if heatmap_path.exists():
            clone_alterations_images[sample_id] = os.path.abspath(str(heatmap_path))
        else:
            print(f"[WARN] VAF heatmap not found for {sample_id}: {heatmap_path}")
    
    
    return {
        'samples_df': samples_df,
        'study_clonal_composition': study_clonal_composition,
        'drug_summary': drug_sum_top,
        'gene_alterations_df': gene_alterations_top,
        'drug_prioritization_df': drugs_df_compact,
        'clonal_tree_image': clonal_tree_image,
        'clonal_histogram_images': clonal_histogram_images,
        'clonal_proportions_images': clonal_proportions_images,
        'clone_alterations_images': clone_alterations_images
    }


def render_report_to_pdf(study_path, drug_filter, output_path, template_path="template.html", logo_path=None):
    """
    Render template to PDF directly
    """

    # Set path to source files

    study_name = os.path.basename(study_path)
    components_path = Path(study_path) / "report" / "components"
    panels_path = components_path / "report_panel*.tsv"
    study_composition_path = components_path / "study_name}_clonal_composition.png"
    drug_summary_path = components_path / "drug_summary.tsv"
    clonal_tree_path = components_path / "clonal_tree.png"
    clonal_histogram_path = components_path / "mutation_contribution" / "{sample_id}_mutation_contribution.png"
    spheres_path = components_path / "spheres_of_clones" / "{sample_id}_sphere_of_clones.png"
    heatmaps_path = components_path / "vaf_heatmaps" / "sampled" / "{sample_id}_sampled_vaf_heatmap.png"
    gene_alterations_path = components_path / "gene_alterations.tsv"
    drug_prioritization_path = components_path / "drug_prioritization.tsv"
    
    # GET EXECUTION INFORMATION
    clonucopya_version = get_clonucopya_version().split('-')[0]
    report_timestamp = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
    
    # Format data
    data = prepare_report_data(study_path, drug_filter)
    
    # Prepare logo as base64
    logo_data_uri = None
    if logo_path and os.path.exists(logo_path):
        import base64
        with open(logo_path, 'rb') as f:
            logo_bytes = f.read()
        
        # Detect format
        if logo_path.lower().endswith('.png'):
            mime_type = 'image/png'
        elif logo_path.lower().endswith(('.jpg', '.jpeg')):
            mime_type = 'image/jpeg'
        else:
            mime_type = 'image/png'
        
        logo_base64 = base64.b64encode(logo_bytes).decode('utf-8')
        logo_data_uri = f"data:{mime_type};base64,{logo_base64}"
    
    # Initialize jinja environment
    template_dir = os.path.dirname(os.path.abspath(template_path)) or "."
    template_filename = os.path.basename(template_path)
    
    env = create_jinja_env(template_dir)
    
    try:
        template = env.get_template(template_filename)
        print(f"Template correctly loaded: {template_path}")
    except Exception as e:
        print(f"Error loading template: {e}")
        return False
    
    # Template's data
    template_data = {
        'system_name': 'Clonucopya',
        'report_title': 'Clonal Evolution & Treatment Report',
        'report_footer_line': f"Clonucopya v{clonucopya_version} | Generated: {report_timestamp}",
        'study_id': study_name,
        'study_description': f"""
        <p>Overview of study {study_name}'s key statistics by sample, including the number of sample's name, sex, number of mutations such as Single Nucleotide Variation (SNVs) or small Indels, Copy Number Variations (CNVs), the number of intersections between mutations and CNVs, total number of drugs, and the Best Therapeutic Candidates (BTCs).</p>
        <p>This report was generated in {drug_filter} mode. In clinical mode, results are filtered more stringently, excluding drugs with experimental status and pathway member interaction type. In discovery mode, no filters are applied, and all drug hits identified are displayed regardless of their experimental status or interaction type.</p>

        <p>The source files are available at: {panels_path}.</p>
""",
        'logo_path': logo_data_uri,
    
        'samples_df': data['samples_df'],
        'study_clonal_composition': data['study_clonal_composition'],
        'drug_summary': data['drug_summary'],
        'gene_alterations_df': data['gene_alterations_df'],
        'drug_prioritization_df': data['drug_prioritization_df'],
        'clonal_histogram_images': data['clonal_histogram_images'],
        'clonal_proportions_images': data['clonal_proportions_images'],
        'clone_alterations_images': data['clone_alterations_images'],


            'study_composition_description': f"""
            <p> This is a rose of nightgale plot. It represents clonal and mutational composition at study level. Each bar represents a clone and its genetics features: the width of the bar is the clonal prevalence (weight of the clone in the clonal structure of the study), and the height is related to the number of mutations contributing to the clone identity.</p>
            <p>It is important to remark that there is no direct correlation between the number of mutations and the clonal prevalence. The impact of each mutation leads to different grades of importance in the representacion at the clonal structure of the study.</p>
            <p>Note that sometimes the total sum of the clonal prevalence proportions it not 1 (100%). This missing percentage of clonal prevelence represents the unknown or undetermined part of the clonal phylogeny due to diferent reasons: Contamination of the sample with normal (nontumor) cells, which PhyClone already correct for using the tumor purity parameter, statistical uncertainty or rounding in the estimate of prevalences, or other causes.</p>
            <ul>           
            <p>The source files are available at: {study_composition_path}.</p>
            """,
        
            'image': data['clonal_proportions_images'],

        'subclonal_tree': {
            'title': 'Tumor Phylogeny',
            'description': f"""
                <p>Based on the variant allele frequency (VAF) and changes in copy number, the evolutionary hierarchy of the tumor has been inferred. The phylogenetic tree illustrates the relationships among the different cellular subpopulations (clones).</p>
                <p>The source files are available at: {clonal_tree_path}.</p>
            """,
            'image': data['clonal_tree_image'],
        },
         'drug_summary_description': f"""
            <p>The Drug Summary provides a high-level overview of the therapeutic candidates identified across the entire study. For each drug, the table reports the number of genetic alterations supporting its prioritization, its regulatory approval status (APPROVED, CLINICAL_TRIALS, or EXPERIMENTAL), and the type of interaction with the affected genes (DIRECT TARGET, BIOMARKER, or PATHWAY_MEMBER).</p>
            <p>This report was generated in {drug_filter} mode. In clinical mode, results are filtered more stringently, excluding drugs with experimental status and pathway member interaction type. In discovery mode, no filters are applied, and all drug hits identified are displayed regardless of their experimental status or interaction type.</p>
            <p>This table summarizes up to 25 top-ranked drug candidates derived from the mutational landscape of all samples included in the study. A detailed per-clone breakdown is available in the Drug Prioritization section.</p>
            <p>The response of each drug is evaluated using the DScore (range -1 to +1):</p>
            <ul>
            <li>Sensitivity (DScore > 0): Positive values indicate that the mutations confer sensitivity to the drug (optimal candidates with a DScore > 0.7 are highlighted in bold).</li>
            <li>Resistance (DScore < 0): Negative values suggest that the clone possesses mutations associated with therapeutic resistance.</li>
            </ul>
            <p>Drug-Gene Interactions:</p>
            <ul>
            <li>DT (Direct Target): The drug acts directly on the mutated protein.</li>
            <li>BM (Biomarker): The alteration acts as a biomarker predictive of response.</li>
            <li>PM (Pathway Member): The mutation affects the signaling pathway targeted by the drug.</li>
            </ul>            
            <p>The table below provides a simplified overview of the drugs targeting the affected clones and the specific genes. Drugs are prioritizited using the following criteria: (1) clone coverage, (2) maximum dScore of the drug regarding the targeted genes (dScore), and (3) interaction type importance (DT>BM>PM).</p>
            <p>The source files are available at: {drug_summary_path}.</p>
            """,

        'comparison_section': {
            'title': 'Sample-Level Clonal Structure Analysis',
            'description': f"""
            <p>Analysis of clonal evolution patterns and proportions across samples.</p>
            <p>The source files are available at: {clonal_histogram_path}.</p>
            """,
            'clonal_histogram_images': data['clonal_histogram_images'],

            'clonal_proportions': {
                'title': 'Clonal Proportions',
                'description': f"""
                    <p>Based on the phylogeny results obtained by PhyClone, the spheres representations show the distribution and percentage of each clone across the different biopsies or time points, allowing us to visualize which clones dominate the tumor tissue. </p>
                    <p>The source files are available at: {spheres_path}.</p>
                """
            },
            'clone_alterations': {
                'title': 'Clonal Alterations',
                'description': f"""
                    <p>Heat maps illustrate the spatial or temporal evolution of mutations within clones. Observing how the frequency of a mutation varies across different samples allows us to distinguish early “trunk” events (present in all cells and responsible for the origin of the tumor) from later “branch” events. Identifying these branches is essential for understanding intratumoral heterogeneity and detecting emerging clones that may be driving disease progression or treatment resistance. Only variants with moderate or high predictive impact are shown.</p>
                    <p>The source files are available at: {heatmaps_path}.</p>
                """
            }
        },
        'alteration_section': {
            'title': 'Pharmacogenomics',
            'gene_alterations': {
                'title': 'Gene Alterations',
                'description': f"""
                <p>Both SNVs and small indels (if selected at the beginning of the workflow) used in clone inference are displayed with the most relevant information at the clone and sample level. To sum up the results, this section only shows top 10 mutations per clone with moderate or high impact.</p>
                <p> Mutations are sorted using the following criteria:</p>
                <ul>
                    <li>Clone.</li>
                    <li>Mutation positon.</li>
                    <li>Impact: High and Moderate.</li>
                </ul>
                <p>The source files are available at: {gene_alterations_path}.</p>
                """
            },
            'drug_prioritization': {
                'title': 'Drug Prioritization',
                'description': f"""
                    <p>The small variant analysis performed by PanDrugs2 displays only mutations classified as clinically relevant.</p>
                    <p>The vulnerability of each tumor clone to different drugs is assessed by considering only mutations with a relevant functional impact. Therapeutic options are prioritized by evaluating two fundamental dimensions:</p>
                    <ul>
                        <li>GScore (Biological Target): Quantifies the relevance of the mutated gene in cancer development and its viability as a pharmacological target (druggability).</li>
                        <li>DScore (Clinical Evidence): Estimates the suitability of the drug based on the current level of evidence, the type of interaction, and the phase of clinical development.</li>
                    </ul>
                    <p>We consider drugs with high clinical potential (DScore ≥ 0.7) that target key driver genes (GScore ≥ 0.6) to be Top Therapeutic Candidates (BTC), highlighted in bold.</p>
                    <p>The table below provides a simplified overview of the drugs targeting the affected genes, with the top 3 drugs selected per genetic alteration and ranked by Status and dScore. Drugs are sorted using the following criteria:</p>
                    <ul>
                        <li>Clone.</li>
                        <li>Mutation positon.</li>
                        <li>Status: APPROVED, CLINICAL_TRIALS, and EXPERIMENTAL.</li>
                        <li>Interaction type: DIRECT_TARGET, BIOMARKER, and PATHWAY_MEMBER.</li>
                    </ul>
                    <p>The source files are available at: {gene_alterations_path}.</p>
"""
            }
        },
        
        'key_concepts': [
    {
        'term': 'Clone',
        'definition': 'A population of tumor cells that all share a specific set of mutations inherited from a common ancestor cell. Ideally, each cluster of mutations identified by PyClone-VI would correspond to a biologically distinct clone.'
    },
    {
        'term': 'Cluster',
        'definition': 'Groups of cancer cells that share a common set of somatic mutations, which are inferred by clustering mutations with similar cellular prevalence estimates across one or more tumor samples.'
    },

    {
        'term': 'Driver Gene',
        'definition': 'Genes whose mutations are responsible for the development and progression of cancer.'
    },

        {
        'term': 'dScore',
        'definition': 'Pandrugs’ indicator which measures the suitability of the treatment for a particular patient. It ranges from -1 to 1, with the negative values corresponding to resistance and the positive values corresponding to sensitivity.'
    },

        {
        'term': 'gScore',
        'definition': 'Pandrugs’ indicator which measures the biological relevance of a gene in the tumoral process and its druggability. It ranges from 0 to 1, with higher values corresponding to more relevant and actionable targets.'
    },
        
        {
        'term': 'Interaction Type',
        'definition': 'Drugs approved for cancer treatment were manually classified into 7 different groups: chemotherapy, targeted therapy, hormone therapy, immunotherapy, photodynamic therapy, combination therapy, and other.'
    },

        {
        'term': 'Pandrugs',
        'definition': 'A bioinformatics platform to prioritize anticancer drug treatments according to individual multi-omics data.'
    },
    
    
    {
        'term': 'Parent cells',
        'definition': 'Group of cells that are set to be the origin of the tumor’s phylogeny while performing the clustering and the clonal tree. The estimation of the proportions of parent cells may vary depending on the samples provided.'
    },

    {
         'term': 'Study',   
         'definition': 'Group of biologically related samples that are analyzed together in Clonucopy. The experimental design and grouping criteria are defined according to specific research objectives, allowing multiple use cases adapted to the needs of the comparative analysis.'
    },
        
    {
        'term': 'Variant Allele Frequencies (VAF)',
        'definition': 'A biomarker that expresses the proportion of sequencing reads that support a specific variant allele relative to the total number of reads within a genomic locus. It ranges from 0 to 1.'
    }

            
]
    }

    
    try:
        html_content = template.render(**template_data)
        print("HTML successfully rendered")
    except Exception as e:
        print(f"Error rendering template: {e}")
        return False
    
    # CSS styles
    css_string = '''
    @page {
        size: A4;
        margin: 0.8cm;
        margin-top: 0.3cm;
    }

    img {
        max-width: 100% !important;
        height: auto !important;
    }

    /* Two-per-row layout for clonal proportions spheres. */
    .sphere-grid {
        width: 80% !important;
        margin: 4px auto !important;
        padding: 0 !important;
        font-size: 0 !important;
        line-height: 0 !important;
    
        page-break-inside: auto !important;
        break-inside: auto !important;
    }
    
    .sphere-item {
        display: inline-block !important;
        width: 48% !important;
        margin: 0 1% 8px 1% !important;
        padding: 0 !important;
        vertical-align: top !important;
        box-sizing: border-box !important;
    
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }
    
    .sphere-image {
        display: block !important;
        width: 100% !important;
        max-width: 100% !important;
        height: auto !important;
        margin: 0 auto !important;
    }



    .image-grid .image-container {
        display: inline-block !important;
        width: 48% !important;
        max-width: 48% !important;
        margin: 1% !important;
        vertical-align: top !important;
        box-sizing: border-box !important;
        font-size: 12px !important;
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }

    .image-grid .image-container img {
        max-width: 100% !important;
        height: auto !important;
    }

    .histogram-grid {
    display: block !important;
    width: 100% !important;
    text-align: center !important;
    font-size: 0 !important;
    margin-top: 12px !important;
    }
    
    .histogram-container {
        display: block !important;
        width: 100% !important;
        max-width: 100% !important;
        margin: 0 0 10px 0 !important;
        text-align: center !important;
        box-sizing: border-box !important;
        font-size: 12px !important;
        page-break-inside: avoid !important;
        break-inside: avoid !important;
        page-break-before: auto !important;
        page-break-after: auto !important;
    }
    
    .histogram-image {
        display: block !important;
        width: 65% !important;
        max-width: 65% !important;
        max-height: 95mm !important;
        height: auto !important;
        margin: 0 auto !important;
        object-fit: contain !important;
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }

    .heatmap-grid .image-container {
    display: inline-block !important;
    width: 98% !important;
    max-width: 98% !important;
    margin: 1% !important;
    vertical-align: top !important;
    box-sizing: border-box !important;
    font-size: 12px !important;
    page-break-inside: avoid !important;
    break-inside: avoid !important;
    }
    .heatmap-grid .image-container img {
        max-width: 100% !important;
        height: auto !important;
    }

    .data-table {
        width: 100% !important;
        font-size: 8px !important;
        table-layout: fixed !important;
    }
    .data-table th,
    .data-table td {
        padding: 2px 3px !important;
        word-wrap: break-word !important;
        overflow-wrap: break-word !important;
        font-size: 8px !important;
    }
    .data-table tr {
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }
    
    .panel {
        page-break-inside: auto !important;
        break-inside: auto !important;
        overflow: visible !important;
    }
    .stats-card {
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }

    .clonal-composition-figure {
        display: flex;
        justify-content: center;
        align-items: center;
        margin: 18px auto 8px auto;
        width: 100%;
        page-break-inside: avoid;
    }
    
    .clonal-composition-figure img {
        display: block !important;
        width: 75% !important;
        max-width: 650px !important;
        height: auto !important;
        margin: 0 auto !important;
        object-fit: contain !important;
    }

    h2, h3 {
        page-break-after: avoid !important;
        break-after: avoid !important;
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }
    '''
    
    try:
        base_url = os.path.dirname(os.path.abspath(template_path))
        html_doc = weasyprint.HTML(string=html_content, base_url=base_url)
        css_doc = weasyprint.CSS(string=css_string)
        
        pdf = html_doc.write_pdf(stylesheets=[css_doc])
        
        with open(output_path, 'wb') as f:
            f.write(pdf)
            
        print(f"Report successfully generated: {output_path}")
        return True
        
    except Exception as e:
        print(f"Error generating report: {str(e)}")
        return False




if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--study", action='store', required=True)
    input_parser.add_argument("--drug_filter", action='store', required=True, choices=["clinical", "discovery"], help="Filter drugs: 'clinical' or 'discovery'.")
    input_parser.add_argument("--output_pdf", action='store', required=True)
    input_parser.add_argument("--template", action='store', required=True)
    input_parser.add_argument("--logo", action='store', required=True)

    args = input_parser.parse_args()

    render_report_to_pdf(args.study, args.drug_filter, args.output_pdf, args.template, args.logo)
