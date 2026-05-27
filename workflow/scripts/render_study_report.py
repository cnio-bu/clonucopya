import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, Template
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




def prepare_report_data(study_path):
    """
    Prepare data for the template injection with absolute paths
    """
    study_name = os.path.basename(study_path)
    components_path = Path(study_path) / "report" / "components"
    
    # Load report panel
    samples_df = pd.read_csv(components_path / "report_panel.tsv", sep='\t')

    # Drug Summary
    drug_sum_path = components_path / "drug_summary.tsv"
    drug_sum_top = None
    if drug_sum_path.exists():
        try:
            drug_sum = pd.read_csv(drug_sum_path, sep='\t')
            if not drug_sum.empty:
                drug_sum_top = drug_sum.head(25)
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
                
                # Grouping key
                group_keys = ["Clone", "Mutation ID", "Gene Symbol", "VAF"]
                
                drugs_df_sorted = drug_hits.sort_values(
                    by=group_keys + ["Status", "dScore", "Interaction Type"],
                    ascending=[True, True, True, True, True, False, True]
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
                is_dup = ~drugs_df_compact[collapse_cols].ne(drugs_df_compact[collapse_cols].shift()).any(axis=1)
                drugs_df_compact.loc[is_dup, collapse_cols] = ""
            
            
        except:
            drug_hits = None
    
    # Absolute paths to images
    clonal_tree_image = os.path.abspath(str(components_path / "clonal_tree.png"))
    
    clonal_proportions_images = {}
    clone_alterations_images = {}

    # Candidate locations / filename patterns for VAF heatmaps.
    heatmap_candidates = [
        ("vaf_heatmaps", "sampled", "{sample_id}_sampled_vaf_heatmap.png"),
        ("vaf_heatmaps", "{sample_id}_sampled_vaf_heatmap.png"),
        ("vaf_heatmaps", "{sample_id}_vaf_heatmap.png"),
        ("vaf_heatmaps", "sampled", "{sample_id}_vaf_heatmap.png"),
    ]

    for _, sample in samples_df.iterrows():
        sample_id = sample['sample_id']

        # Sphere of clones
        sphere_path = components_path / "spheres_of_clones" / f"{sample_id}_sphere_of_clones.png"
        if sphere_path.exists():
            clonal_proportions_images[sample_id] = os.path.abspath(str(sphere_path))
        else:
            print(f"[WARN] sphere of clones not found for {sample_id}: {sphere_path}")

        # VAF heatmap — try several known locations/patterns
        found_heatmap = None
        tried = []
        for parts in heatmap_candidates:
            resolved_parts = [p.format(sample_id=sample_id) for p in parts]
            candidate = components_path.joinpath(*resolved_parts)
            tried.append(str(candidate))
            if candidate.exists():
                found_heatmap = candidate
                break

        if found_heatmap is not None:
            clone_alterations_images[sample_id] = os.path.abspath(str(found_heatmap))
        else:
            print(f"[WARN] VAF heatmap not found for {sample_id}. Tried:")
            for t in tried:
                print(f"        - {t}")
    
    return {
        'samples_df': samples_df,
        'drug_summary': drug_sum_top,
        'gene_alterations_df': gene_alterations_top,
        'drug_prioritization_df': drugs_df_compact,
        'clonal_tree_image': clonal_tree_image,
        'clonal_proportions_images': clonal_proportions_images,
        'clone_alterations_images': clone_alterations_images
    }


def render_report_to_pdf(study_path, output_path, template_path="template.html", logo_path=None):
    """
    Render template to PDF directly
    """

    # Set path to source files

    study_name = os.path.basename(study_path)
    components_path = Path(study_path) / "report" / "components"
    
    panels_path = components_path / "report_panel.tsv"
    drug_summary_path = components_path / "drug_summary.tsv"
    clonal_tree_path = components_path / "clonal_tree.png"
    spheres_path = components_path / "spheres_of_clones" / "{sample_id}_sphere_of_clones.png"
    heatmaps_path = components_path / "vaf_heatmaps" / "sampled" / "{sample_id}_sampled_vaf_heatmap.png"
    gene_alterations_path = components_path / "gene_alterations.tsv"
    drug_prioritization_path = components_path / "drug_prioritization.tsv"
    
    
    
    # Format data
    data = prepare_report_data(study_path)
    
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
        'report_title': 'Clonal Evolution Report',
        'study_id': study_name,
        'study_description': f"""
        <p>Overview of study {study_name}'s key statistics by sample, including the number of sample's name, sex, number of mutations such as Single Nucleotide Variation (SNVs) or small Indels, Copy Number Variations (CNVs), and the number of intersections between mutations and CNVs.</p>
        <p>The intersection statistics indicate the number of matches that occurred between SNVs and/or small indels contained in CNVs during preprocessing prior to clonal inference using Pyclone-VI and Phyclone. Intersections that have passed the Pyclone-VI filter can be found in the Clonal Alterations and Gene Alterations sections.</p>
        <p>The Pyclone-VI filter is based on the following aspects:</p>
        <ul>
            <li>Mutations that have a major copy number equals to 0.</li>
            <li>Mutations with missing coverage data (ref and/or alt counts).</li>
        </ul>

        <p>The source files are available at: {panels_path}.</p>
""",
        'logo_path': logo_data_uri,
        
    
        'samples_df': data['samples_df'],
        'drug_summary': data['drug_summary'],
        'gene_alterations_df': data['gene_alterations_df'],
        'drug_prioritization_df': data['drug_prioritization_df'],
        'clonal_proportions_images': data['clonal_proportions_images'],
        'clone_alterations_images': data['clone_alterations_images'],

            'drug_summary_description': f"""
            <p>The Drug Summary provides a high-level overview of the therapeutic candidates identified across the entire study. For each drug, the table reports the number of genetic alterations supporting its prioritization, its regulatory approval status (APPROVED, CLINICAL_TRIALS, or EXPERIMENTAL), and the type of interaction with the affected genes
            (DIRECT_TARGET, BIOMARKER, or PATHWAY_MEMBER).</p>
            <p>This table summarizes up to 25 top-ranked drug candidates derived from the mutational landscape of all samples included in the study. A detailed per-clone breakdown is available in the Drug Prioritization section.</p>
            
            <p>The source files are available at: {drug_summary_path}.</p>
            """,
        
        'comparison_section': {
            'title': 'Clonal Evolution Analysis',
            'description': 'Analysis of clonal evolution patterns and proportions across samples.',
            'subclonal_tree': {
                'title': 'Clonal Tree',
                'description': f"""
                <p>The tree is inferred using Phyclone from Pyclone-VI results. Pyclone-VI results provide an initial clustering of the clones but have no order, so they do not follow an established phylogeny.</p>
                <p>The source files are available at: {clonal_tree_path}.</p>
                """,
                'image': data['clonal_tree_image']
            },
            'clonal_proportions': {
                'title': 'Clonal Proportions',
                'description': f"""
                <p>Based on the phylogeny results obtained by Phyclone, the spheres of clones show the distribution of the clonal population in each sample.</p>
                <p>The source files are available at: {spheres_path}.</p>
                """
            },
            'clone_alterations': {
                'title': 'Clonal Alterations',
                'description': f"""
                <p>The Variant Allele Frequencies (VAF) heatmaps of each sample represent the evolutionary dynamics within the tumor(s). The distribution of VAF intensity patterns for each mutation in the different clones evidences how certain mutations are shared in the same clone from different samples or produce divergence and give rise to different evolutionary branches, representing the temporal sequence of mutational events. It allows to understand tumor heterogeneity among cell populations, with the aim of highlighting which clones may influence cancer progression or therapy resistance.</p>
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
                    <p>The small variant analysis performed by PanDrugs2 displays only mutations deemed clinically relevant.</p>
                    <p>The filtering criteria are as follows:</p>
                    <ul>
                        <li>GMAF/gnomAD population frequency less than 0.01.</li>
                        <li>Predicted moderate or high functional impact, including variant types such as missense, nonsense, frameshift, and splice site mutations.</li>
                        <li>Affection of relevant isoforms. Priority is given to canonical or unknown isoforms.</li>
                    </ul>
                    <p>The table below provides a simplified overview of the drugs targeting the affected genes, with the top 3 drugs selected per genetic alteration and ranked by Status and dScore.</p>
                    <p> Drugs are sorted using the following criteria:</p>
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
    }

    img {
        max-width: 100% !important;
        height: auto !important;
    }

    /* Two-per-row layout for clonal proportions spheres. */
    .sphere-container {
        display: block !important;
        width: 100% !important;
        text-align: center !important;
        margin: 20px 0 !important;
        font-size: 0 !important; /* removes whitespace gaps between inline-blocks */
    }

    .sphere-image {
        display: inline-block !important;
        width: 48% !important;
        max-width: 48% !important;
        height: auto !important;
        margin: 1% !important;
        vertical-align: top !important;
        box-sizing: border-box !important;
        page-break-inside: avoid !important;
        break-inside: avoid !important;
    }
    
    .image-grid {
        display: block !important;
        width: 100% !important;
        text-align: center !important;
        font-size: 0 !important;
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
    
    /* Allow panels to split across pages */
    .panel {
        page-break-inside: auto !important;
        break-inside: auto !important;
        overflow: visible !important;
    }
    .stats-card {
        page-break-inside: avoid !important;
        break-inside: avoid !important;
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
    input_parser.add_argument("--output_pdf", action='store', required=True)
    input_parser.add_argument("--template", action='store', required=True)
    input_parser.add_argument("--logo", action='store', required=True)

    args = input_parser.parse_args()

    renderization = render_report_to_pdf(args.study, args.output_pdf, args.template, args.logo)
