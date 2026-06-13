import pandas as pd
import glob
import re
import os
import ast
import argparse
import numpy as np
from clonucopya_tools import chr_to_num


def is_file_empty(file):
    try:
        df = pd.read_csv(file)
        return df.empty
    except Exception:
        return True


# Format samples' VAF per mutation
def vaf_string_for_mutation(mutation_id, vaf_dict):
    vaf_entries = []
    for sample_id, mut_vafs in vaf_dict.items():
        vaf_value = mut_vafs.get(mutation_id, 'NA')
        vaf_entries.append(f"{sample_id}: {vaf_value}")
    return "; ".join(vaf_entries)


# Build Gene Alterations dataframe with annotated information
def build_gene_alterations(pvi_input, pandrugs_dir, out_dir):

    """
    Build Daframe for sample panels of the study.

    Args:
        pvi_input (str): Path to the PyClone-VI input TSV
        pandrugs_dir (str): Path to the query_pandrugs study
        out_dir (str): Path to the output directory

    Return:
        DataFrame of the annotated gene alterations of the clonal inference

    """

    # Load Phyclone input
    pvi_df = pd.read_table(pvi_input)

    # Calculate VAF from ref_counts and alt_counts
    total_counts = pvi_df['ref_counts'] + pvi_df['alt_counts']
    pvi_df = pvi_df.copy()
    pvi_df['VAF'] = np.where(total_counts > 0, pvi_df['alt_counts'] / total_counts, 0)

    # Keep only needed columns from PyClone-VI input
    pvi_df = pvi_df[['sample_id', 'mutation_id', 'VAF']].copy()

    vaf_dict = {}
    for _, row in pvi_df.iterrows():
        sample = row['sample_id']
        mut_id = row['mutation_id']
        vaf = row['VAF']
        vaf_dict.setdefault(sample, {})[mut_id] = f"{vaf:.3f}"
    
    # Obtain relevant information about gene alterations fed to Pandrugs
    vscore_path = f"{pandrugs_dir}/clone_*/*_vscore.vcf"
    
    # Search all files that match the vscore pattern
    vscore_files = glob.glob(vscore_path)
    
    # Create list of dataframes to store all samples from the same study
    clone_dfs = []
    
    for file in vscore_files:
        file_name = os.path.basename(file)

        # Pass clone -1 (outlier)
        if "clone_-1" in file_name:
            continue

        match = re.search(r'clone_(-?\d+)', file)
        if match:
            clone = int(match.group(1))
        else:
            raise ValueError(f"No clone number found in filename: {file_name}")
        
        sample_df = pd.read_table(file)
        
        subset_df = sample_df[['ID', 'gene_hgnc', 'gene', 'consequence', 'impact']].copy()
        
        subset_df['clone_id'] = clone
    
        subset_df.drop_duplicates(subset=['ID'], inplace = True, ignore_index = True)
    
        clone_dfs.append(subset_df)
        
    study_df = pd.concat(clone_dfs, ignore_index=True)

    study_df['VAF'] = study_df['ID'].apply(
    lambda mut_id: vaf_string_for_mutation(mut_id, vaf_dict)
     )
    
    # Format Table: rename columns, sort rows and columns
    study_df.columns = ['Mutation ID', 'Gene Symbol', 'Ensembl ID', 'Consequence', 'Impact', 'Clone', 'VAF']
    
    
    study_df  = study_df.set_index('Mutation ID')
    temp_df = pd.DataFrame(index=study_df.index)
    temp_df['Clone'] = study_df['Clone']
    temp_df['chr_num'] = temp_df.index.to_series().apply(lambda x: chr_to_num(x.split(':')[0]))
    temp_df['pos'] = temp_df.index.to_series().apply(lambda x: int(x.split(':')[1]))
    temp_df_sorted = temp_df.sort_values(by=['Clone','chr_num', 'pos'])
    study_sorted = study_df.loc[temp_df_sorted.index]
    study_sorted = study_sorted.reset_index()
    study_sorted = study_sorted[['Clone', 'Mutation ID', 'Gene Symbol', 'Ensembl ID', 'VAF', 'Consequence', 'Impact']]

    # Save to TSV table
    study_sorted.to_csv(f"{out_dir}/gene_alterations.tsv", sep='\t', index=False)

    return study_sorted


def build_drug_prioritization(gene_alterations, pandrugs_dir, out_dir):

    """
    Build Daframe for sample panels of the study.

    Args:
        gene_alterations (str): Path to the dataframe of Phyclone results of the study
        pandrugs_dir (str): Path to the query_pandrugs study
        out_dir (str): Path to the output directory

    Return:
        DataFrame of Drug Prioritization from Pandrugs' query with extended clonal/gene/VAF information

    """

    # Search all files that match the gene-drug file pattern
    gd_path = f"{pandrugs_dir}/clone_*/*_gene-drug.csv"
    
    # Globbing the pattern
    gd_files = glob.glob(gd_path)

    
    # Filter out empty files
    non_empty_files = [file for file in gd_files if not is_file_empty(file)]
    
    # Load and concatenate non-empty files
    if non_empty_files:
        concatenated_df = pd.concat([pd.read_csv(file) for file in non_empty_files], ignore_index=True)
    else:
        # Save empy dataframe if all sample's drug-gene interactions files are empty
        concatenated_df = pd.DataFrame()

    
    if concatenated_df.empty:
        concatenated_df.to_csv(f"{out_dir}/drug_prioritization.tsv", sep='\t', index=False)
    else:
        # Subset Gene Alterations dataframe to get relevant columns for Drug Prioritizaton datataframe
        study_subset = gene_alterations[['Clone', 'Mutation ID', 'Gene Symbol', 'VAF']]
        gene_drugs = concatenated_df.copy()

        # Subset gen-drug concat file of all clones to resume information of the each query
        genalt_subset = gene_drugs[['gene', 'drug', 'status', 'interactionType', 'dScore', 'gScore' , 'cancer', 'source']].copy()

        # Format source column to clear the dataframe
        genalt_subset['source'] = genalt_subset['source'].apply(lambda x: '; '.join(ast.literal_eval(x)) if isinstance(x, str) and x.startswith("[") else None)
        genalt_subset['cancer'] = genalt_subset['cancer'].apply(lambda x: '; '.join(ast.literal_eval(x)) if isinstance(x, str) and x.startswith("[") else None)

        # Format gene column to obtain driverGene and geneSymbol information in new columns
        genalt_subset[['driverGene', 'geneSymbol']] = genalt_subset['gene'].apply(
        lambda x: pd.Series(ast.literal_eval(x)[0]) if isinstance(x, str) and x.startswith("[{") else pd.Series({'driverGene': None, 'geneSymbol': None})
    )[['driverGene', 'geneSymbol']]
    
        #  Sort and Rename columns of the subset datraframe from gene-drugs files
        genalt_subset_formatted = genalt_subset[['geneSymbol', 'drug', 'status', 'interactionType', 'driverGene', 'dScore', 'gScore' , 'cancer', 'source']]
        genalt_subset_formatted.columns = ['Gene Symbol', 'Drug', 'Status', 'Interaction Type', 'Driver Gene', 'dScore', 'gScore' , 'Cancer', 'Source']

        # Merge subsetted Gene alterations info with drug-gene interactions files
        drug_prioritization = (
        study_subset
        .merge(genalt_subset_formatted, on='Gene Symbol', how='inner')
    )
        # Remove duplicate queries
        drug_prioritization = drug_prioritization.drop_duplicates(subset=['Clone', 'Mutation ID', 'Gene Symbol', 'Drug', 'VAF'])        
        
        drug_prioritization.to_csv(f"{out_dir}/drug_prioritization.tsv", sep='\t', index=False, na_rep='-')

        # DRUG SUMMARY
        drug_hits_no_outliers = drug_prioritization[drug_prioritization["Clone"] != -1] 

        status_order = ["APPROVED", "CLINICAL_TRIALS", "EXPERIMENTAL"]
        
        drug_summary = (
        drug_hits_no_outliers
        .groupby('Drug')
        .agg(
            Status=('Status', 'first'),
            max_dScore=('dScore', lambda x: round(x.max(), 4)),
            Clones=('Clone', lambda x: ', '.join(sorted(x.unique().astype(str)))),
            Genes=('Gene Symbol', lambda x: ', '.join(sorted(x.unique().astype(str)))),
            n_clones=('Clone', 'nunique'),
        )
        .reset_index()
        )

        drug_summary["Status"] = pd.Categorical(
            drug_summary["Status"],
            categories=status_order,
            ordered=True
        )
        
        drug_summary = drug_summary.sort_values(
            ['n_clones', 'Status'],
            ascending=[False, True]
        )

        
        # Drop n_clones column
        drug_summary.drop(columns = ["n_clones"], axis=1, inplace=True)

        drug_summary.to_csv(f"{out_dir}/drug_summary.tsv", sep='\t', index=False)

    

    
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--pvi_input", action='store', required=True)
    input_parser.add_argument("--pandrugs_dir", action='store', required=True)
    input_parser.add_argument("--out_dir", action='store', required=True)

    args = input_parser.parse_args()


    gene_alterations = build_gene_alterations(args.pvi_input, args.pandrugs_dir, args.out_dir)
    build_drug_prioritization(gene_alterations, args.pandrugs_dir, args.out_dir)
