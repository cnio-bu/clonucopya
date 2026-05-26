import pandas as pd
import itertools
import sys
import argparse




def format_intersect(intersect_list, samplesheet, pvi_prep, phyclone):
    """
        Build DataFrames for Pyclone-VI and Phyclone, one per sample.

    Args:
        intersect_list (str): spaced list of intersect dataframes of each sample (TSV)
        samplesheet (str): Path to the study samplesheet
        pvi_prep (str): Path to formatted dataframe for Pyclone-VI (TSV).
        phyclone (str):  Path to formatted dataframe for Phyclone (TSV).

    Return:
        Dictionary of samples' dataframes to plot the VAF Heatmap
    """
    # List of pvi prep samples

    pvi_preps = [pd.read_csv(file, sep='\t') for file in intersect_list]

    # Concatenate
    combined_pvi = pd.concat(pvi_preps, ignore_index=True)

    # Drop artifactual duplicates
    combined_pvi_dedup = combined_pvi.drop_duplicates()

    # Add mock mutations

    ## Scan available mutations
    mutations = set(combined_pvi_dedup['mutation_id'].unique())
    samples = set(combined_pvi_dedup['sample_id'].unique())

    ## Get all hypothetical combinations 
    mutations_complete = set(itertools.product(mutations, samples))

    ## Get existing mutations
    true_mutations = set(zip(combined_pvi_dedup['mutation_id'], combined_pvi_dedup['sample_id']))

    ## Get missing mutations
    missing_mutations = mutations_complete - true_mutations

    ## Sort completed df by mutation_id and sample_id
    combined_pvi_dedup.sort_values(by=['mutation_id', 'sample_id'], inplace=True)

    # Add tumour_content aka purity
    samplesheet = pd.read_csv(samplesheet)
    tumour_content_dict = dict(zip(samplesheet['sample_id'], samplesheet['tumour_content']))
    combined_pvi_dedup['tumour_content'] = combined_pvi_dedup['sample_id'].map(tumour_content_dict)

    # Save pyclone-vi formatted intersect df
    combined_pvi_dedup.to_csv(pvi_prep, sep='\t', index=False)

    # Subset and save intersect df for phyclone input format
    phyclone_df = combined_pvi_dedup[['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'major_cn', 'minor_cn', 'normal_cn', 'tumour_content']]
    phyclone_df.to_csv(phyclone, sep='\t', index=False)

    return combined_pvi_dedup, phyclone_df


if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--intersect_list", action='store', nargs='+', required=True)
    input_parser.add_argument("--samplesheet", action='store', required=True)
    input_parser.add_argument("--pvi_prep", action='store', required=True)
    input_parser.add_argument("--phyclone", action='store', required=True)
    args = input_parser.parse_args() 
    
    format_intersect(args.intersect_list, args.samplesheet, args.pvi_prep, args.phyclone)
