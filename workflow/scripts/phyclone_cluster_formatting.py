import pandas as pd
import argparse


def format_clusters_pvi(output_pvi, clusters_chrom):
    """
    Process VCF file to extract mutation information and read counts.
    
    Args:
        input_vcf (str): Path to input TSV file, pyclone-vi output
        output_file (str, optional): Path to output TSV file. If None, won't save to file
        
    Returns:
        pandas.DataFrame: Processed clusters data with chromosomal information
    """
    # Read TSV file
    try:
        pvi_clusters = pd.read_csv(output_pvi, sep='\t')
    except Exception as e:
        raise ValueError(f"Error reading TSV file: {e}")

    # Extract mutation ids column with chromosomal info of the mutation
    mutation_ids = pvi_clusters["mutation_id"].copy()

    # Parse chrom info
    try:
        pvi_clusters['chrom'] = mutation_ids.str.split(':').str[0]
        
    except Exception as e:
        raise ValueError(f"Error parsing chromosomal info: {e}")

    # Save to output_file path

    pvi_clusters.to_csv(clusters_chrom, sep='\t', index=False)
    
    return pvi_clusters


    
if __name__ == '__main__':
    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input", action='store', required=True)
    input_parser.add_argument("--output", action='store', required=True)
    args = input_parser.parse_args()
    
    # Process the clusters file from pyclone-vi
    format_clusters_pvi(args.input, args.output)
