import pandas as pd
import argparse


def format_clusters_pvi(output_pvi, clusters_chrom):
    """
    Process Pyclone-VI file to adapt it to Phyclone input.
    
    Args:
        output_pvi (str): Path to input TSV file, pyclone-vi output
        clusters_chrom (str): Path to output TSV file.
        
    Returns:
        pandas.DataFrame: Processed clusters data with chromosomal information
    """
    # Load Pyclone-vi's clustered samples
    try:
        pvi_clusters = pd.read_csv(output_pvi, sep='\t')
    except Exception as e:
        raise ValueError(f"Error reading TSV file: {e}")

    # Select mutation ids with chromosomal info of the mutations
    mutation_ids = pvi_clusters["mutation_id"].copy()

    # Parse chrom info
    try:
        pvi_clusters['chrom'] = mutation_ids.str.split(':').str[0]
        
    except Exception as e:
        raise ValueError(f"Error parsing chromosomal info: {e}")


    pvi_clusters.to_csv(clusters_chrom, sep='\t', index=False)
    
    return pvi_clusters


    
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input", action='store', required=True)
    input_parser.add_argument("--output", action='store', required=True)
    args = input_parser.parse_args()
    
    format_clusters_pvi(args.input, args.output)
