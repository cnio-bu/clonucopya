import os
import pandas as pd
import argparse
from typing import Dict, List
from pathlib import Path

def process_pyclone_snps_clones(pvi_file, sample_id, out_dir):
    """
    Format output from PyClone-VI to VEP standard as input for ensembl-vep.
    
    Args:
        pvi_file: Path to input TSV file containing PyClone-VI results
        sample_id: Name of the sample
        out_dir: Path to output directory
    
    Returns:
        Dictionary mapping cluster IDs to DataFrames containing SNP information
    """
    # Create output directory if it doesn't exist
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    # Read PyClone-VI data
    try:
        pvi_data = pd.read_csv(pvi_file, sep='\t')
    except Exception as e:
        raise ValueError(f"Failed to read PyClone-VI file: {e}")

    # Initialize storage for cluster data
    clone_dataframes: Dict[int, List[Dict]] = {}
    
    # Process mutations
    for mut in pvi_data['mutation_id'].unique():
        try:
            # Clean chromosome prefix and split mutation components
            chrom, pos, ref, alt = mut.split(':')
            if ref == '-':
                # Insertion
                end = int(pos) - 1
            elif alt == '-':
                # Deletion
                end = int(pos) + len(ref) - 1
            else:
                # Substitution
                end = int(pos) + len(ref) - 1

            clone_id = pvi_data.loc[pvi_data['mutation_id'] == mut, 'cluster_id'].iloc[0]

            if clone_id not in clone_dataframes:
                clone_dataframes[clone_id] = []
                
            clone_dataframes[clone_id].append({
                  'chr': chrom,
                  'start': pos,
                  'end': end,
                  'allele': f"{ref}/{alt}"
              })
        except Exception as e:
            raise ValueError(f"Error processing mutation {mut}: {e}")
    
    # Convert to DataFrames and save
    result: Dict[int, pd.DataFrame] = {}
    for clone_id, variants in clone_dataframes.items():
        df = pd.DataFrame(variants)
        output_file = out_path / f"{sample_id}_cluster_{clone_id}.tsv"
        df.to_csv(output_file, sep='\t', index=False, header=False)
        result[clone_id] = df
    
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvi_data",  action='store', required=True)
    parser.add_argument("--sample_id", action='store', required=True)
    parser.add_argument("--out_dir", action='store', required=True)
    
    args = parser.parse_args()
    process_pyclone_snps_clones(args.pvi_data, args.sample_id, args.out_dir)




