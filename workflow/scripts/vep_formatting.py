import os
import pandas as pd
import argparse
from typing import Dict, List
from pathlib import Path
from clonucopya_tools import chr_to_num



def process_pyclone_muts_clones(phy_out, study, out_dir):
    """
    Format output from Phyclone to VEP standard as input for ensembl-vep.
    
    Args:
        phy_out: Path to input TSV file containing Phyclone results
        study: Name of the sample
        out_dir: Path to output directory
    
    Returns:
        Dictionary mapping clone ids to dataframes containing mutation information
    """

    
    # Check and create output directory if it doesn't exist
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    
    # Load Phyclone output
    try:
        mut_info = pd.read_csv(phy_out, sep='\t')
    except Exception as e:
        raise ValueError(f"Failed to read Phyclone output file: {e}")


    # Initialize cluster dictionary
    clone_dataframes: Dict[int, List[Dict]] = {}
    
    # Process mutations
    for mut in mut_info['mutation_id'].unique():
        try:
            # Process chromosome prefix and split mutation components
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

            clone_id = mut_info.loc[mut_info['mutation_id'] == mut, 'clone_id'].iloc[0]

            if clone_id not in clone_dataframes:
                clone_dataframes[clone_id] = []
                
            clone_dataframes[clone_id].append({
                  'chr': chrom,
                  'start': pos,
                  'end': end,
                  'allele': f"{ref}/{alt}",
                  'strand': '+',
                  'mutation_id': mut
              })
        except Exception as e:
            raise ValueError(f"Error processing mutation {mut}: {e}")

    # Sort clone dfs and save them in the same output directory
    result: Dict[int, pd.DataFrame] = {}
    for clone_id, variants in clone_dataframes.items():
        df = pd.DataFrame(variants)
        df['start'] = df['start'].astype(int)
        df['chr_num'] = df['chr'].apply(chr_to_num)
        df.sort_values(by =['chr_num','start'], inplace=True)
        df.drop(columns=['chr_num'], inplace=True)
        output_file = out_path / f"{study}_cluster_{clone_id}.tsv"
        df.to_csv(output_file, sep='\t', index=False, header=False)
        result[clone_id] = df
    
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mut_data",  action='store', required=True)
    parser.add_argument("--study", action='store', required=True)
    parser.add_argument("--out_dir", action='store', required=True)
    
    args = parser.parse_args()
    process_pyclone_muts_clones(args.mut_data, args.study, args.out_dir)




