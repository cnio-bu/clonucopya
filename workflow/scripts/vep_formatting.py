import os
import pandas as pd
import argparse
from typing import Dict, List
from pathlib import Path
from clonucopya_tools import chr_to_num



def filt_mock_mutations(pvi_in, pvi_out):

    """
    Filter out mock mutations before formatting vep input.
    
    Args:
        pvi_int: Path to input TSV file containing PyClone-VI input
        pvi_out: Path to input TSV file containing PyClone-VI results
    
    Returns:
        DataFrame contained filtered Pyclone-VI output witout mock mutations
    """

    # Filter mock mutations out of Pyclone-VI input
    pvi_mock = pvi_in[((pvi_in.ref_counts == 0) & (pvi_in.alt_counts == 0))]


    # Filter mock mutations out of Pyclone-VI output using filtered pvi input
    pvi_filt = (
        pvi_out
        .merge(
            pvi_mock[['mutation_id', 'sample_id']],
            on=['mutation_id', 'sample_id'],
            how='left',
            indicator=True
        )
        .query('_merge == "left_only"')
        .drop(columns=['_merge'])
    )
    return pvi_filt




def process_pyclone_muts_clones(pvi_in, pvi_out, study, out_dir):
    """
    Format output from PyClone-VI to VEP standard as input for ensembl-vep.
    
    Args:
        pvi_int: Path to input TSV file containing PyClone-VI input
        pvi_out: Path to input TSV file containing PyClone-VI results
        study: Name of the sample
        out_dir: Path to output directory
    
    Returns:
        Dictionary mapping cluster ids to dataframes containing mutation information
    """

    
    # Check and create output directory if it doesn't exist
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Load PyClone-VI input
    try:
        pvi_prep = pd.read_csv(pvi_in, sep='\t')
    except Exception as e:
        raise ValueError(f"Failed to read PyClone-VI input file: {e}")
    
    # Load PyClone-VI output
    try:
        pvi_results = pd.read_csv(pvi_out, sep='\t')
    except Exception as e:
        raise ValueError(f"Failed to read PyClone-VI output file: {e}")

    # Filter mock mutations out of Pyclone-VI output
    pvi_filt = filt_mock_mutations(pvi_prep, pvi_results)

    # Initialize cluster dictionary
    clone_dataframes: Dict[int, List[Dict]] = {}
    
    # Process mutations
    for mut in pvi_filt['mutation_id'].unique():
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

            clone_id = pvi_filt.loc[pvi_filt['mutation_id'] == mut, 'cluster_id'].iloc[0]

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
    parser.add_argument("--pvi_prep",  action='store', required=True)
    parser.add_argument("--pvi_data",  action='store', required=True)
    parser.add_argument("--study", action='store', required=True)
    parser.add_argument("--out_dir", action='store', required=True)
    
    args = parser.parse_args()
    process_pyclone_muts_clones(args.pvi_prep, args.pvi_data, args.study, args.out_dir)




