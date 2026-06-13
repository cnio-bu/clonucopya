import glob
import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from clonucopya_tools import chr_to_num



def even_distribution_tolerant(df_heatmap, n_samples=100):
    """
    Sampling using even distribution that deals with al sizes of heatmaps.
    If the Dataframe has less rows than n_sampes, it returns the whole dataframe.

    Args:
    df_heatmap: DataFrame of the VAF heatmap
    Return:
    Sampled Dataframe
    """
    total_rows = len(df_heatmap)
    
    # If DataFrame equal or smaller than the threshold
    if total_rows <= n_samples:
        return df_heatmap.copy()
    
    # Sampling for larger dataframes
    indices = np.linspace(0, total_rows-1, n_samples, dtype=int)
    
    # Drop duplicates preserving the order
    unique_indices = []
    seen = set()
    for idx in indices:
        if idx not in seen:
            unique_indices.append(idx)
            seen.add(idx)

    return df_heatmap.iloc[unique_indices]




def build_heatmap_df(tree_df, pvi_input, gene_alterations):

    """
    Build DataFrame for VAF heatmap, one per sample.

    Args:
        tree_df (str): Path to the PhyClone TSV output.
        pvi_input (str): Path to the PyClone-VI input TSV.

    Returns:
        dict: Dictionary with complete and sampled heatmap dataframes per sample.
    """

    # Load PhyClone output
    try:
        phy_df = pd.read_table(tree_df)
    except Exception as e:
        raise ValueError(f"Error reading PhyClone file: {e}")

    # Remove outlier mutations
    phy_df = phy_df[phy_df["clone_id"] != -1].copy()

    # Keep only needed columns from PhyClone
    phy_df = phy_df[['sample_id', 'mutation_id', 'clone_id']].drop_duplicates()

    # Load PyClone-VI input
    try:
        pvi_df = pd.read_table(pvi_input)
    except Exception as e:
        raise ValueError(f"Error reading PyClone-VI input file: {e}")

    # Calculate VAF from ref_counts and alt_counts
    total_counts = pvi_df['ref_counts'] + pvi_df['alt_counts']
    pvi_df = pvi_df.copy()
    pvi_df['VAF'] = np.where(total_counts > 0, pvi_df['alt_counts'] / total_counts, 0)

    # Keep only needed columns from PyClone-VI input
    pvi_df = pvi_df[['sample_id', 'mutation_id', 'VAF']].copy()

    # Merge PyClone-VI input with PhyClone output to recover clone_id
    merged_df = pvi_df.merge(
        phy_df,
        on=['sample_id', 'mutation_id'],
        how='inner'
    )

    # Group by sample_id directly from merged dataframe
    pvi_dict = {
        sample_id: sub_df.copy()
        for sample_id, sub_df in merged_df.groupby('sample_id')
    }

    # Track all mutations included after merge
    all_muts_list = merged_df['mutation_id'].drop_duplicates().tolist()


        # Load gene alterations table
    try:
        gene_alt = pd.read_table(gene_alterations)
    except Exception as e:
        raise ValueError(f"Error reading {gene_alterations}: {e}")
        
    # Keep only HIGH and MODERATE impact mutations
    gene_alt = gene_alt[gene_alt['Impact'].isin(['HIGH', 'MODERATE'])].copy()

    heatmap_dict = {}
    heatmap_dict_sampled = {}

    for sample_id, sample_df in pvi_dict.items():

        # Keep only required columns
        df = sample_df[['mutation_id', 'clone_id', 'VAF']].copy()

        # If repeated mutation_id/clone_id pairs exist, aggregate them
        df = df.groupby(['mutation_id', 'clone_id'], as_index=False)['VAF'].mean()

        # Build heatmap matrix
        pivoted = df.pivot(index='mutation_id', columns='clone_id', values='VAF')
        pivoted = pivoted.fillna(0)

        # Add absent mutations as zero-VAF rows
        df_reindexed = pivoted.reindex(all_muts_list, fill_value=0)

        # Sort mutations by chromosome and position
        temp_df = pd.DataFrame(index=df_reindexed.index)
        temp_df['chr_num'] = temp_df.index.to_series().apply(lambda x: chr_to_num(x.split(':')[0]))
        temp_df['pos'] = temp_df.index.to_series().apply(lambda x: int(x.split(':')[1]))
        temp_df_sorted = temp_df.sort_values(by=['chr_num', 'pos'])
        df_reindexed = df_reindexed.loc[temp_df_sorted.index]

        # Remove mutations with VAF 0 in all clones
        df_no_empty_muts = df_reindexed.loc[(df_reindexed != 0).any(axis=1)]

        # Merge with gene alteration annotations
        temp_mut_idx_df = df_no_empty_muts.reset_index()
        mut_gene_idx_df = temp_mut_idx_df.merge(
            gene_alt[['Mutation ID', 'Gene Symbol']],
            left_on='mutation_id',
            right_on='Mutation ID',
            how='inner'
        )

        # Create new index with mutation ID + gene symbol
        mut_gene_idx_df['mut_gene_idx'] = mut_gene_idx_df.apply(
            lambda row: f"{row['mutation_id']} - {row['Gene Symbol']}"
            if pd.notna(row['Gene Symbol']) else row['mutation_id'],
            axis=1
        )

        # Set new index and drop auxiliary columns
        mut_gene_idx_df = mut_gene_idx_df.set_index('mut_gene_idx').drop(
            columns=['mutation_id', 'Mutation ID', 'Gene Symbol']
        )

        # Reorder clones
        ordered_cols = sorted(mut_gene_idx_df.columns)

        # Store complete and sampled heatmaps
        heatmap_dict[sample_id] = mut_gene_idx_df[ordered_cols]
        heatmap_dict_sampled[sample_id] = even_distribution_tolerant(mut_gene_idx_df[ordered_cols])

    heatmap_data = {
        'complete': heatmap_dict,
        'sampled': heatmap_dict_sampled
    }

    return heatmap_data




def plot_heatmaps(heatmap_dicts, out_dir):
    """
    Plot VAF Heatmap, one per sample.

    Args:
        heatmap_dicts (str): dictionary of dictionaries of complete and sampled heatmap of samples' dataframes to plot the VAF Heatmap
        out_file (str): Path to output file of the study.
        
    """
    for data_type, heatmap_dict in heatmap_dicts.items():
        # Create directory if it does not exist
        type_dir = os.path.join(out_dir, data_type)
        os.makedirs(type_dir, exist_ok=True)
        
        print(f"Plotting {data_type} heatmaps...")
        
        for sample_id, df in heatmap_dict.items():
            n_rows, n_cols = df.shape
            fig, ax = plt.subplots(figsize=(max(12, n_cols*4.5), max(5, n_rows*0.2)))
            sns.set(font_scale=0.9)
            
            # Set color palette
            cmap = sns.light_palette("#009c8c", as_cmap=True)
            
            # Plot heatmap
            heatmap = sns.heatmap(
                df,
                annot=False,
                fmt=".3f",
                cmap=cmap,
                linewidths=0.5,
                vmax=0.6
            )
            
            heatmap.xaxis.tick_top()
            heatmap.tick_params(axis='x', which='both', pad=10, top=True, bottom=False, length=0)
            heatmap.tick_params(axis='y', which='both', left=False, right=False, length=0)
            heatmap.xaxis.set_label_position('top') 
            
            # Adjust color bar params
            cbar = heatmap.collections[0].colorbar
            cbar.ax.tick_params(axis='x', labelrotation=0) 
            cbar.ax.set_xlabel("VAF", fontsize=16, labelpad=10)
            cbar.ax.xaxis.set_label_position('top')
        
            # Set axis labels and title
            title = f"{sample_id} ({data_type.capitalize()})"
            plt.title(title, fontsize=20, pad=25)
            plt.xlabel("Clones", fontsize=16, labelpad=15)
            plt.ylabel("Mutations", fontsize=16, labelpad=15)
            
            plt.savefig(f"{type_dir}/{sample_id}_{data_type}_vaf_heatmap.png", bbox_inches='tight')
            plt.close()



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--tree_df", action='store', required=True)
    input_parser.add_argument("--pvi_input", action='store', required=True)
    input_parser.add_argument("--gene_alterations", action='store', required=True)
    input_parser.add_argument("--out_dir", action='store', required=True)


    args = input_parser.parse_args()


    heatmap_data = build_heatmap_df(args.tree_df, args.pvi_input, args.gene_alterations)
    plot_heatmaps(heatmap_data, args.out_dir)
