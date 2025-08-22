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




def build_heatmap_df(tree_df, pvi_out, mut_dir, gene_alterations):

    """
    Build Daframe for VAF Heatmap, one per sample.

    Args:
        tree_df (str): Path to the phyclone TSV output
        pvi_out (str): Path to the pyclone-vi TSV output
        mut_dir (str): Path to the files from mutation_prep study.
        gene_alterations (str): Path to the gene alterations file of the study

    Return:
        Dictionary of dictioraries with complete and sampled data of samples' dataframes to plot the VAF Heatmap
        
    """
    # Load Phyclone output
    try:
        phy_df = pd.read_table(tree_df)
    except Exception as e:
        raise ValueError(f"Error reading Phyclone file: {e}")
    
    # Obtain cluter-clone equivalences
    phy_clones = phy_df[['clone_id', 'cluster_id']].drop_duplicates()
    cluster_to_clone = dict(zip(phy_clones['cluster_id'], phy_clones['clone_id']))
    
    # PYCLONE-VI INFO
    # Load Pyclone-VI output
    try:
        pvi_df = pd.read_table(pvi_out)
    except Exception as e:
        raise ValueError(f"Error reading Pyclone-VI file: {e}")
    
    # Subset Pyclone-VI dataframe by samples
    pvi_dict = {sample_id: sub_df for sample_id, sub_df in pvi_df.groupby('sample_id')}
    
    # Track all possible mutations that have been taken into account in Pyclone-VI inference
    all_muts = pvi_df[['mutation_id']].drop_duplicates()
    all_muts_list = all_muts['mutation_id'].tolist()
    
    
    # MUTATIONS INFO
    # Load samples' mutations
    muts_path = f"{mut_dir}/*.tsv"
    mut_files = glob.glob(muts_path)
    
    # Store all mutations dataframe in a dict
    mut_dict = {}

    for file in mut_files:
        file_name = os.path.basename(file)
        sampleid = file_name.split('_prep.mut.tsv')[0]
        
        try:
            sample_df = pd.read_table(file)
        except Exception as e:
            raise ValueError(f"Error reading {file_name}: {e}")
        
        mut_dict[sampleid] = sample_df
    
    # GENE ALTERATIONS INFO
    # Load gene alterations table from report components
    try:
        gene_alt = pd.read_table(gene_alterations)
    except Exception as e:
        raise ValueError(f"Error reading {gene_alterations}: {e}")
        
    
    # INNER JOINT OF MUT AND PVI
    # Obtain a dict of samples' dataframe with ['mutation_id','cluster_id','VAF']
    heatmap_dict = {}
    heatmap_dict_sampled = {}
    
    for sample_id in mut_dict:
        merged = mut_dict[sample_id][['mutation_id', 'VAF']].merge(
            pvi_dict[sample_id][['mutation_id', 'cluster_id']],
            on='mutation_id',
            how='inner')
        heatmap_dict[sample_id] = merged[['mutation_id', 'cluster_id', 'VAF']]
    
    
    # FORMAT DATAFRAMES AND INPUT VALUES
    
    for sample_id in heatmap_dict:
    
        # Reshape dataframe and input NA values with 0
        df = heatmap_dict[sample_id]
        pivoted = df.pivot(index='mutation_id', columns='cluster_id', values='VAF')
        pivoted.columns = [cluster for cluster in pivoted.columns]
        pivoted = pivoted.fillna(0)
        
        # Add mutations that are absent in each sample and reindex
        df_reindexed = pivoted.reindex(all_muts_list, fill_value=0)
        
    	# Sort mutations by clone, chr, and pos
        df_reindexed = df_reindexed.sort_index()
        temp_df = pd.DataFrame(index=df_reindexed.index)
        temp_df['chr_num'] = temp_df.index.to_series().apply(lambda x: chr_to_num(x.split(':')[0]))
        temp_df['pos'] = temp_df.index.to_series().apply(lambda x: int(x.split(':')[1]))
    
        # Sort temporary Dataframe mutations by chromosome and position
        temp_df_sorted = temp_df.sort_values(by=['chr_num', 'pos'])
        
        # Sort the main Dataframe using sorted index from temporary DataFrame
        df_reindexed = df_reindexed.loc[temp_df_sorted.index]

        ## Remove with 0 VAF in all clones

        df_no_empty_muts = df_reindexed.loc[(df_reindexed != 0).any(axis=1)]

        # Reset index to convert mutation_id into a column
        temp_mut_idx_df = df_no_empty_muts.reset_index()

        # Merge sorted index DataFrame with gene_alteration Dataframe to obtain Gene Symbol information
        mut_gene_idx_df = temp_mut_idx_df.merge(gene_alt[['Mutation ID', 'Gene Symbol', 'Impact']], 
            left_on='mutation_id', 
            right_on='Mutation ID', 
            how='left')

        # Only keep mutations with High or Moderate Impact annotated by vep
        mut_gene_idx_df = mut_gene_idx_df[(mut_gene_idx_df['Impact'] == 'MODERATE') | (mut_gene_idx_df['Impact'] == 'HIGH')]
    
        # Create a new mutation_id index combining mutation_id and its gene symbol
        mut_gene_idx_df['mut_gene_idx'] = mut_gene_idx_df.apply(
            lambda row: f"{row['mutation_id']} - {row['Gene Symbol']}" 
            if pd.notna(row['Gene Symbol']) 
            else row['mutation_id'], 
            axis=1)
    
        # Set mut_gene_idx as the new index and drop auxiliary columns
        mut_gene_idx_df = mut_gene_idx_df.set_index('mut_gene_idx').drop(['mutation_id', 'Mutation ID', 'Gene Symbol', 'Impact'], axis=1)

        
        # Transform cluster to clones
        df_renamed = mut_gene_idx_df.rename(columns=cluster_to_clone)
    
        # Reorder clones. Sometimes, clusters and clones dont have the same order
        ordered_cols = sorted(df_renamed.columns)
    
        # Update original dataframes
        heatmap_dict[sample_id] = df_renamed[ordered_cols]
        heatmap_dict_sampled[sample_id] = even_distribution_tolerant(df_renamed[ordered_cols])

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
    input_parser.add_argument("--pvi_out", action='store', required=True)
    input_parser.add_argument("--mut_dir", action='store', required=True)
    input_parser.add_argument("--gene_alterations", action='store', required=True)
    input_parser.add_argument("--out_dir", action='store', required=True)


    args = input_parser.parse_args()


    heatmap_data = build_heatmap_df(args.tree_df, args.pvi_out, args.mut_dir, args.gene_alterations)
    plot_heatmaps(heatmap_data, args.out_dir)
