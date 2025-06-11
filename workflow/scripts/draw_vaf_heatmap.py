import glob
import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


def build_heatmap_df(tree_df,pvi_out, mut_project):

    """
    Build Daframe for VAF Heatmap, one per sample.

    Args:
        tree_df (str): Path to the phyclone TSV output
        pvi_out (str): Path to the pyclone-vi TSV output
        mut_project (str, optional): Path to the mutation_prep project.

    Return:
        Dictionary of samples' dataframes to plot the VAF Heatmap
        
    """
    # Load Phyclone output
    phy_df = pd.read_table(tree_df)
    
    # Obtain cluter-clone equivalences
    phy_clones = phy_df[['clone_id', 'cluster_id']].drop_duplicates()
    cluster_to_clone = dict(zip(phy_clones['cluster_id'], phy_clones['clone_id']))
    
    # PYCLONE-VI INFO
    # Load Pyclone-VI output
    pvi_df = pd.read_table(pvi_out)
    
    # Subset Pyclone-VI dataframe by samples
    pvi_dict = {sample_id: sub_df for sample_id, sub_df in pvi_df.groupby('sample_id')}
    
    # Track all possible mutations that have been taken into account in Pyclone-VI inference
    all_muts = pvi_df[['mutation_id']].drop_duplicates()
    all_muts_list = all_muts['mutation_id'].tolist()
    
    
    # MUTATIONS INFO
    # Load samples' mutations
    muts_path = f"{mut_project}/*.tsv"
    mut_files = glob.glob(muts_path)
    
    # Store all mutations dataframe in a dict
    mut_dict = {}
    
    for file in mut_files:
        file_name = os.path.basename(file)
        sampleid = file_name.split('_prep.mut.tsv')[0]
        sample_df = pd.read_table(file)
        
        mut_dict[sampleid] = sample_df
    
    
    # INNER JOINT OF MUT AND PVI
    # Obtain a dict of samples' dataframe with ['mutation_id','cluster_id','VAF']
    heatmap_dict = {}
    
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
        # Sort mutations alpha-numerically
        df_reindexed = df_reindexed.sort_index(ascending=False)
    
        # Transform cluster to clones
        df_renamed = df_reindexed.rename(columns=cluster_to_clone)
    
        # Reorder clones. Sometimes, clusters and clones dont have the same order
        ordered_cols = sorted(df_renamed.columns)
    
        # Update original dataframe
        heatmap_dict[sample_id] = df_renamed[ordered_cols]

    return heatmap_dict




def plot_heatmaps(heatmap_dict, out_dir):
    """
    Plot VAF Heatmap, one per sample.

    Args:
        heatmap_dict (str): dictionary of samples' dataframes to plot the VAF Heatmap
        out_dir (str, optional): Path to output directory of the project.
        
    """
    for sample_id, df in heatmap_dict.items():
        n_rows, n_cols = df.shape
        fig, ax = plt.subplots(figsize=(max(12, n_cols*4.5), max(5, n_rows*0.2)))
        sns.set(font_scale=0.9)
        
        # Set color palette
        cmap = sns.light_palette("#009c8c", as_cmap=True)
        
        # plt.figure(figsize=(12, 9))  
    
        
        
        # Plot heatmap
        heatmap = sns.heatmap(
            df,
            annot=False,
            fmt=".3f",
            cmap=cmap,
            linewidths=0.5
        )
        
        heatmap.invert_yaxis()
        heatmap.xaxis.tick_top()
        heatmap.tick_params(axis='x', which='both', pad=10, top=True, bottom=False, length=0)
        
        heatmap.xaxis.set_label_position('top') 
        
        
        # Adjust color bar params
        cbar = heatmap.collections[0].colorbar
        cbar.ax.tick_params(axis='x', labelrotation=0) 
        cbar.ax.set_xlabel("VAF", fontsize=16, labelpad=10)
        
        # Adjust the position of the color bar label
        cbar.ax.xaxis.set_label_position('top')
    
        # Set axis labels and title
        plt.title(sample_id, fontsize=20, pad=25)
        plt.xlabel("Clones", fontsize=16, labelpad=15)
        plt.ylabel("Alterations", fontsize=16, labelpad=15)
        
        plt.savefig(f"{out_dir}/{sample_id}_VAF_heatmap.png", bbox_inches='tight')
        plt.close()



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--tree_df", action='store', required=True)
    input_parser.add_argument("--pvi_out", action='store', required=True)
    input_parser.add_argument("--mut_dir", action='store', required=True)
    input_parser.add_argument("--out_dir", action='store', required=True)

    args = input_parser.parse_args()


    heatmap_dict = build_heatmap_df(args.tree_df,args.pvi_out, args.mut_dir)
    plot_heatmaps(heatmap_dict,args.out_dir)