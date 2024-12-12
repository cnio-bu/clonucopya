import pandas as pd
import argparse


import pandas as pd
import argparse


def prep_pvi_to_vep(mutations, pvi_file, sample_id, out_dir):
    """
    Format output from PyClone-VI to vcf standard as input 
    for ensembl-vep
    Args:
        mutations (str): Path to input TSV file
        pvi_file (str): Path to input TSV file 
        sample_id (str): Name of the sample 
        out_dir (str): Path to output directory
    Returns:
        dict: Dictionary containing DataFrames for each cluster
    """
    # Load mutation's vcf file and pyclone-vi output
    try:
        mutation_vcf = pd.read_csv(mutations,
                     sep='\t', comment='#', header=None)
    except Exception as e:
        raise ValueError(f"Error reading VCF file: {e}")

    try:
        pvi_data = pd.read_csv(pvi_file, sep='\t')
    except Exception as e:
        raise ValueError(f"Error reading TSV file: {e}")

    sample = sample_id
    
    # List clone clusters 
    unique_clusters = pvi_data['cluster_id'].unique()
    
    # Dictionary to store DataFrames for each cluster
    cluster_dataframes = {}

    try:
        # Process each cluster
        for cluster in unique_clusters:
            # Create empty dataframe for this cluster
            mut_vcf_vep = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
            
            # Process each mutation in cluster
            cluster_df = pvi_data[pvi_data['cluster_id'] == cluster]
            
            for _, row in cluster_df.iterrows():
                mutation_parts = row['mutation_id'].split(':')
                chrom = mutation_parts[0]
                pos = int(mutation_parts[1])
                
                matching_vcf_rows = mutation_vcf[
                    (mutation_vcf[0] == chrom) & 
                    (mutation_vcf[1] == pos)
                ]
                if not matching_vcf_rows.empty:
                    vcf_row = matching_vcf_rows.iloc[0]
                    
                    new_row = pd.DataFrame({
                        "#CHROM": [vcf_row[0]],
                        "POS": [vcf_row[1]],
                        "ID": [row['mutation_id']],
                        "REF": [vcf_row[3]],
                        "ALT": [vcf_row[4]],
                        "QUAL": [vcf_row[5]],
                        "FILTER": [vcf_row[6]],
                        "INFO": [vcf_row[7]],
                        "FORMAT": [vcf_row[8]]
                    })
                
                    mut_vcf_vep = pd.concat([mut_vcf_vep, new_row], ignore_index=True)
            
            # Save the cluster DataFrame
            mut_vcf_vep.to_csv(f'{out_dir}/{sample}_cluster_{cluster}.tsv', sep='\t', index=False)
            cluster_dataframes[cluster] = mut_vcf_vep

        return cluster_dataframes

    except Exception as e:
        raise ValueError(f"Error creating formatting p-vi to vep format: {e}")




    
if __name__ == '__main__':
    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--mutations", action='store', required=True)
    input_parser.add_argument("--pvi_output", action='store', required=True)
    input_parser.add_argument("--sample_id", action='store', required=True)
    input_parser.add_argument("--out_dir", action='store', required=True)
    args = input_parser.parse_args() 
    
    # Process the vcf file
    prep_pvi_to_vep(args.mutations, args.pvi_output, args.sample_id, args.out_dir)
