import pandas as pd
import argparse

def create_pyclone_vi_input(sample_id, mutations_file, cnv_file, output_file):
    """
    Intersect mutations with copy number regions. Copy number must contain mutation to be valid.

    Args:
        sample_id (str): Sample ID of the pair mutation-cnv files
        mutations_file (str): Path to mutation file of the processed sample
        cnv_file (str): Path to copy number file of the processed sample
        output_file (str, optional): Path to output TSV file.

    Returns:
        pandas.DataFrame: Intersected sample's mutations-cnvs pair
    """

    # Load mutations and copy number files
    mutations = pd.read_csv(mutations_file, sep='\t')
    cnv = pd.read_csv(cnv_file, sep='\t')

    result = []
    # Iterate through mutations
    for _, mut_row in mutations.iterrows(): 
        chr_mut = mut_row['CHROM']  
        pos_mut = mut_row['POS']
        
        # Iterate through CNV segments
        for _, cnv_row in cnv.iterrows():
            chr_cnv = cnv_row['Chrom']
            start = cnv_row['Start']
            end = cnv_row['End']
            
            if (chr_mut == chr_cnv) and (start <= pos_mut <= end):
                result.append({
                    'mutation_id': mut_row['mutation_id'],
                    'sample_id': sample_id,
                    'ref_counts': mut_row['ref_counts'], 
                    'alt_counts': mut_row['alt_counts'],
                    'major_cn': int(cnv_row['major_cn']),
                    'minor_cn': int(cnv_row['minor_cn']),
                    'normal_cn': int(cnv_row['normal_cn'])
                })

    
    # Convert dictionary to dataframe and save it
    output = pd.DataFrame(result)
    output.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--sample_id", action='store', required=True)
    input_parser.add_argument("--mutations", action='store', required=True)
    input_parser.add_argument("--cnvs", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args()


    create_pyclone_vi_input(args.sample_id,args.mutations, args.cnvs, args.output_file)
