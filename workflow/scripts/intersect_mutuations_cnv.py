import pandas as pd
import argparse



# i.e.
# mutations_file = 'mutation_B2237D1M1_formatted_111124.tsv'
# cnv_file = 'TCGA-LAML_allele_cnv_formatted_SAMPLEid.tsv'
# output_file = 'TCGA-LAML_allele_cnv_formatted_SAMPLEid.tsv'



def create_pyclone_vi_input(mutations_file, cnv_file, output_file):
    # Load mutations and CNV data
    mutations = pd.read_csv(mutations_file, sep='\t')
    cnv = pd.read_csv(cnv_file, sep='\t')
    
    result = []
    
    # Iterate through mutations
    for _, mutation in mutations.iterrows():
        chrom = mutation['CHROM']
        pos = mutation['POS']
        
        # Find matching CNV segment
        for _, cnv_row in cnv[cnv['Chrom'] == chrom].iterrows():
            start = cnv_row['Start']
            end = cnv_row['End']
            
            # Only append if there is a match
            if start <= pos <= end:
                result.append({
                    'mutation_id': mutation['mutation_id'],
                    'sample_id': cnv_row['sample'],
                    'ref_counts': mutation['ref_counts'],
                    'alt_counts': mutation['alt_counts'],
                    'major_cn': int(cnv_row['major_cn']),
                    'minor_cn': int(cnv_row['minor_cn']),
                    'normal_cn': int(cnv_row['normal_cn'])
                })
                break
    
    # Create DataFrame and save to file
    output = pd.DataFrame(result)
    output.to_csv(output_file, sep='\t', index=False)



if __name__ == '__main__':
    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--mutations", action='store', required=True)
    input_parser.add_argument("--cnvs", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args() 

    # Process the data
    create_pyclone_vi_input(args.mutations, args.cnvs, args.output_file)
