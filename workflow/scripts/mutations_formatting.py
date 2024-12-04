import pandas as pd
import argparse

# i.e.
# input_vcf = '01_B2237D1M1.final.vcf'
# output_file = 'Mutation_CLEAN_TEST.tsv'


    
    
def process_vcf_mutations(input_vcf, output_file):
    """
    Process VCF file to extract mutation information and read counts.
    
    Args:
        input_vcf (str): Path to input VCF file
        output_file (str, optional): Path to output TSV file. If None, won't save to file
        
    Returns:
        pandas.DataFrame: Processed mutations data
    """
    # Read VCF file
    try:
        mut_vcf = pd.read_csv(input_vcf, sep='\t', comment='#', header=None)
    except Exception as e:
        raise ValueError(f"Error reading VCF file: {e}")
    
    # Extract relevant columns
    mut_vcf_filt = mut_vcf.iloc[:,[0,1,3,4]].copy()
    mut_vcf_filt.columns = ['CHROM', 'POS', 'REF', 'ALT']
    
    # Extract genotype information
    genotype_column = mut_vcf.iloc[:,[9]].squeeze().copy()
    
    # Parse read counts
    try:
        mut_vcf_filt['ref_counts'] = genotype_column.str.split(':').str[1].str.split(',').str[0]
        mut_vcf_filt['alt_counts'] = genotype_column.str.split(':').str[1].str.split(',').str[1]
        
        # Convert counts to numeric
        mut_vcf_filt['ref_counts'] = pd.to_numeric(mut_vcf_filt['ref_counts'], errors='coerce')
        mut_vcf_filt['alt_counts'] = pd.to_numeric(mut_vcf_filt['alt_counts'], errors='coerce')
    except Exception as e:
        raise ValueError(f"Error parsing read counts: {e}")
    
    # Create mutation ID
    mut_vcf_filt['mutation_id'] = (mut_vcf_filt['CHROM'] + ':' + 
                                  mut_vcf_filt['POS'].astype(str) + ':' + 
                                  mut_vcf_filt['REF'] + ':' + 
                                  mut_vcf_filt['ALT'])
    
    # Save to file if output_file is specified
    if output_file:
        mut_vcf_filt.to_csv(output_file, sep='\t', index=True)
    
    return mut_vcf_filt


    
if __name__ == '__main__':
    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_vcf", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)
    args = input_parser.parse_args() 
    
    # Process the vcf file
    process_vcf_mutations(args.input_vcf, args.output_file)