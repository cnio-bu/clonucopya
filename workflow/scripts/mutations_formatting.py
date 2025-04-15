import pandas as pd
import argparse



    
    
def process_vcf_mutations(input_vcf, output_file):
    """
    Process VCF file to extract mutation information and read counts.
    
    Args:
        input_vcf (str): Path to input VCF file
        output_file (str, optional): Path to output TSV file. If None, won't save to file
        
    Returns:
        pandas.DataFrame: Processed mutations data
    """
    # Load mutations
    try:
        mut_vcf = pd.read_csv(input_vcf, sep='\t', comment='#', header=None)
    except Exception as e:
        raise ValueError(f"Error reading VCF file: {e}")
    
    # Select relevant columns
    mut_vcf_filt = mut_vcf.iloc[:,[0,1,3,4]].copy()
    mut_vcf_filt.columns = ['CHROM', 'POS', 'REF', 'ALT']
    
    # Adapt format to VEP Standard
    mut_vcf_filt['REF'] = mut_vcf_filt['REF'].apply(lambda x: x.replace('.', '-') if '.' in x else x)
    mut_vcf_filt['ALT'] = mut_vcf_filt['ALT'].apply(lambda x: x.replace('.', '-') if '.' in x else x)

    # Extract genotype info
    genotype_column = mut_vcf.iloc[:,[-1]].squeeze().copy()

    # Parse read counts and variant alelle frequencies
    try:
        mut_vcf_filt['ref_counts'] = genotype_column.str.split(':').str[1].str.split(',').str[0]
        mut_vcf_filt['alt_counts'] = genotype_column.str.split(':').str[1].str.split(',').str[1]
        mut_vcf_filt['VAF'] = genotype_column.str.split(':').str[2]

        mut_vcf_filt['ref_counts'] = pd.to_numeric(mut_vcf_filt['ref_counts'], errors='coerce')
        mut_vcf_filt['alt_counts'] = pd.to_numeric(mut_vcf_filt['alt_counts'], errors='coerce')
        mut_vcf_filt['VAF'] = pd.to_numeric(mut_vcf_filt['VAF'], errors='coerce')
    except Exception as e:
        raise ValueError(f"Error parsing read counts: {e}")
    
    # Build mutation ID
    mut_vcf_filt['mutation_id'] = (mut_vcf_filt['CHROM'] + ':' + 
                                  mut_vcf_filt['POS'].astype(str) + ':' + 
                                  mut_vcf_filt['REF'] + ':' + 
                                  mut_vcf_filt['ALT'])

    if output_file:
        mut_vcf_filt.to_csv(output_file, sep='\t', index=False)

    return mut_vcf_filt


    
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_vcf", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)
    args = input_parser.parse_args() 
    
    process_vcf_mutations(args.input_vcf, args.output_file)
