import pandas as pd
import argparse


def process_cna_facets(input_file, output_file, sex):
    """
    Process CNA data from facets and calculate copy numbers.
    
    Args:
        input_file (str): Path to input Facets CNA calling, VCF file (TSV) 
        output_file (str): Path to output TSV file
        sex (str): Sample sex ('male' or 'female')
    
    Returns:
        pandas.DataFrame: Processed data with copy numbers
    """
    try:
        cna_vcf = pd.read_csv(input_file, sep='\t', comment='#', header=None)
    except Exception as e:
        raise ValueError(f"Error reading VCF file: {e}")

    # Select atributes of interest
    try:
        cna_vcf_filt = cna_vcf.iloc[:,[0,1]].copy()
    except IndexError:
        raise ValueError("Incorrect cna dataframe format.")

    cna_vcf_filt.columns = ['Chrom', 'Start']
    info_column = cna_vcf.iloc[:,[7]].squeeze().copy()

    # Parse atributes
    cna_vcf_filt['End'] = info_column.str.split('END=').str[1].str.split(';').str[0]
    cna_vcf_filt['total_cn'] = info_column.str.split('TCN_EM=').str[1].str.split(';').str[0]
    cna_vcf_filt['minor_cn'] = info_column.str.split('LCN_EM=').str[1].str.split(';').str[0]

    # Drop missing minor_cn
    cna_vcf_filt = cna_vcf_filt[cna_vcf_filt.minor_cn != '.']

    cna_vcf_filt['End'] = cna_vcf_filt['End'].astype(int)
    cna_vcf_filt['total_cn'] = cna_vcf_filt['total_cn'].astype(int)
    cna_vcf_filt['minor_cn'] = cna_vcf_filt['minor_cn'].astype(int)

    # Calculate major_cn = total_cn - minor_cn
    cna_vcf_filt['major_cn'] = cna_vcf_filt.apply(lambda row: row['total_cn'] - row['minor_cn'] if pd.notnull(row['minor_cn']) else np.nan, axis=1)

    # Calculate normal_cn per cn
    cna_vcf_filt['normal_cn'] = cna_vcf_filt['Chrom'].apply(
    lambda chrom: 1 if (sex.lower() == 'male' and str(chrom) in ['chrX', 'chrY']) else 2)

    # Remove the total_cn column
    cna_vcf_filt =cna_vcf_filt.drop('total_cn', axis=1)

    # Reset the index if needed
    cna_vcf_filt = cna_vcf_filt.reset_index(drop=True)

    # Rearrange dataframe
    correct_order = ['Chrom', 'Start', 'End', 'major_cn', 'minor_cn', 'normal_cn']
    cna_vcf_filt = cna_vcf_filt[correct_order]

    cna_vcf_filt.to_csv(output_file, index=False, sep='\t')

    return cna_vcf_filt



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_file", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)
    input_parser.add_argument("--sex", action='store', required=True)

    args = input_parser.parse_args() 

    process_cna_facets(args.input_file, args.output_file, args.sex)
