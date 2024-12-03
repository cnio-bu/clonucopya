import pandas as pd
import argparse

# i.e.
# sex = "female"
# input_file = 'TCGA-LAML_allele_cnv_ascat3.tsv'
# output_file = 'CLEAN_TEST.tsv'


def calculate_copy_numbers(total_cn, chromosome, sex):
    """
    Calculate copy numbers based on total CN, chromosome, and sex.
    
    Args:
        total_cn (int): Total copy number
        chromosome (str): Chromosome name
        sex (str): Sample sex ('male' or 'female')
    
    Returns:
        tuple: (major_cn, minor_cn, normal_cn)
    """
    # Determine normal_cn based on chromosome
    if sex.lower() == 'male' and chromosome in ['chrX', 'chrY']:
        normal_cn = 1  # cn_normal value for sexual chromosome in male
    else:
        normal_cn = 2

    # Calculate major_cn and minor_cn
    if total_cn % 2 == 0:
        major_cn = minor_cn = total_cn // 2
    else:
        major_cn = (total_cn + 1) // 2
        minor_cn = (total_cn - 1) // 2

    return major_cn, minor_cn, normal_cn

def process_data(input_file, output_file, sex):
    """
    Process CNV data and calculate copy numbers.
    
    Args:
        input_file (str): Path to input TSV file
        output_file (str): Path to output TSV file
        sex (str): Sample sex
    
    Returns:
        pandas.DataFrame: Processed data with copy numbers
    """
    # Read the input file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')
    cnv_df = df.sort_values(by=['Chrom', 'Start', 'End'])
    
    # Calculate copy numbers for each row
    cnv_df['major_cn'], cnv_df['minor_cn'], cnv_df['normal_cn'] = zip(*cnv_df.apply(
        lambda row: calculate_copy_numbers(int(row['value']), row['Chrom'], sex), 
        axis=1
    ))
    
    # Select columns to save
    columns_to_save = ['Chrom', 'Start', 'End', 'value', 'major_cn', 'minor_cn', 
                      'normal_cn', 'sample']
    df_to_save = cnv_df[columns_to_save]

    # Save the selected columns to a TSV file
    df_to_save.to_csv(output_file, index=False, sep='\t')

    return df_to_save



if __name__ == '__main__':
    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_file", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)
    input_parser.add_argument("--sex", action='store', required=True)

    args = input_parser.parse_args() 

    # Process the data
    process_data(args.input_file, args.output_file, args.sex)
