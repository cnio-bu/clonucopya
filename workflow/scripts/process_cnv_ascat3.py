import pandas as pd
import argparse


def process_data(input_file, output_file, sex):
    """
    Process CNV data from ascat3 and calculate copy numbers.
    
    Args:
        input_file (str): Path to input TSV file
        output_file (str): Path to output TSV file
        sex (str): Sample sex ('male' or 'female')
    
    Returns:
        pandas.DataFrame: Processed data with copy numbers
    """
    # Load copy number output from ascat3
    df = pd.read_csv(input_file, sep='\t', header=None )

    # Select atributes of interest
    try:
        cnv_filt = df.iloc[:,[6,7,8,4,5]].copy()
    except IndexError:
        raise ValueError("Incorrect cnv dataframe format.")

    cnv_filt.columns = ['Chrom','Start', 'End', 'major_cn', 'minor_cn']

    cnv_filt['major_cn'] = cnv_filt['major_cn'].astype(int)
    cnv_filt['minor_cn'] = cnv_filt['minor_cn'].astype(int)

    # Calculate normal_cn per cn
    cnv_filt['normal_cn'] = cnv_filt['Chrom'].apply(
    lambda chrom: 1 if (sex.lower() == 'male' and str(chrom) in ['chrX', 'chrY']) else 2)

    cnv_filt.to_csv(output_file, index=False, sep='\t')

    return cnv_filt



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_file", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)
    input_parser.add_argument("--sex", action='store', required=True)

    args = input_parser.parse_args() 

    process_data(args.input_file, args.output_file, args.sex)
