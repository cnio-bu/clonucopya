import pandas as pd
import argparse


def process_data(input_file, output_file):
    """
    Process CNV data from battemberg and calculate copy numbers.
    
    Args:
        input_file (str): Path to input txt file
        output_file (str): Path to output TSV file
    
    Returns:
        pandas.DataFrame: Processed data with copy numbers
    """
    # Load copy number output from battemberg
    df = pd.read_csv(input_file, sep='\t')

    # Select atributes of interest
    try:
        cnv_filt = df[['chr', 'startpos', 'endpos','nMaj1_A','nMin1_A']].copy()
        cnv_filt
    except IndexError:
        raise ValueError("Incorrect cnv dataframe format.")

    cnv_filt.columns = ['Chrom','Start', 'End', 'major_cn', 'minor_cn']

    cnv_filt['major_cn'] = cnv_filt['major_cn'].astype(int)
    cnv_filt['minor_cn'] = cnv_filt['minor_cn'].astype(int)
    cnv_filt['normal_cn'] = 2

    cnv_filt.to_csv(output_file, index=False, sep='\t')

    return cnv_filt



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_file", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args() 

    process_data(args.input_file, args.output_file)
