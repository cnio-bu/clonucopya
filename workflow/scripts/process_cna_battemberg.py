import pandas as pd
import argparse


def process_data(input_file, output_file):
    """
    Process CNA data from battemberg and calculate copy numbers.
    
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
        cna_filt = df[['chr', 'startpos', 'endpos','nMaj1_A','nMin1_A']].copy()
        cna_filt
    except IndexError:
        raise ValueError("Incorrect cna dataframe format.")

    cna_filt.columns = ['Chrom','Start', 'End', 'major_cn', 'minor_cn']

    cna_filt['major_cn'] = cna_filt['major_cn'].astype(int)
    cna_filt['minor_cn'] = cna_filt['minor_cn'].astype(int)
    cna_filt['normal_cn'] = 2

    cna_filt.to_csv(output_file, index=False, sep='\t')

    return cna_filt



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_file", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args() 

    process_data(args.input_file, args.output_file)
