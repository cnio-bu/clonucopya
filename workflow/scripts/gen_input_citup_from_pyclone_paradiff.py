import pandas as pd
import numpy as np
import argparse

freq_assignment_dict = {}

if __name__ == '__main__':

    frames = []
    df_frames = pd.DataFrame()

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--loci", nargs='+', action='append')
    input_parser.add_argument("--out_dir", action='store')
    args = input_parser.parse_args()
    loci_tsv = args.loci

    output_dir = args.out_dir

    arr = np.array(loci_tsv)
    # iterate all input loci.tsv filenames
    for idx_loci_tsv, row_loci_tsv in np.ndenumerate(arr):
        l_freq = []
        df_loci = pd.read_csv(row_loci_tsv, sep='\t')

        for index, row in df_loci.iterrows():
            sample_id = row['sample_id']
            vaf = row['variant_allele_frequency']
            l_freq.append(vaf)

        freq_assignment_dict[sample_id] = l_freq

    print(freq_assignment_dict)
    dict_df = pd.DataFrame({key: pd.Series(value) for key, value in freq_assignment_dict.items()})
    dict_df.to_csv(output_dir + '/freq.txt', index=False, header=False, sep="\t")

