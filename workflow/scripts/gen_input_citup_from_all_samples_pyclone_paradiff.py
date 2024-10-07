import pandas as pd
import numpy as np
import argparse

freq_assignment_dict = {}

if __name__ == '__main__':

    frames = []
    df_frames = pd.DataFrame()

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--loci", action='store')
    input_parser.add_argument("--out_dir", action='store')
    args = input_parser.parse_args()
    loci_tsv = args.loci
    output_dir = args.out_dir

    l_freq = []
    df_loci = pd.read_csv(loci_tsv, sep='\t')

    # array of unique samples
    arr_samples = np.array(df_loci['sample_id'].unique())

    # iterate by samples
    for idx_patient, row_sample in np.ndenumerate(arr_samples):
        # filter by sample_id
        df_sample = pd.DataFrame(df_loci.loc[df_loci['sample_id'] == row_sample])

        l_freq = []
        for index, row in df_sample.iterrows():
            sample_id = row['sample_id']
            vaf = row['variant_allele_frequency']
            l_freq.append(vaf)

        freq_assignment_dict[row_sample] = l_freq

    # print(freq_assignment_dict)
    dict_df = pd.DataFrame({key: pd.Series(value) for key, value in freq_assignment_dict.items()})
    dict_df.to_csv(output_dir + '/freq.txt', index=False, header=False, sep="\t")

