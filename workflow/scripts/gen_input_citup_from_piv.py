import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import os

freq_assignment_dict = {}

if __name__ == '__main__':

    frames = []
    df_frames = pd.DataFrame()

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--prevalence", action='store')
    input_parser.add_argument("--samples", action='store')
    input_parser.add_argument("--out_dir", action='store')
    args = input_parser.parse_args()
    prevalence_tsv = args.prevalence
    output_dir = args.out_dir

    l_freq = []

    # i.e.: /home/cleon/PycharmProjects/MSc-Project/out/pyclone-vi/output_933124.tsv
    df_prevalence = pd.read_csv(prevalence_tsv, sep='\t')

    # i.e.: 'UPN933124_454'
    samples = args.samples

    # array of unique samples
    arr_samples = np.array(df_prevalence['sample_id'].unique())

    # iterate by samples
    for idx_patient, row_sample in np.ndenumerate(arr_samples):
        # filter by sample_id
        df_sample = pd.DataFrame(df_prevalence.loc[df_prevalence['sample_id'] == row_sample])

        l_freq = []
        for index, row in df_sample.iterrows():
            sample_id = row['sample_id']
            vaf = row['cellular_prevalence']
            l_freq.append(vaf)

        freq_assignment_dict[row_sample] = l_freq

    if not Path(output_dir).exists():
        os.mkdir(output_dir)

    dict_df = pd.DataFrame({key: pd.Series(value) for key, value in freq_assignment_dict.items()})
    dict_df.to_csv(output_dir + '/freq_{}.txt'.format(samples), index=False, header=False, sep="\t")
