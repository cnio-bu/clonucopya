# It is also valid for the AML's dataset single-sample
import numpy as np
import pandas as pd
import argparse
from pathlib import Path


if __name__ == '__main__':

    frames = []

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--mutations", action='store')
    input_parser.add_argument("--loci", nargs='+', action='append')
    input_parser.add_argument("--out_dir", action='store')
    input_parser.add_argument("--samples", action='store')
    args = input_parser.parse_args()
    loci_tsv = args.loci
    samples_tsv = pd.read_csv(args.samples, sep='\t')

    output_dir = args.out_dir

    # file containing all mutations
    mutations_tsv = pd.read_csv(args.mutations, sep='\t')

    # array of unique patients
    arr_patients = np.array(samples_tsv['patient'].unique())

    arr = np.array(loci_tsv)
    # iterate by loci.tsv of its regarding case
    for idx_loci_tsv, row_loci_tsv in np.ndenumerate(arr):

        # file containing all mutations group by cluster_id
        df_clusters = pd.read_csv(row_loci_tsv, sep='\t')

        # i.e: 14T126
        sample_id = df_clusters['sample_id'].iat[0]

        arr_cases_by_patient = []

        # TODO: limitation. Currently, it is only working with one patient.
        for idx_patient, row_patient in np.ndenumerate(arr_patients):
            # filter by Patient
            df_patient = pd.DataFrame(samples_tsv.loc[samples_tsv['patient'] == row_patient])
            # array of unique cases by patient
            arr_cases_by_patient = np.array(df_patient['case'].unique())

        # get all cluster ids of the case
        arr = np.array(df_clusters['cluster_id'].unique())
        for idx, x in np.ndenumerate(arr):
            df_cluster = pd.DataFrame(df_clusters.loc[df_clusters['cluster_id'] == x])
            frames = []

            for index, row in df_cluster.iterrows():
                mutation_id = row['mutation_id']
                chro = mutation_id.split(':')[0]
                pos = mutation_id.split(':')[1]

                mutation_tsv = pd.DataFrame(
                    mutations_tsv.loc[mutations_tsv['case'].apply(lambda y: y in arr_cases_by_patient)])
                for idx_mut, row_mut in mutation_tsv.iterrows():

                    if row_mut['loc'] == int(pos):
                        mut = row_mut['mut']
                        mut_split = str(mut).split('/')
                        if mut_split[0] == '-':
                            # insertion
                            df2 = pd.DataFrame(np.array([[chro, pos, int(pos) - 1, mut]]),
                                               columns=['chro', 'start', 'end', 'mut'])
                            frames.append(df2)
                        elif mut_split[1] == '-':
                            # deletion
                            df2 = pd.DataFrame(np.array([[chro, pos, int(pos) + 1, mut]]),
                                               columns=['chro', 'start', 'end', 'mut'])
                            frames.append(df2)
                        else:
                            # substitution
                            df2 = pd.DataFrame(np.array([[chro, pos, pos, mut]]),
                                               columns=['chro', 'start', 'end', 'mut'])
                            frames.append(df2)

            if not Path(output_dir + '/' + sample_id).exists():
                Path(output_dir + '/' + sample_id).mkdir()
            output_file = output_dir + '/' + sample_id + '/ensembl_{}.{}.tsv'.format(sample_id, x)

            if len(frames) > 0:
                df_frames = pd.concat(frames)
                df_frames.to_csv(output_file,
                                 sep='\t',
                                 index=False,
                                 header=False,
                                 mode='w')
