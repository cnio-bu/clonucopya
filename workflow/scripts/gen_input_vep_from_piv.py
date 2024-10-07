import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import os


if __name__ == '__main__':

    frames = []

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--mutations", action='store')
    input_parser.add_argument("--prevalence", action='store')
    input_parser.add_argument("--out_dir", action='store')
    input_parser.add_argument("--patient", action='store')
    args = input_parser.parse_args()

    # i.e.: /home/cleon/PycharmProjects/MSc-Project/out/pyclone-vi/output_933124.tsv
    prevalence_tsv = pd.read_csv(args.prevalence, sep='\t')

    # i.e.: /home/cleon/vep_data/933124
    output_dir = args.out_dir

    # i.e.: '933124'
    patient_id = args.patient

    # file containing all mutations
    # i.e.: /home/cleon/TFM/datos_AML/mutations/mutations.tsv
    mutations_tsv = pd.read_csv(args.mutations, sep='\t')

    # ['tumor', 'relapse']
    # array of unique sample_ids
    arr_samples = np.array(prevalence_tsv['sample_id'].unique())

    # loop for samples in the prevalence_tsv file
    for idx_sample, row_sample in np.ndenumerate(arr_samples):
        # print('row_sample={}'.format(row_sample))

        # filter by sample
        df_sample = pd.DataFrame(prevalence_tsv.loc[prevalence_tsv['sample_id'] == row_sample])

        # array of unique clones of the sample
        arr_sample_clones = np.array(df_sample['cluster_id'].unique())
        # print('arr_sample_clones={}'.format(arr_sample_clones))

        # loop for sample's clones in the prevalence_tsv file
        for idx_smp_clone, row_smp_clone in np.ndenumerate(arr_sample_clones):
            # i.e: /home/cleon/vep_data/933124/tumor/ensemble_tumor.0.tsv
            sample_id = row_sample
            cluster_id = row_smp_clone
            # print('sample_id={} and cluster_id{}'.format(sample_id, cluster_id))

            # get rows by assigned clone
            df_sample_clone = pd.DataFrame(df_sample.loc[df_sample['cluster_id'] == cluster_id])
            # print('df_sample_clone={}'.format(df_sample_clone))

            frames = []

            for index, row in df_sample_clone.iterrows():
                mutation_id = row['mutation_id']
                chro = mutation_id.split(':')[0]
                pos = mutation_id.split(':')[1]
                # print('chro={} and pos={}'.format(chro, pos))

                # get mutations by assigned case
                mutation_case_tsv = pd.DataFrame(mutations_tsv.loc[mutations_tsv['case'] == sample_id])
                # print('mutation_case_tsv={}'.format(mutation_case_tsv))
                for idx_mut, row_mut in mutation_case_tsv.iterrows():

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

            if not Path(output_dir).exists():
                os.mkdir(output_dir)
                os.chmod(output_dir, 0o777)

            if not Path(output_dir + '/' + sample_id).exists():
                os.mkdir(output_dir + '/' + sample_id)
                os.chmod(output_dir + '/' + sample_id, 0o777)

            output_file = output_dir + '/' + sample_id + '/ensembl_{}.{}.tsv'.format(sample_id, cluster_id)

            if len(frames) > 0:
                # print(output_file)
                df_frames = pd.concat(frames)
                df_frames.to_csv(output_file,
                                 sep='\t',
                                 index=False,
                                 header=False,
                                 mode='w')
            else:
                print('frames has no elements')
