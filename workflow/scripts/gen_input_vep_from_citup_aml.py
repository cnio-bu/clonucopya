import numpy as np
import pandas as pd
import argparse
from pathlib import Path


if __name__ == '__main__':

    frames = []

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--samples", nargs='+', action='append')
    input_parser.add_argument("--mutations", action='store')
    input_parser.add_argument("--results_mutations", action='store')
    input_parser.add_argument("--results_clonal_prev", action='store')
    input_parser.add_argument("--vep_default_mutations", action='store')
    input_parser.add_argument("--out_dir", action='store')

    args = input_parser.parse_args()

    output_dir = args.out_dir

    # file containing all mutations
    mutations_tsv = pd.read_csv(args.mutations, sep='\t')

    # file containing all mutations in clusters
    # i.e.:
    results_mutations_tsv = pd.read_csv(args.results_mutations, sep=',')

    # ['tumor', 'relapse']
    arr_samples = np.array(args.samples)
    # loop for samples in the results_mutations.tsv file
    for idx_sample, row_sample in np.ndenumerate(arr_samples):
        # filter by sample
        df_sample = pd.DataFrame(results_mutations_tsv.loc[results_mutations_tsv['sample_id'] == row_sample])

        # array of unique clones of the sample
        arr_sample_clones = np.array(df_sample['clone_id'].unique())

        # loop for sample's clones in the results_mutations.tsv file
        for idx_smp_clone, row_smp_clone in np.ndenumerate(arr_sample_clones):
            # filter by clone
            df_sample_clone = pd.DataFrame(results_mutations_tsv.loc[results_mutations_tsv['clone_id'] == row_smp_clone])
            frames = []

            # iterate each row for sample's clone
            for index, row_smp in df_sample_clone.iterrows():
                chro = row_smp['chrom']
                pos = row_smp['coord']

                mutation_tsv = pd.DataFrame(
                    mutations_tsv.loc[mutations_tsv['case'].apply(lambda y: y in arr_samples)])
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

            if not Path(output_dir).exists():
                Path(output_dir).mkdir()

            if not Path(output_dir + '/' + row_sample).exists():
                Path(output_dir + '/' + row_sample).mkdir()
            output_file = output_dir + '/' + row_sample + '/ensembl_{}.{}.tsv'.format(row_sample, row_smp_clone)

            if len(frames) > 0:
                df_frames = pd.concat(frames)
                df_frames.to_csv(output_file,
                                 sep='\t',
                                 index=False,
                                 header=False,
                                 mode='w')
