import numpy as np
import pandas as pd
import tables
import argparse
import os
import h5py

if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_dir", action='store')
    input_parser.add_argument("--in_file", nargs='+', action='append')
    input_parser.add_argument("--out_dir", action='store')
    input_parser.add_argument("--prevalence", action='store')
    args = input_parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.out_dir
    in_file_h5 = args.in_file
    arr_h5 = np.array(in_file_h5)
    prevalence_tsv = args.prevalence

    # iterate all input h5 filenames
    for idx_h5, row_h5 in np.ndenumerate(arr_h5):

        file_name = row_h5

        os.chmod(file_name, 0o777)

        # h5 file
        f = tables.open_file(file_name, mode='r')
        # f = h5py.File(file_name, 'r')

        table1 = f.root.results.optimal
        # select the best tree
        best = list(np.array(table1.index))[list(np.array(table1.values)).index(True)]
        print(best)
        array = np.array(f.root.trees._f_get_child(str(best)).adjacency_list.block0_values)
        np.savetxt(file_name.replace(".h5", "_tree.csv"), array, delimiter=',', newline='\n', header='source,target',
                   fmt='%d', comments='')

        # for mapscape
        clones = np.array(f.root.trees._f_get_child(str(best)).clone_freq.block0_items)
        freq = np.array(f.root.trees._f_get_child(str(best)).clone_freq.block0_values)

        sample_list = []
        for x in np.array(freq):
            sample_list.append(x)

        clone_str = pd.DataFrame(
            {'sample_id': np.repeat(np.array(f.root.trees._f_get_child(str(best)).clone_freq.axis1), len(clones)),
             'clone_id': np.tile(clones, len(sample_list)),
             'clonal_prev': np.concatenate(freq, axis=0)}, columns=['sample_id', 'clone_id', 'clonal_prev'])

        clone_str.to_csv(file_name.replace(".h5", "_clonal_prev.csv"), index=False, header=True, sep=",")
        # get mutations for both timescape and mapscape
        variant_assignment = pd.DataFrame(
            {'index': np.array(f.root.trees._f_get_child(str(best)).variant_assignment.index),
             'values': np.array(f.root.trees._f_get_child(str(best)).variant_assignment.values)},
            columns=['index', 'values'])

        variant_assignment.to_csv(file_name.replace(".h5", "_variant_assignment.csv"), index=False, header=True,
                                  sep=",")

        va_values = np.array(f.root.trees._f_get_child(str(best)).variant_assignment.values)

        df_frames = []
        df_mutations = pd.DataFrame()

        arr_prevalence = np.array(prevalence_tsv)
        # iterate all input prevalence tsv filenames (/out/pyclone-vi/upn933124/output_upn933124.tsv)
        for idx_prev_tsv, row_prev_tsv in np.ndenumerate(arr_prevalence):

            mutations = pd.read_csv(row_prev_tsv, sep="\t", header=0, low_memory=False)
            for i in variant_assignment["index"].tolist():
                mutations.loc[mutations['cluster_id'] == i, 'cluster_citup'] = str(
                    variant_assignment[variant_assignment["index"] == i]["values"].tolist()[0])

            df_frames.append(pd.DataFrame(
                {'chrom': mutations["mutation_id"].str.split(':', expand=True)[0],
                 'coord': mutations["mutation_id"].str.split(':', expand=True)[1],
                 'clone_id': mutations["cluster_citup"],
                 'sample_id': mutations["sample_id"],
                 'VAF': mutations["cellular_prevalence"]},
                columns=['chrom', 'coord', 'clone_id', 'sample_id', 'VAF']))

        # append all dataframes of mutations
        for df_frame in df_frames:
            df_mutations.append(df_frame)

        df_mutations = pd.concat(df_frames, ignore_index=True)

        # export to csv
        df_mutations.to_csv(file_name.replace(".h5", "_mutations.csv"), index=False, header=True, sep=",")

    f.close()
