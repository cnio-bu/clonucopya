import pandas as pd
import argparse
import numpy as np
import os
import json


if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--samples", nargs='+', action='append')
    input_parser.add_argument("--results_mutations", action='store')
    input_parser.add_argument("--pandrugs_dir", action='store')
    input_parser.add_argument("--out_dir", action='store')
    args = input_parser.parse_args()

    # i.e: /home/cleon/PycharmProjects/MSc-Project/out/mapscape/933124/results_mutations.tsv file
    results_mutations_tsv = pd.read_csv(args.results_mutations, sep=',')

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/pandrugs/933124'
    pandrugs_dir = args.pandrugs_dir

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap/933124'
    output_dir = args.out_dir

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
            # i.e: /home/cleon/vep_data/933124/tumor/ensemble_tumor.0.tsv
            sample_id = row_sample
            cluster_id = row_smp_clone

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/933124/tumor/0
            path_computation_tsv = pandrugs_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/computation.tsv'
            computation_tsv = pd.read_csv(path_computation_tsv, sep='\t')
            url_computation = ''
            computation_id = ''
            for idx_comp_tsv, row_comp_tsv in computation_tsv.iterrows():
                url_computation = row_comp_tsv['url_computation']
                url_computation_words = url_computation.split('/')
                computation_id = url_computation_words[len(url_computation_words) - 1]

            path_vscore_tsv = pandrugs_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/' \
                              + str(computation_id) + '-vscore.tsv'
            vscore_tsv = pd.read_csv(path_vscore_tsv, sep='\t')

            clone_column = []
            VAF_column = []

            for idx_vscore, row_vscore in vscore_tsv.iterrows():
                # compound columns 'clone', y 'VAF' for adding to vscore file later
                vscore_loc = row_vscore['Loc']
                vscore_chr = 'chr' + str(row_vscore['Chr'])
                df_results_mutations_pos = pd.DataFrame(results_mutations_tsv.loc[(results_mutations_tsv['chrom'] == vscore_chr) &
                                                                                  (results_mutations_tsv['coord'] == vscore_loc) &
                                                                                  (results_mutations_tsv['sample_id'] == sample_id)])

                clone_column.append(df_results_mutations_pos['clone_id'].values[0])
                VAF_column.append(df_results_mutations_pos['VAF'].values[0])

            vscore_tsv['clone'] = clone_column
            vscore_tsv['VAF'] = VAF_column

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/933124/tumor
            if not os.path.exists(output_dir + '/' + str(sample_id)):
                os.makedirs(output_dir + '/' + str(sample_id))

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/933124/tumor/0
            if not os.path.exists(output_dir + '/' + str(sample_id) + '/' + str(cluster_id)):
                os.makedirs(output_dir + '/' + str(sample_id) + '/' + str(cluster_id))

            vscore_tsv.fillna(0)
            vscore_tsv.to_csv(
                str(output_dir + '/{}/{}/{}-vscore.ann.tsv').format(sample_id, cluster_id, computation_id),
                sep='\t',
                index=False)

            # do the same for gene/drugs interaction
            # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap'
            output_dir = args.out_dir

            for idx_smp_clone_gd, row_smp_clone_gd in np.ndenumerate(arr_sample_clones):
                # i.e: /home/cleon/vep_data/933124/tumor/ensemble_tumor.0.tsv
                sample_id = row_sample
                cluster_id = row_smp_clone

                # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/933124/tumor/0
                path_computation_tsv = pandrugs_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/computation.tsv'
                computation_tsv = pd.read_csv(path_computation_tsv, sep='\t')
                url_computation = ''
                computation_id = ''
                for idx_comp_tsv, row_comp_tsv in computation_tsv.iterrows():
                    url_computation = row_comp_tsv['url_computation']
                    url_computation_words = url_computation.split('/')
                    computation_id = url_computation_words[len(url_computation_words) - 1]

                path_gen_drug_csv = pandrugs_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/' + 'gene-drug.csv'
                if os.path.exists(path_gen_drug_csv):
                    gen_drug_csv = pd.read_csv(path_gen_drug_csv, sep=',')

                    gd_clone_column = []
                    gd_VAF_column = []

                    path_vscore_ann_tsv = output_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/' \
                                          + str(computation_id) + '-vscore.ann.tsv'
                    vscore_ann_tsv = pd.read_csv(path_vscore_ann_tsv, sep='\t')

                    for idx_gd, row_gd in gen_drug_csv.iterrows():
                        # compound columns 'clone', y 'VAF' for adding to gene-drug file later
                        gd_gene_obj_json_str = str(row_gd['gene']).replace('[', '').replace(']', '').replace("'", "\"")
                        gd_gene_json_object = json.loads(gd_gene_obj_json_str)
                        gd_gene_symbol = gd_gene_json_object["geneSymbol"]

                        vscore_gene_hgnc = pd.DataFrame(vscore_ann_tsv.loc[vscore_ann_tsv['gene_hgnc'] == gd_gene_symbol])
                        gd_clone_column.append(np.array(vscore_gene_hgnc['clone'].unique())[0])
                        gd_VAF_column.append(np.array(vscore_gene_hgnc['VAF'].unique())[0])

                    gen_drug_csv['clone'] = gd_clone_column
                    gen_drug_csv['VAF'] = gd_VAF_column

                    # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/933124/tumor
                    if not os.path.exists(output_dir + '/' + str(sample_id)):
                        os.makedirs(output_dir + '/' + str(sample_id))

                    # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/933124/tumor/0
                    if not os.path.exists(output_dir + '/' + str(sample_id) + '/' + str(cluster_id)):
                        os.makedirs(output_dir + '/' + str(sample_id) + '/' + str(cluster_id))

                    gen_drug_csv.to_csv(str(output_dir + '/{}/{}/gene-drug.ann.tsv').format(sample_id, cluster_id),
                                        sep='\t',
                                        index=False)
