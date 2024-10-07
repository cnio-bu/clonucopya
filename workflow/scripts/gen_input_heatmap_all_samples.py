import pandas as pd
import argparse
import numpy as np
import os
import json


if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--locis", nargs='+', action='append')
    input_parser.add_argument("--pandrugs_dir", action='store')
    input_parser.add_argument("--out_dir", action='store')
    args = input_parser.parse_args()

    # i.e: /home/cleon/vep_data/18T59-I/tables/cluster.tsv /home/cleon/vep_data/16T320/tables/cluster.tsv
    locis_tsv = args.locis

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/pandrugs'
    pandrugs_dir = args.pandrugs_dir

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap'
    output_dir = args.out_dir

    vscore_tsv = pd.DataFrame()

    arr_locis_tsv = np.array(locis_tsv)
    for idx_loci_tsv, row_loci_tsv in np.ndenumerate(arr_locis_tsv):

        # loop for getting PyClone's clusters in order to know the number of sub-directories/clones
        # i.e: /home/cleon/vep_data/18T59-I/tables/loci.tsv
        loci_tsv = pd.read_csv(row_loci_tsv, sep='\t')
        arr = np.array(loci_tsv['cluster_id'].unique())
        for idx, x in np.ndenumerate(arr):
            df_loci = pd.DataFrame(loci_tsv.loc[loci_tsv['cluster_id'] == x])

            # recorrer por fila el fichero vscore (no sé si sacar el id de computación del tsv)
            sample_id = np.array(loci_tsv['sample_id'].unique())[0]
            cluster_id = x

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/16T320/0
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
                # componer columnas 'clone', y 'VAF' que se añadirá después al fichero de vscore
                vscore_loc = row_vscore['Loc']
                vscore_chr = row_vscore['Chr']
                loci_tsv_mutation_id = 'chr' + str(vscore_chr) + ':' + str(vscore_loc)
                df_loci_pos = pd.DataFrame(df_loci.loc[df_loci['mutation_id'] == loci_tsv_mutation_id])
                clone_column.append(df_loci_pos['cluster_id'].values[0])
                VAF_column.append(df_loci_pos['variant_allele_frequency'].values[0])

            vscore_tsv['clone'] = clone_column
            vscore_tsv['VAF'] = VAF_column

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/tumor
            if not os.path.exists(output_dir + '/' + str(sample_id)):
                os.makedirs(output_dir + '/' + str(sample_id))

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/tumor/0
            if not os.path.exists(output_dir + '/' + str(sample_id) + '/' + str(cluster_id)):
                os.makedirs(output_dir + '/' + str(sample_id) + '/' + str(cluster_id))

            vscore_tsv.to_csv(str(output_dir + '/{}/{}/{}-vscore.ann.tsv').format(sample_id, cluster_id, computation_id),
                              sep='\t',
                              index=False)

        # do the same for gene/drugs interaction
        # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap'
        output_dir = args.out_dir

        for idx, x in np.ndenumerate(arr):
            df_loci = pd.DataFrame(loci_tsv.loc[loci_tsv['cluster_id'] == x])

            # recorrer por fila el fichero vscore (no sé si sacar el id de computación del tsv)
            sample_id = np.array(loci_tsv['sample_id'].unique())[0]
            cluster_id = x

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/16T320/0
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
                    # componer columnas 'clone', y 'VAF' que se añadirá después al fichero de gene-drug
                    gd_gene_obj_json_str = str(row_gd['gene']).replace('[', '').replace(']', '').replace("'", "\"")
                    gd_gene_json_object = json.loads(gd_gene_obj_json_str)
                    gd_gene_symbol = gd_gene_json_object["geneSymbol"]

                    vscore_gene_hgnc = pd.DataFrame(vscore_ann_tsv.loc[vscore_ann_tsv['gene_hgnc'] == gd_gene_symbol])
                    gd_clone_column.append(np.array(vscore_gene_hgnc['clone'].unique())[0])
                    gd_VAF_column.append(np.array(vscore_gene_hgnc['VAF'].unique())[0])

                gen_drug_csv['clone'] = gd_clone_column
                gen_drug_csv['VAF'] = gd_VAF_column

                # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/tumor
                if not os.path.exists(output_dir + '/' + str(sample_id)):
                    os.makedirs(output_dir + '/' + str(sample_id))

                # i.e: /home/cleon/PycharmProjects/MSc-Project/out/heatmap/tumor/0
                if not os.path.exists(output_dir + '/' + str(sample_id) + '/' + str(cluster_id)):
                    os.makedirs(output_dir + '/' + str(sample_id) + '/' + str(cluster_id))

                gen_drug_csv.to_csv(str(output_dir + '/{}/{}/gene-drug.ann.tsv').format(sample_id, cluster_id),
                                    sep='\t',
                                    index=False)
