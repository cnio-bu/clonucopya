import docker
from pathlib import Path
import os
import pandas as pd
import numpy as np
import argparse


if __name__ == '__main__':

    client = docker.from_env()

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--clusters", nargs='+', action='append')
    input_parser.add_argument("--patient", action='store')
    input_parser.add_argument("--vep_dir", action='store')
    args = input_parser.parse_args()

    # i.e: /home/cleon/vep_data/18T59-I/tables/cluster.tsv /home/cleon/vep_data/16T320/tables/cluster.tsv
    clusters_tsv = args.clusters

    # i.e: '/home/cleon/vep_data'
    vep_dir = args.vep_dir

    arr = np.array(clusters_tsv)
    for idx_cluster_tsv, row_cluster_tsv in np.ndenumerate(arr):

        # i.e: /home/cleon/vep_data/18T59-I/tables/cluster.tsv
        cluster_tsv = pd.read_csv(row_cluster_tsv, sep='\t')
        arr = np.array(cluster_tsv['cluster_id'].unique())
        for idx, x in np.ndenumerate(arr):
            df_cluster = pd.DataFrame(cluster_tsv.loc[cluster_tsv['cluster_id'] == x])

            for index, row in df_cluster.iterrows():
                sample_id = row['sample_id']
                cluster_id = row['cluster_id']

                os.chmod(vep_dir + '/' + sample_id, 0o777)

                # i.e: ensembl_16T320.0.tsv
                vep_file = 'ensembl_{}.{}.tsv'.format(sample_id, cluster_id)
                # print('>>>>>>>>>>>>>>{}'.format(vep_file))

                # i.e: /opt/vep/.vep/18T59-I/ensembl_16T320.0.tsv
                input_ensembl = "/opt/vep/.vep/" + sample_id + "/" + vep_file

                # i.e: /opt/vep/.vep/18T59-I/ensembl_16T320.0.vcf
                Path(vep_dir + '/' + sample_id + '/' + vep_file.replace('.tsv', '.vcf')).write_text("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO")
                os.chmod(vep_dir + '/' + sample_id + '/' + vep_file.replace('.tsv', '.vcf'), 0o777)
                output_ensembl = "/opt/vep/.vep/" + sample_id + '/' + vep_file.replace('.tsv', '.vcf')

                # run docker
                client.containers.run("ensemblorg/ensembl-vep:release_104.3",
                                      "./vep --cache --cache_version 104 --offline --format ensembl --vcf --force_overwrite -dir_cache /opt/vep/.vep --input_file {} --output_file {} --assembly GRCh38".format(input_ensembl, output_ensembl),
                                      stdout=True,
                                      stderr=False,
                                      remove=True,
                                      volumes=[f'{vep_dir}:/opt/vep/.vep'])

    # do it also for multi-sample (patient)
    # i.e: /home/cleon/vep_data/47/tables/cluster.tsv
    path_patient_cluster_tsv = vep_dir + '/' + args.patient + '/tables/cluster.tsv'
    df_patient_cluster_tsv = pd.read_csv(path_patient_cluster_tsv, sep='\t')

    for index_pc, row_pc in df_patient_cluster_tsv.iterrows():
        sample_id = row_pc['sample_id']
        cluster_id = row_pc['cluster_id']

        os.chmod(vep_dir + '/' + sample_id, 0o777)

        # i.e: ensembl_47.0.tsv
        vep_file = 'ensembl_{}.{}.tsv'.format(sample_id, cluster_id)
        # print('>>>>>>>>>>>>>>{}'.format(vep_file))

        # i.e: /opt/vep/.vep/47/ensembl_47.0.tsv
        input_ensembl = "/opt/vep/.vep/" + sample_id + "/" + vep_file

        # i.e: /opt/vep/.vep/47/ensembl_47.0.vcf
        Path(vep_dir + '/' + sample_id + '/' + vep_file.replace('.tsv', '.vcf')).write_text(
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO")
        os.chmod(vep_dir + '/' + sample_id + '/' + vep_file.replace('.tsv', '.vcf'), 0o777)
        output_ensembl = "/opt/vep/.vep/" + sample_id + '/' + vep_file.replace('.tsv', '.vcf')

        # run docker
        client.containers.run("ensemblorg/ensembl-vep:release_104.3",
                              "./vep --cache --cache_version 104 --offline --format ensembl --vcf --force_overwrite -dir_cache /opt/vep/.vep --input_file {} --output_file {} --assembly GRCh38".format(
                                  input_ensembl, output_ensembl),
                              stdout=True,
                              stderr=False,
                              remove=True,
                              volumes=[f'{vep_dir}:/opt/vep/.vep'])
