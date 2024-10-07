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
    input_parser.add_argument("--samples", nargs='+', action='append')
    input_parser.add_argument("--patient", action='store')
    input_parser.add_argument("--results_mutations", action='store')
    input_parser.add_argument("--vep_dir", action='store')
    args = input_parser.parse_args()

    # i.e: /home/cleon/PycharmProjects/MSc-Project/out/mapscape/933124/results_mutations.tsv file
    results_mutations_tsv = pd.read_csv(args.results_mutations, sep=',')

    # i.e: '933124'
    patient_id = args.patient

    # i.e: '/home/cleon/vep_data'
    vep_dir = args.vep_dir

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

            os.chmod(vep_dir + '/' + str(patient_id), 0o777)

            os.chmod(vep_dir + '/' + str(patient_id) + '/' + str(sample_id), 0o777)

            # i.e: ensembl_tumor.0.tsv
            vep_file = 'ensembl_{}.{}.tsv'.format(sample_id, cluster_id)
            # print('>>>>>>>>>>>>>>{}'.format(vep_file))

            # i.e: /opt/vep/.vep/933124/tumor/ensembl_tumor.0.tsv
            input_ensembl = "/opt/vep/.vep/" + str(patient_id) + "/" + str(sample_id) + "/" + vep_file

            # i.e: /opt/vep/.vep/933124/tumor/ensembl_tumor.0.vcf
            Path(vep_dir + '/' + str(patient_id) + '/' + str(sample_id) + '/' + vep_file.replace('.tsv', '.vcf')).write_text(
                "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO")
            os.chmod(vep_dir + '/' + str(patient_id) + '/' + str(sample_id) + '/' + vep_file.replace('.tsv', '.vcf'), 0o777)
            output_ensembl = "/opt/vep/.vep/" + str(patient_id) + '/' + sample_id + '/' + vep_file.replace('.tsv', '.vcf')

            # run docker
            client.containers.run("ensemblorg/ensembl-vep:release_104.3",
                                  "./vep --cache --cache_version 104 --offline --format ensembl --vcf --force_overwrite -dir_cache /opt/vep/.vep/{} --input_file {} --output_file {} --assembly GRCh38".format(
                                      patient_id, input_ensembl, output_ensembl),
                                  stdout=True,
                                  stderr=False,
                                  remove=True,
                                  volumes=[f'{vep_dir}/{patient_id}:/opt/vep/.vep/{patient_id}'])
