import subprocess
import argparse
import numpy as np
import pandas as pd
import os

if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--in_files", nargs='+', action='append')
    input_parser.add_argument("--in_dir", action='store')
    input_parser.add_argument("--out_dir", action='store')
    input_parser.add_argument("--samples", action='store')
    args = input_parser.parse_args()
    input_dir = args.in_dir
    input_files = args.in_files
    config_dataset = pd.read_csv(args.samples, sep='\t')

    # array of unique patients
    arr_patients = np.array(config_dataset['patient'].unique())

    # iterate by patient
    patient_id = ''
    all_sample_ids = []
    all_purities = []
    for idx_patient, row_patient in np.ndenumerate(arr_patients):
        # filter by Patient
        df_patient = pd.DataFrame(config_dataset.loc[config_dataset['patient'] == row_patient])

        for idx_patient_sample_tsv, row_patient_sample_tsv in df_patient.iterrows():
            patient_id = row_patient_sample_tsv['patient']
            case = row_patient_sample_tsv['case']
            purity = row_patient_sample_tsv['Purity']
            all_sample_ids.append(case)
            all_purities.append(purity)

    all_samples_tsv = []
    sample_tsv = ''
    '''arr = np.array(input_files)
    for idx, row_sample_tsv in np.ndenumerate(arr):

        # reset variables
        output_dir = args.out_dir

        path_sample_tsv_splited = str(row_sample_tsv).split('/')
        name_file_sample_tsv = path_sample_tsv_splited[len(path_sample_tsv_splited)-1]
        sample_id = name_file_sample_tsv.replace('.tsv', '')

        sample_tsv = input_dir + '/' + str(row_sample_tsv)
        all_samples_tsv.append(sample_tsv)
        output_dir = output_dir + '/{}'.format(sample_id)

        process = subprocess.Popen(["PyClone", "run_analysis_pipeline", "--in_files", sample_tsv,
                                    "--working_dir", output_dir,
                                    "--max_clusters", "10", "--num_iters", "1000"])
        process.communicate()'''

    # do it all together now once again by patient
    output_dir = args.out_dir
    output_dir = output_dir + "/patient_{}".format(patient_id)

    # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pyclone/patient_47
    # if not os.path.exists(output_dir + '/patient_' + str(patient_id)):
    #     os.makedirs(output_dir + '/patient' + str(patient_id))

    command = ["PyClone", "run_analysis_pipeline", "--in_files", "/home/cleon/PycharmProjects/MSc-Project/out/pyclone/UPN933124/input_UPN933124.tsv"]
    for idx_com, row_com_sample_tsv in np.ndenumerate(all_samples_tsv):
        command.append(str(row_com_sample_tsv))

    '''command.append("--tumour_contents")
    for idx_pur, row_pur in np.ndenumerate(all_purities):
        command.append(str(row_pur))'''

    command.append("--working_dir")
    command.append(str(output_dir))
    command.append("--max_clusters")
    # it seems CITUP num-nodes=8
    command.append("8")
    command.append("--num_iters")
    command.append("1000")

    process = subprocess.Popen(command)
    process.communicate()
