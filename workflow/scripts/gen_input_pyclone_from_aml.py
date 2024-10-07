import pandas as pd
import numpy as np
import argparse


def index_in_array(arr, index, loc):
    try:
        a = arr[index]
        if len(a) > 0:
            return True
        else:
            return False
    except KeyError:
        return False


if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_dir", action='store')
    input_parser.add_argument("--samples", action='store')
    input_parser.add_argument("--mutations", action='store')
    input_parser.add_argument("--output_dir", action='store')
    args = input_parser.parse_args()
    input_dir = args.input_dir
    samples_tsv = pd.read_csv(args.samples, sep='\t')

    # file containing all mutations
    mutations_tsv = pd.read_csv(args.mutations, sep='\t')

    # array of unique patients
    arr_patients = np.array(samples_tsv['patient'].unique())

    # iterate by patient
    for idx_patient, row_patient in np.ndenumerate(arr_patients):
        # filter by Patient
        df_patient = pd.DataFrame(samples_tsv.loc[samples_tsv['patient'] == row_patient])

        # array of unique cases by patient
        arr_cases_by_patient = np.array(df_patient['case'].unique())

        for idx_patient_sample_tsv, row_patient_sample_tsv in df_patient.iterrows():

            # reset all lists
            l_mutation_id = []
            l_ref_counts = []
            l_var_counts = []
            l_normal_cn = []
            l_minor_cn = []
            l_major_cn = []
            l_variant_case = []

            dict2pyclone = {}

            case = str(row_patient_sample_tsv['case'])
            patient = str(row_patient_sample_tsv['patient'])
            gender = str(row_patient_sample_tsv['gender'])

            reads_by_base_tsv = pd.read_csv(input_dir + '/' + row_patient_sample_tsv['mutations'], sep='\t')

            # filter by all cases of the patient
            # the point is to get the same mutations also in the next round of the loop, when the other mutations file
            # will be used for all cases. It is a requirement of PyClone
            # df_case = pd.DataFrame(mutations_tsv.loc[mutations_tsv['case'] in arr_cases_by_patient])
            df_case = pd.DataFrame(mutations_tsv.loc[mutations_tsv['case'].apply(lambda x: x in arr_cases_by_patient)])

            # iterate by each SNP
            for idx_case, row_case in df_case.iterrows():
                ref_counts = 0
                var_counts = 0

                chr = row_case['chr']
                pos = row_case['loc']

                mutation = row_case['mut']
                ref_allele = mutation.split('/')[0]
                var_allele = mutation.split('/')[1]

                variable_id = 'chr' + str(chr) + '_' + str(pos) + '_' + str(ref_allele)

                df_mut = reads_by_base_tsv.loc[reads_by_base_tsv['Variable_ID'] == variable_id]

                if df_mut.empty:
                    variable_id = 'chr' + str(chr) + '_' + str(pos) + '_' + str(ref_allele).lower()
                    df_mut = reads_by_base_tsv.loc[reads_by_base_tsv['Variable_ID'] == variable_id]

                case_mut_a = 0
                if index_in_array(df_mut, patient + '_A', pos) and len(df_mut[patient + '_A']) > 0:
                    case_mut_a = df_mut[patient + '_A'].iat[0]
                elif index_in_array(df_mut, patient + '_a', pos) and len(df_mut[patient + '_a'] > 0):
                    case_mut_a = df_mut[patient + '_a'].iat[0]

                case_mut_c = 0
                if index_in_array(df_mut, patient + '_C', pos) and len(df_mut[patient + '_C']) > 0:
                    case_mut_c = df_mut[patient + '_C'].iat[0]
                elif index_in_array(df_mut, patient + '_c', pos) and len(df_mut[patient + '_c'] > 0):
                    case_mut_c = df_mut[patient + '_c'].iat[0]

                case_mut_g = 0
                if index_in_array(df_mut, patient + '_G', pos) and len(df_mut[patient + '_G']) > 0:
                    case_mut_g = df_mut[patient + '_G'].iat[0]
                elif index_in_array(df_mut, patient + '_g', pos) and len(df_mut[patient + '_g'] > 0):
                    case_mut_g = df_mut[patient + '_g'].iat[0]

                case_mut_t = 0
                if index_in_array(df_mut, patient + '_T', pos):
                    case_mut_t = df_mut[patient + '_T'].iat[0]
                if index_in_array(df_mut, patient + '_t', pos):
                    case_mut_t = df_mut[patient + '_t'].iat[0]

                if ref_allele == 'A' or ref_allele == 'a':
                    ref_counts = case_mut_a
                elif ref_allele == 'C' or ref_allele == 'c':
                    ref_counts = case_mut_c
                elif ref_allele == 'G' or ref_allele == 'g':
                    ref_counts = case_mut_g
                elif ref_allele == 'T' or ref_allele == 't':
                    ref_counts = case_mut_t
                ref_counts = ref_counts if not pd.isna(ref_counts) else 0

                if var_allele == 'A' or var_allele == 'a':
                    var_counts = case_mut_a
                elif var_allele == 'C' or var_allele == 'c':
                    var_counts = case_mut_c
                elif var_allele == 'G' or var_allele == 'g':
                    var_counts = case_mut_g
                elif var_allele == 'T' or var_allele == 't':
                    var_counts = case_mut_t
                var_counts = var_counts if not pd.isna(var_counts) else 0

                # calculate minor_cn
                minor_cn = 0

                # TODO we assume diploid(2N)
                # calculate normal_cn, minor_cn, major_cn
                if str.lower(gender) == 'male' and (chr == 'X' or chr == 'x' or chr == 'Y' or chr == 'y'):
                    normal_cn = 1
                    major_cn = 1
                    minor_cn = 0
                else:
                    normal_cn = 2
                    major_cn = 2
                    minor_cn = 0

                # mutation_id ref_counts var_counts normal_cn minor_cn major_cn variant_case
                list.append(l_mutation_id, 'chr' + str(chr) + ':' + str(pos))
                list.append(l_ref_counts, int(ref_counts))
                list.append(l_var_counts, int(var_counts))
                list.append(l_normal_cn, normal_cn)
                list.append(l_minor_cn, minor_cn)
                list.append(l_major_cn, major_cn)
                list.append(l_variant_case, case)

            # dictionary of lists
            # mutation_id ref_counts var_counts normal_cn minor_cn major_cn variant_case
            dict2pyclone = {'mutation_id': l_mutation_id,
                            'ref_counts': l_ref_counts,
                            'var_counts': l_var_counts,
                            'normal_cn': l_normal_cn,
                            'minor_cn': l_minor_cn,
                            'major_cn': l_major_cn,
                            'variant_case': l_variant_case}

            df_output = pd.DataFrame(dict2pyclone)

            # saving the dataframe
            df_output.to_csv("{}/{}.tsv".format(args.output_dir, case), sep='\t', index=False, mode='a')
