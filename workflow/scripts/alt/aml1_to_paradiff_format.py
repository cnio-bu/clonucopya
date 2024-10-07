import pandas as pd
import numpy as np
import argparse


if __name__ == '__main__':

    dict_mutations_normal = {}
    dict_mutations_tumor = {}
    dict_mutations_relapse = {}

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_dir", action='store')
    input_parser.add_argument("--mutations", nargs='+', action='append')
    input_parser.add_argument("--output_dir", action='store')
    args = input_parser.parse_args()
    input_dir = args.input_dir
    input_files = args.mutations

    # array of input mutations files
    arr_mutations_files = np.array(input_files)

    # iterate by mutations files
    for idx_mut_file, row_mut_file in np.ndenumerate(arr_mutations_files):

        # reset variables
        l_variable_id = []

        l_nrm_933124_a = []
        l_nrm_933124_c = []
        l_nrm_933124_g = []
        l_nrm_933124_t = []

        l_tum_933124_a = []
        l_tum_933124_c = []
        l_tum_933124_g = []
        l_tum_933124_t = []

        l_rel_933124_a = []
        l_rel_933124_c = []
        l_rel_933124_g = []
        l_rel_933124_t = []

        dict_normal = {}
        dict_tumor = {}
        dict_relapse = {}

        l_case = []
        l_chr = []
        l_loc = []
        l_mut = []
        dict_mutations = {}

        mut_file_path = row_mut_file

        # file containing all mutations
        mutations_tsv = pd.read_csv(mut_file_path, sep=',')

        # iterate by mutations
        for idx_mutation_tsv, row_mutation_tsv in mutations_tsv.iterrows():

            chromosome = row_mutation_tsv['Chr']
            start_GRCh37 = row_mutation_tsv['Start (GRCh37)']
            ref_allele = row_mutation_tsv['Reference allele']
            var_allele = row_mutation_tsv['Variant allele']

            # normal roche 454
            roche_454_nrm_reads_ref = row_mutation_tsv['454.nrm.reads.ref']
            roche_454_nrm_reads_var = row_mutation_tsv['454.nrm.reads.var']
            roche_454_nrm_var_freq = row_mutation_tsv['454.nrm.var.freq']

            # tumor roche 454
            roche_454_tum_reads_ref = row_mutation_tsv['454.tum.reads.ref']
            roche_454_tum_reads_var = row_mutation_tsv['454.tum.reads.var']
            roche_454_tum_var_freq = row_mutation_tsv['454.tum.var.freq']

            # relapse roche 454
            roche_454_rel_reads_ref = row_mutation_tsv['454.rel.reads.ref']
            roche_454_rel_reads_var = row_mutation_tsv['454.rel.reads.var']
            roche_454_rel_var_freq = row_mutation_tsv['454.rel.var.freq']

            # variable_id l_933124_a l_933124_c l_933124_g l_933124_t
            variable_id = 'chr{}_{}_{}'.format(chromosome, start_GRCh37, ref_allele)
            list.append(l_variable_id, variable_id)

            # Reference
            if str.upper(ref_allele) == 'A':
                list.append(l_nrm_933124_a, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_a, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_a, int(roche_454_rel_reads_ref))
            else:
                list.append(l_nrm_933124_a, 0)
                list.append(l_tum_933124_a, 0)
                list.append(l_rel_933124_a, 0)

            if str.upper(ref_allele) == 'C':
                list.append(l_nrm_933124_c, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_c, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_c, int(roche_454_rel_reads_ref))
            else:
                list.append(l_nrm_933124_c, 0)
                list.append(l_tum_933124_c, 0)
                list.append(l_rel_933124_c, 0)

            if str.upper(ref_allele) == 'G':
                list.append(l_nrm_933124_g, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_g, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_g, int(roche_454_rel_reads_ref))
            else:
                list.append(l_nrm_933124_g, 0)
                list.append(l_tum_933124_g, 0)
                list.append(l_rel_933124_g, 0)

            if str.upper(ref_allele) == 'T':
                list.append(l_nrm_933124_t, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_t, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_t, int(roche_454_rel_reads_ref))
            else:
                list.append(l_nrm_933124_t, 0)
                list.append(l_tum_933124_t, 0)
                list.append(l_rel_933124_t, 0)

            # Variant
            if str.upper(var_allele) == 'A':
                list.append(l_nrm_933124_a, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_a, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_a, int(roche_454_rel_reads_ref))
            elif not str.upper(ref_allele) == 'A':
                list.append(l_nrm_933124_a, 0)
                list.append(l_tum_933124_a, 0)
                list.append(l_rel_933124_a, 0)

            if str.upper(var_allele) == 'C':
                list.append(l_nrm_933124_c, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_c, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_c, int(roche_454_rel_reads_ref))
            else:
                list.append(l_nrm_933124_c, 0)
                list.append(l_tum_933124_c, 0)
                list.append(l_rel_933124_c, 0)

            if str.upper(var_allele) == 'G':
                list.append(l_nrm_933124_g, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_g, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_g, int(roche_454_rel_reads_ref))
            else:
                list.append(l_nrm_933124_g, 0)
                list.append(l_tum_933124_g, 0)
                list.append(l_rel_933124_g, 0)

            if str.upper(var_allele) == 'T':
                list.append(l_nrm_933124_t, int(roche_454_nrm_reads_ref))
                list.append(l_tum_933124_t, int(roche_454_tum_reads_ref))
                list.append(l_rel_933124_t, int(roche_454_rel_reads_ref))
            else:
                list.append(l_nrm_933124_t, 0)
                list.append(l_tum_933124_t, 0)
                list.append(l_rel_933124_t, 0)

            # dictionary of lists
            dict_normal = {'l_variable_id': l_variable_id,
                           'l_933124_a': l_nrm_933124_a,
                           'l_933124_c': l_nrm_933124_c,
                           'l_933124_g': l_nrm_933124_g,
                           'l_933124_t': l_nrm_933124_t}

            dict_tumor = {'l_variable_id': l_variable_id,
                          'l_933124_a': l_tum_933124_a,
                          'l_933124_c': l_tum_933124_c,
                          'l_933124_g': l_tum_933124_g,
                          'l_933124_t': l_tum_933124_t}

            dict_relapse = {'l_variable_id': l_variable_id,
                            'l_933124_a': l_rel_933124_a,
                            'l_933124_c': l_rel_933124_c,
                            'l_933124_g': l_rel_933124_g,
                            'l_933124_t': l_rel_933124_t}

            # saving the mutations
            if str.__contains__(mut_file_path, 'normal'):
                list.append(l_case, 'normal')
            elif str.__contains__(mut_file_path, 'tumor'):
                list.append(l_case, 'tumor')
            elif str.__contains__(mut_file_path, 'relapse'):
                list.append(l_case, 'relapse')

            list.append(l_chr, chromosome)
            list.append(l_loc, start_GRCh37)
            list.append(l_mut, str.upper(ref_allele) + '/' + str.upper(var_allele))

        # saving the dataframe
        dict_mutations = {'case': l_case,
                          'chr': l_chr,
                          'loc': l_loc,
                          'mut': l_mut}
        df_mutations = pd.DataFrame(dict_mutations)
        df_mutations.to_csv('{}/{}'.format(args.output_dir, 'mutations.tsv'), sep='\t', index=False, mode='a')

        if str.__contains__(mut_file_path, 'normal'):
            df_output = pd.DataFrame(dict_normal)
            df_output.to_csv('{}/{}'.format(args.output_dir, '933124_normal.pars.txt'), sep='\t',
                             header=["Variable_ID", "933124_A", "933124_C", "933124_G", "933124_T"], index=False,
                             mode='a')

        elif str.__contains__(mut_file_path, 'tumor'):
            df_output = pd.DataFrame(dict_tumor)
            df_output.to_csv('{}/{}'.format(args.output_dir, '933124_tumor.pars.txt'), sep='\t',
                             header=["Variable_ID", "933124_A", "933124_C", "933124_G", "933124_T"], index=False,
                             mode='a')

        elif str.__contains__(mut_file_path, 'relapse'):
            df_output = pd.DataFrame(dict_relapse)
            df_output.to_csv('{}/{}'.format(args.output_dir, '933124_relapse.pars.txt'), sep='\t',
                             header=["Variable_ID", "933124_A", "933124_C", "933124_G", "933124_T"], index=False,
                             mode='a')
