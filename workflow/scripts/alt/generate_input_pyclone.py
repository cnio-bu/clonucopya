# This is a sample Python script.

# Press Mayús+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import io
import pandas as pd

l_mutation_id = []
l_ref_counts = []
l_var_counts = []
l_normal_cn = []
l_minor_cn = []
l_major_cn = []


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def print_row(df_row):
    if len(df_row["REF"]) == 1 == len(df_row["ALT"]):
        sample_format = str.split(df_row["TCGA-05-4244-01A-01D-1105-08"], sep=":")
        sample_ad = str.split(sample_format[1], sep=",")
        """print(df_row["CHROM"] + " " + df_row["REF"] + " " + str(df_row["ALT"]) +
              " " + str(df_row["TCGA-05-4244-01A-01D-1105-08"]) +
              " " + str(sample_ad[0]) +
              " " + str(sample_ad[1]))"""

        # list of mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn
        list.append(l_mutation_id, "TCGA-05-4244-01A-01D-1105-08" + ":" + str(df_row["CHROM"]) + ":" + str(df_row["POS"]))
        list.append(l_ref_counts, str(sample_ad[0]))
        list.append(l_var_counts, str(sample_ad[1]))
        list.append(l_normal_cn, "2")
        list.append(l_minor_cn, "0")
        list.append(l_major_cn, "2")


if __name__ == '__main__':

    df_input = read_vcf('/everest/vep_data/output/TCGA-05-4244-01A-01D-1105-08.vcf')

    for index, row in df_input.iterrows():
        if len(row.values[3]) == len(row.values[4]) == 1:
            print_row(row)

    # dictionary of lists
    # mutation_id ref_counts var_counts normal_cn minor_cn major_cn
    dict2pyclone = {'mutation_id': l_mutation_id,
                    'ref_counts': l_ref_counts,
                    'var_counts': l_var_counts,
                    'normal_cn': l_normal_cn,
                    'minor_cn': l_minor_cn,
                    'major_cn': l_major_cn}

    df_output = pd.DataFrame(dict2pyclone)

    # saving the dataframe
    df_output.to_csv('/everest/vep_data/output/TCGA-05-4244-01A-01D-1105-08_pyclone.tsv',
                     sep='\t',
                     index=False)

    # preparing data for citup
    '''f = open("/home/cleon/workspace_pandrugs/varca/out/input_citup.txt", "w+")
    f.write("Num_mutations: %d\r\n" % 5)
    f.write("Num_samples: %d\r\n" % 1)
    f.write("Error_rate: %d\r\n" % 0.03)

    dict2citup = {}

    f.close()'''
