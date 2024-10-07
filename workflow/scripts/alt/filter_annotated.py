
import io
import pandas as pd


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
    info_csq = str.split(df_row["INFO"], sep=";")
    #print(info_csq)
    csq_fields = str.split(info_csq[12], sep="|")
    csq_fields[33]
    print(csq_fields[33])
    #print(f'AF={csq_fields[33]} gnomAD_AF={csq_fields[42]}')


if __name__ == '__main__':

    print('start')

    df_input = read_vcf('/home/cleon/workspace_pandrugs/aux/20ID00738_mutect.vep.parted.vcf')

    for index, row in df_input.iterrows():
        print_row(row)

    print('end')
