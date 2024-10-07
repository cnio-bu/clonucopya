import re

import vcf
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


n = 0
info_values = []
csq_groups = []
r = re.compile('^CSQ')
if __name__ == '__main__':

    print('start')
    df_input = read_vcf('/home/cleon/workspace_pandrugs/aux/current/all.vep.vcf')
    for index, row in df_input.iterrows():
        #print(row['INFO'])
        info_values = str.split(row['INFO'], sep=';')
        for iv in info_values:
            if re.search('^CSQ', iv):
                #print(iv)
                csq_groups = str.split(iv, sep=',')
                for group in csq_groups:
                    csq_group_values = str.split(group, sep='|')
                    if csq_group_values[2] == 'HIGH' or csq_group_values[2] == 'MODERATE':
                        n = n + 1
                        print(csq_group_values[2])
                        break
    print(n)
    print('end')

    """i = 0
    
    vcf_reader2write = vcf.Reader(open('/home/cleon/workspace_pandrugs/aux/current/all.vep.reduced.vcf'))
    
    
    vcf_reader = vcf.Reader(open('/home/cleon/workspace_pandrugs/aux/current/all.vep.vcf', 'r'))
    vcf_writer = vcf.Writer(open('/dev/null', 'w'), vcf_reader2write)
    for record in vcf_reader:
        if record.is_snp:
            csq = record.INFO['CSQ']
            for csq_record in csq:
                #print(csq_record)
                csq_record_fields = str.split(csq_record, sep='|')
                if (csq_record_fields[37] != '' and float(csq_record_fields[37]) < 0.01) or \
                        (csq_record_fields[45] != '' and float(csq_record_fields[45]) < 0.01):
                    #print(f'Allele={csq_record_fields[0]} AF={csq_record_fields[37]} gnomAD_AF={csq_record_fields[45]}')
                    if csq_record_fields[2] == 'HIGH' or csq_record_fields[2] == 'MODERATE':
                        i = i + 1
                        print(f'IMPACT={csq_record_fields[2]}')
                        vcf_writer.write_record(record)
                        record.add_filter()
    vcf_writer.close()
    print(i)"""
