# Note: Reshape your data either using array.reshape(-1, 1) if your data has a single feature or array.reshape(1, -1)
#       if it contains a single sample.

import pandas as pd

input_from_pyclone = "/home/cleon/extdata/tumor_TCGA-05-4244-01A-01D-1105-08_pyclone.tsv" #"/everest/pyclone_data/output/tables/loci.tsv"
output_to_citup = "/home/cleon/extdata/tumor_TCGA-05-4244-01A-01D-1105-08_cituped.tsv" #"/everest/pyclone_data/output/tables/freq.txt"
freq_assignment_dict = {}
l_freqs = []

if __name__ == '__main__':

    loci = pd.read_csv(input_from_pyclone, sep="\t", low_memory=False)
    for index, row in loci.iterrows():
        sample_id = row['sample_id']
        t_alt_count = row['alt_counts']
        if t_alt_count == 0:
            break
        t_ref_count = row['ref_counts']
        vaf = t_alt_count / (t_ref_count + t_alt_count)
        l_freqs.append(vaf)

    freq_assignment_dict[sample_id] = l_freqs
    print(freq_assignment_dict)
    dict_df = pd.DataFrame({key: pd.Series(value) for key, value in freq_assignment_dict.items()})
    dict_df.to_csv(output_to_citup, index=False, header=False, sep="\t")
    # print(dict_df)
