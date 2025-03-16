import pandas as pd
import pybedtools
import argparse

def create_pyclone_vi_input(sample_id, mutations_file, cnv_file, output_file):
    # Load mutations
    mutations = pd.read_csv(mutations_file, sep='\t')
    
    # Load CNV data
    cnv = pd.read_csv(cnv_file, sep='\t')
    
    # Create BedTool objects with required columns
    mut_bed = pybedtools.BedTool.from_dataframe(
        mutations[['CHROM', 'POS', 'POS', 'mutation_id', 'REF', 'ALT']]
        .rename(columns={'POS': 'start', 'CHROM': 'chrom'})
    )
    
    cnv_bed = pybedtools.BedTool.from_dataframe(
        cnv[['Chrom', 'Start', 'End', 'major_cn', 'minor_cn', 'normal_cn']]
        .rename(columns={'Chrom': 'chrom'})
    )
    
    # Intersect mutations with CNV regions
    intersect = mut_bed.intersect(cnv_bed, wa=True, wb=True)

    # Process intersected data preserving mutation_id
    result = []
    for item in intersect:
        mutation_id = item[3]
        result.append({
            'mutation_id': mutation_id,
            'sample_id': sample_id,
            'ref_counts': mutations.loc[mutations['mutation_id'] == mutation_id, 'ref_counts'].iloc[0],
            'alt_counts': mutations.loc[mutations['mutation_id'] == mutation_id, 'alt_counts'].iloc[0],
            'major_cn': int(item[9]),
            'minor_cn': int(item[10]),
            'normal_cn': int(item[11])
        })
    
    # Create DataFrame and save to file
    output = pd.DataFrame(result)
    output.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--sample_id", action='store', required=True)
    input_parser.add_argument("--mutations", action='store', required=True)
    input_parser.add_argument("--cnvs", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args()

    # Process the data
    create_pyclone_vi_input(args.sample_id,args.mutations, args.cnvs, args.output_file)
