import pandas as pd
import pybedtools
import argparse

def create_pyclone_vi_input(sample_id, mutations_file, cnv_file, output_file):
    """
    Intersect mutations with copy number regions. Copy number must contain mutation to be valid.

    Args:
        sample_id (str): Sample ID of the pair mutation-cnv files
        mutations_file (str): Path to mutation file of the processed sample
        cnv_file (str): Path to copy number file of the processed sample
        output_file (str, optional): Path to output TSV file.

    Returns:
        pandas.DataFrame: Intersected sample's mutations-cnvs pair
    """

    # Load mutations and copy number files
    mutations = pd.read_csv(mutations_file, sep='\t')
    
    cnv = pd.read_csv(cnv_file, sep='\t')
    
    # Create bedtool objects
    mut_bed = pybedtools.BedTool.from_dataframe(
        mutations[['CHROM', 'POS', 'POS', 'mutation_id', 'REF', 'ALT']]
        .rename(columns={'POS': 'start', 'CHROM': 'chrom'})
    )
    
    cnv_bed = pybedtools.BedTool.from_dataframe(
        cnv[['Chrom', 'Start', 'End', 'major_cn', 'minor_cn', 'normal_cn']]
        .rename(columns={'Chrom': 'chrom'})
    )
    
    # Intersect mutations with copy number regions
    intersect = mut_bed.intersect(cnv_bed, wa=True, wb=True)

    # Format intersected object into a dictionary
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
    
    # Convert dictionary to dataframe and save it
    output = pd.DataFrame(result)
    output.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--sample_id", action='store', required=True)
    input_parser.add_argument("--mutations", action='store', required=True)
    input_parser.add_argument("--cnvs", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args()


    create_pyclone_vi_input(args.sample_id,args.mutations, args.cnvs, args.output_file)
