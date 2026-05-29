import pandas as pd
import argparse



def get_study_panels(phy_out, out_file):

    """
    Build Dataframe for sample panels of the study.

    Args:
        phy_out (str): Path to the Phyclone of the study
        out_file (str): Path to the output file (TSV)

    Return:
        DataFrame of the main statistic of the samples of the same study
        
    """

    # Load Phyclone tree clone dataframe
    try:
       tree_df = pd.read_table(phy_out, sep='\t')
    except Exception as e:
        raise ValueError(f"Error reading Phyclone output file file: {e}")

    study_filt = tree_df[['mutation_id', 'clone_id', 'sample_id']]


    panel = (
    study_filt
    .groupby('sample_id')
    .agg(
        num_mutations=('mutation_id', 'nunique'),
        num_clones=('clone_id', 'nunique')
    )
    .reset_index()
    )

    
    panel.to_csv(out_file, sep='\t', index=False)
    
    return panel



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--phy_out", action='store', required=True)
    input_parser.add_argument("--out_file", action='store', required=True)

    args = input_parser.parse_args()


    get_study_panels(args.phy_out, args.out_file)
