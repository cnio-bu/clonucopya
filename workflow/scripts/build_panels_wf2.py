import pandas as pd
import argparse



def get_study_panels(study, phy_out, pvi_input, out_file):

    """
    Build Dataframe for sample panels of the study.

    Args:
        study (str): Name of the study of interest
        phy_out (str): Path to the Phyclone of the study
        pvi_input (str): Path to the Pyclone-VI of the study
        out_file (str): Path to the output file (TSV)

    Return:
        DataFrame of the main statistic of the samples of the same study
        
    """

   # Load PyClone-VI input to obtain tumour_content information
    try:
        cluster_df = pd.read_table(pvi_input, sep="\t")
    except Exception as e:
        raise ValueError(f"Error reading PyClone-VI input file: {e}")

    tumour_panel = (
        cluster_df
        .groupby("sample_id", as_index=False)
        .agg(tumour_content=("tumour_content", "first"))
    )

    # Load Phyclone tree/clone dataframe
    try:
        tree_df = pd.read_table(phy_out, sep="\t")
    except Exception as e:
        raise ValueError(f"Error reading PyClone-VI tree file: {e}")

    # sample statistics
    study_filt = tree_df[["mutation_id", "clone_id", "sample_id"]].copy()
    study_filt["clone_id"] = study_filt["clone_id"].astype(int).astype(str)

    # Gathering sex information
    tree_df["chrom"] = (
    tree_df["mutation_id"]
    .astype(str)
    .str.split(":", n=1)
    .str[0]
    )
    tree_df["chrom"] = tree_df["chrom"].str.replace("^chr", "", regex=True).str.upper()



    sex_by_sample = (
        tree_df
        .groupby("sample_id")["chrom"]
        .agg(
            has_X=lambda s: (s == "X").any(),
            has_Y=lambda s: (s == "Y").any(),
        )
        .reset_index()
    )
    
    sex_by_sample["sex"] = "Unknown"
    sex_by_sample.loc[sex_by_sample["has_X"], "sex"] = "female"
    sex_by_sample.loc[sex_by_sample["has_Y"], "sex"] = "male"

    sex_by_sample = sex_by_sample[["sample_id", "sex"]]

    # Aggregate number of mutations and clones per sample
    panel = (
        study_filt
        .groupby("sample_id")
        .agg(
            num_mutations=("mutation_id", "nunique"),
            num_clones=("clone_id", "nunique"),
        )
        .reset_index()
    )

    # Merge sex and tumour content
    panel = panel.merge(sex_by_sample, on="sample_id", how="left")
    panel = panel.merge(tumour_panel, on="sample_id", how="left")

    # Add study column and reorder columns
    panel["study"] = study
    panel = panel[["study", "sample_id", "sex",
                   "tumour_content", "num_mutations", "num_clones"]]

    panel.to_csv(out_file, sep="\t", index=False)

    return panel



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--study", action='store', required=True)
    input_parser.add_argument("--phy_out", action='store', required=True)
    input_parser.add_argument("--pvi_input", action='store', required=True)
    input_parser.add_argument("--out_file", action='store', required=True)

    args = input_parser.parse_args()


    get_study_panels(args.study, args.phy_out, args.pvi_input,args.out_file)
