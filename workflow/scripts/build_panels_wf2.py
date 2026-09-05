import pandas as pd
import argparse



def get_study_panels(study, phy_out, pvi_input, drug_prior, drug_filter, out_file):

    """
    Build Dataframe for sample panels of the study.

    Args:
        study (str): Name of the study of interest
        phy_out (str): Path to the Phyclone of the study
        pvi_input (str): Path to the Pyclone-VI of the study
        drug_prior(str): Path to drug prioritization file
        drug_filter(str): Drug Analysis Filter: clinical | discovery
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
        phy_df = pd.read_table(phy_out, sep="\t")
    except Exception as e:
        raise ValueError(f"Error reading PyClone-VI tree file: {e}")


    phy_df = phy_df[(phy_df['clone_id'] != -1) & (phy_df['clonal_prev'] != 0)].reset_index(drop=True)
    
    # sample statistics 
    panel = (
        phy_df.groupby("sample_id", as_index=False)
        .agg(
            num_mutations=("mutation_id", "nunique"),
            num_clones=("clone_id", "nunique")
        )
        .reset_index()
    )

    # Gathering sex information
    phy_df["chrom"] = (
        phy_df["mutation_id"]
        .astype(str)
        .str.split(":", n=1)
        .str[0]
    )

    # Gathering sex information
    phy_df["chrom"] = (
    phy_df["mutation_id"]
    .astype(str)
    .str.split(":", n=1)
    .str[0]
    )
    phy_df["chrom"] = phy_df["chrom"].str.replace("^chr", "", regex=True).str.upper()



    sex_by_sample = (
        phy_df
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


    # GET DRUG PRIORITIZATION: Total drugs and BTC
    try:
        drugs = pd.read_table(drug_prior)
    except Exception as e:
        raise ValueError(f"Error reading drug prioritization of study {study}: {e}")
    
    drugs = drugs[drugs["Clone"] != -1]
    drugs["dScore"] = pd.to_numeric(drugs["dScore"], errors="coerce")

    # Sort Status and Intereaction type column
    status_order = ["APPROVED", "CLINICAL_TRIALS", "EXPERIMENTAL"]
    drugs["Status"] = pd.Categorical(drugs["Status"], categories=status_order, ordered=True)
    
    interaction_type_order = ["DIRECT_TARGET", "BIOMARKER", "PATHWAY_MEMBER"]
    drugs["Interaction Type"] = pd.Categorical(drugs["Interaction Type"], categories=interaction_type_order, ordered=True)

     # Apply Drug filter if clinical mode was set up
    if drug_filter == 'clinical':
        drugs = drugs[(drugs['Status'] != 'EXPERIMENTAL') & (drugs['Interaction Type'] != 'PATHWAY_MEMBER')]

    # Disaggregate samples by VAF values  
    df_expanded = (
    drugs.assign(VAF=drugs["VAF"].str.split("; "))
      .explode("VAF")
      .assign(
          sample_id=lambda d: d["VAF"].str.split(": ").str[0],
          VAF=lambda d: d["VAF"].str.split(": ").str[1].astype(float))
     )

    # Filter out drug occurrences with zero VAF
    drugs_hits = df_expanded[df_expanded['VAF'] != 0].copy()

    # Count Total drugs and BTCs
    drugs_hits["is_BTC"] = (drugs_hits["dScore"] > 0.7) & (drugs_hits["gScore"] > 0.6)
    drug_stats = (
        drugs_hits.groupby("sample_id", as_index=False)
        .agg(
            total_drugs=("Drug", "nunique"),
            BTCs=("is_BTC", "sum")
        )
    )

    # Merge sex and tumour content
    panel = (panel
             .merge(sex_by_sample, on="sample_id", how="left")
             .merge(tumour_panel, on="sample_id", how="left")
             .merge(drug_stats, on='sample_id', how='inner'))

    # Add study column and reorder columns
    panel["study"] = study
    panel = panel[["study", "sample_id", "sex",
                   "tumour_content", "num_mutations", "num_clones", "total_drugs", "BTCs"]]

    panel.to_csv(out_file, sep="\t", index=False)

    return panel



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--study", action='store', required=True)
    input_parser.add_argument("--phy_out", action='store', required=True)
    input_parser.add_argument("--pvi_input", action='store', required=True)
    input_parser.add_argument("--drug_prioritization", action='store', required=True)
    input_parser.add_argument("--drug_filter", action='store', required=True)
    input_parser.add_argument("--out_file", action='store', required=True)

    args = input_parser.parse_args()


    get_study_panels(args.study, args.phy_out, args.pvi_input, args.drug_prioritization, args.drug_filter, args.out_file)
