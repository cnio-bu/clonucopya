import pandas as pd
import glob
import os
import argparse



def get_study_panels(study, samplesheet, mut_dir, intersect_combined, drug_prior, drug_filter, out_file):

    """
    Build Dataframe for sample panels of the study.

    Args:
        study (str): name of the study
        samplesheet (str): Path to the samplesheet (CSV)
        mut_dir (str): Path to the files from mutation_prep study
        intersect_combined(str): Path to the intersect combined of the study
        drug_prior(str): Path to drug prioritization file
        drug_filter(str): Drug Analysis Filter: clinical | discovery
        out_file (str): Path to the output file (TSV)

    Return:
        DataFrame of the main statistic of the samples of the same study
        
    """

    # Load Samplesheet
    try:
       sheet = pd.read_table(samplesheet, sep=',')
    except Exception as e:
        raise ValueError(f"Error reading samplesheet file: {e}")

    study_filt = sheet.loc[sheet['study'] == study, ['sample_id', 'sex', 'tumour_content', 'cnas']]
    samplesheet_stats = study_filt[['sample_id', 'sex', 'tumour_content']]


    # COUNT MUTATIONS
    # Search all files that match the mut pattern
    muts_path = f"{mut_dir}/*.tsv"
    mut_files = glob.glob(muts_path)
    
    # Create list to build the final dataframe
    sampleids = []
    mut_counts = []
    
    for file in mut_files:
        file_name = os.path.basename(file)
        sampleid = file_name.split('_prep.mut.tsv')[0]
        try:
            sample_df = pd.read_table(file)
        except Exception as e:
            raise ValueError(f"Error reading mutation file of {sampleid} of study {study}: {e}")
            
        mut_count = sample_df.shape[0]
        
        sampleids.append(sampleid)
        mut_counts.append(mut_count)
    
    mut_df = pd.DataFrame({'sample_id': sampleids, 'mutations': mut_counts})


    # COUNT CNVS
    cnas = study_filt[['sample_id', 'cnas']]
    try:
        cna_df = pd.DataFrame(
           [(sampleid, pd.read_table(file).shape[0]) for sampleid, file in cnas.values],
           columns=['sample_id', 'cnas'])
    except Exception as e:
        raise ValueError(f"Error reading cna file of sample {sampleid} of study {study}: {e}")
        
    # COUNT INTERSECTS
    try:
        intersect = pd.read_table(intersect_combined)
    except Exception as e:
        raise ValueError(f"Error reading intersect file of study {study}: {e}")
        
    intersect_counts = intersect.groupby('sample_id').size().to_frame('intersect')

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
    
    # CREATE PANEL DATAFRAME
    panel = (samplesheet_stats
             .merge(mut_df, on='sample_id', how='inner')
             .merge(cna_df, on='sample_id', how='inner')
             .merge(intersect_counts, on='sample_id', how='inner')
             .merge(drug_stats, on='sample_id', how='inner'))

    
    panel.to_csv(out_file, sep='\t', index=False)
    
    return panel



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--study", action='store', required=True)
    input_parser.add_argument("--samplesheet", action='store', required=True)
    input_parser.add_argument("--mut_dir", action='store', required=True)
    input_parser.add_argument("--intersect_combined", action='store', required=True)
    input_parser.add_argument("--drug_prioritization", action='store', required=True)
    input_parser.add_argument("--drug_filter", action='store', required=True)
    input_parser.add_argument("--out_file", action='store', required=True)

    args = input_parser.parse_args()


    get_study_panels(args.study, args.samplesheet, args.mut_dir, args.intersect_combined, args.drug_prioritization, args.drug_filter, args.out_file)
