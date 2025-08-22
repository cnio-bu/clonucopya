import pandas as pd
import glob
import os
import argparse



def get_study_panels(study, samplesheet, mut_dir, intersect_combined, out_file):

    """
    Build Daframe for sample panels of the study.

    Args:
        study (str): name of the study
        samplesheet (str): Path to the samplesheet (CSV)
        mut_dir (str): Path to the files from mutation_prep study
        intersect_combined(str): Path to the intersect combined of the study
        out_file (str): Path to the output file (TSV)

    Return:
        DataFrame of the main statistic of the samples of the same study
        
    """

    # Load Samplesheet
    try:
       sheet = pd.read_table(samplesheet, sep=',')
    except Exception as e:
        raise ValueError(f"Error reading samplesheet file: {e}")

    study_filt = sheet.loc[sheet['study'] == study, ['sample_id', 'sex', 'tumour_content', 'cnvs']]
    samplesheet_stats = study_filt[['sample_id', 'sex', 'tumour_content']]


    # COUNT MUTATIONS
    # Search all files that match the mut pattern
    muts_path = f"{mut_dir}/*.tsv"
    mut_files = glob.glob(muts_path)
    
    # Crear el diccionario
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
    cnvs = study_filt[['sample_id', 'cnvs']]
    try:
        cnv_df = pd.DataFrame(
           [(sampleid, pd.read_table(file).shape[0]) for sampleid, file in cnvs.values],
           columns=['sample_id', 'cnvs'])
    except Exception as e:
        raise ValueError(f"Error reading cnv file of sample {sampleid} of study {study}: {e}")
        
    # COUNT INTERSECTS
    try:
        intersect = pd.read_table(intersect_combined)
    except Exception as e:
        raise ValueError(f"Error reading intersect file of study {study}: {e}")
        
    intersect_counts = intersect.groupby('sample_id').size().to_frame('intersect')

    # CREATE PANEL DATAFRAME

    panel = (samplesheet_stats
             .merge(mut_df, on='sample_id', how='inner')
             .merge(cnv_df, on='sample_id', how='inner')
             .merge(intersect_counts, on='sample_id', how='inner'))

    
    panel.to_csv(out_file, sep='\t', index=False)
    
    return panel



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--study", action='store', required=True)
    input_parser.add_argument("--samplesheet", action='store', required=True)
    input_parser.add_argument("--mut_dir", action='store', required=True)
    input_parser.add_argument("--intersect_combined", action='store', required=True)
    input_parser.add_argument("--out_file", action='store', required=True)

    args = input_parser.parse_args()


    get_study_panels(args.study, args.samplesheet, args.mut_dir, args.intersect_combined, args.out_file)
