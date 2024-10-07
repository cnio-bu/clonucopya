from __future__ import print_function, division, absolute_import
import logging

import pandas as pd
from typechecks import require_string
import numpy as np

TCGA_PATIENT_ID_LENGTH = 12

MAF_COLUMN_NAMES = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "Match_Norm_Seq_Allele1",
    "Match_Norm_Seq_Allele2",
    "Tumor_Validation_Allele1",
    "Tumor_Validation_Allele2",
    "Match_Norm_Validation_Allele1",
    "Match_Norm_Validation_Allele2",
    "Verification_Status",
    "Validation_Status",
    "Mutation_Status",
    "Sequencing_Phase",
    "Sequence_Source",
    "Validation_Method",
    "Score",
    "BAM_File",
    "Sequencer",
    "Tumor_Sample_UUID",
    "Matched_Norm_Sample_UUID",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "Exon_Number",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
    "all_effects",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "One_Consequence",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "ALLELE_NUM",
    "DISTANCE",
    "TRANSCRIPT_STRAND",
    "SYMBOL",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "BIOTYPE",
    "CANONICAL",
    "CCDS",
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "RefSeq",
    "SIFT",
    "PolyPhen",
    "EXON",
    "INTRON",
    "DOMAINS",
    "GMAF",
    "AFR_MAF",
    "AMR_MAF",
    "ASN_MAF",
    "EAS_MAF",
    "EUR_MAF",
    "SAS_MAF",
    "AA_MAF",
    "EA_MAF",
    "CLIN_SIG",
    "SOMATIC",
    "PUBMED",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "IMPACT",
    "PICK",
    "VARIANT_CLASS",
    "TSL",
    "HGVS_OFFSET",
    "PHENO",
    "MINIMISED",
    "ExAC_AF",
    "ExAC_AF_Adj",
    "ExAC_AF_AFR",
    "ExAC_AF_AMR",
    "ExAC_AF_EAS",
    "ExAC_AF_FIN",
    "ExAC_AF_NFE",
    "ExAC_AF_OTH",
    "ExAC_AF_SAS",
    "GENE_PHENO",
    "FILTER",
    "CONTEXT",
    "src_vcf_id",
    "tumor_bam_uuid",
    "normal_bam_uuid",
    "case_id",
    "GDC_FILTER",
    "COSMIC",
    "MC3_Overlap",
    "GDC_Validation_Status"
]


def load_maf_dataframe(path, nrows=None, raise_on_error=True, encoding=None):

    require_string(path, "Path to MAF")

    n_basic_columns = len(MAF_COLUMN_NAMES)

    df = pd.read_csv(
        path,
        comment="#",
        sep="\t",
        low_memory=False,
        skip_blank_lines=True,
        header=0,
        nrows=nrows,
        encoding=encoding)

    if len(df.columns) < n_basic_columns:
        error_message = (
            "Too few columns in MAF file %s, expected %d but got  %d : %s" % (
                path, n_basic_columns, len(df.columns), df.columns))
        if raise_on_error:
            raise ValueError(error_message)
        else:
            logging.warn(error_message)

    for expected, actual in zip(MAF_COLUMN_NAMES, df.columns):
        if expected != actual:
            if expected.lower() == actual.lower():
                df[expected] = df[actual]
                del df[actual]
            else:
                error_message = (
                    "Expected column %s but got %s" % (expected, actual))
                if raise_on_error:
                    raise ValueError(error_message)
                else:
                    logging.warn(error_message)

    return df


if __name__ == '__main__':
    chromosome = []
    start = []
    end = []
    allele = []
    strand = []
    identifier = []

    maf_df = load_maf_dataframe("/home/cleon/extdata/TCGA-38-4629/acc.mutect.maf.TCGA_38_4629_01A_02D_1265_08.tsv")

    arr_barcodes = np.array(maf_df["Tumor_Sample_Barcode"].unique())
    for idx, x in np.ndenumerate(arr_barcodes):
        df_barcode = pd.DataFrame(maf_df.loc[maf_df["Tumor_Sample_Barcode"] == x])

        for index, row in df_barcode.iterrows():
            chromosome.append(str(row["Chromosome"]))
            start.append(row["Start_Position"])
            end.append(row["End_Position"])
            allele_tmp = row["Reference_Allele"]
            if row["Reference_Allele"] != row["Tumor_Seq_Allele1"]:
                allele_tmp = allele_tmp + "/" + str(row["Tumor_Seq_Allele1"])
            if row["Reference_Allele"] != row["Tumor_Seq_Allele2"]:
                allele_tmp = allele_tmp + "/" + str(row["Tumor_Seq_Allele2"])
            allele.append(allele_tmp)
            strand.append(row["Strand"])
            # identifier puedes añadir la info de profundidad, número de lecturas de referencia y variante, ...
            identifier.append(str(row["t_ref_count"]) + " " + str(row["t_alt_count"]))

    '''# iterate through each row and select
    for index, row in maf_df.iterrows():
        chromosome.append(str(row["Chromosome"]))
        start.append(row["Start_Position"])
        end.append(row["End_Position"])
        allele_tmp = row["Reference_Allele"]
        if row["Reference_Allele"] != row["Tumor_Seq_Allele1"]:
            allele_tmp = allele_tmp + "/" + str(row["Tumor_Seq_Allele1"])
        if row["Reference_Allele"] != row["Tumor_Seq_Allele2"]:
            allele_tmp = allele_tmp + "/" + str(row["Tumor_Seq_Allele2"])
        allele.append(allele_tmp)
        strand.append(row["Strand"])
        # identifier puedes añadir la info de profundidad, número de lecturas de referencia y variante, ...
        identifier.append(str(row["t_ref_count"]) + " " + str(row["t_alt_count"]))'''

    df_output = pd.DataFrame({'chromosome': chromosome,
                              'start': start,
                              'end': end,
                              'allele': allele,
                              'strand': strand,
                              'identifier': identifier})

    # saving the dataframe
    df_output.to_csv('/home/cleon/extdata/TCGA-38-4629/acc.mutect.maf.TCGA_38_4629_01A_02D_1265_08.vep_default.tsv',
                     header=False,
                     sep='\t',
                     index=False)
