import argparse
import pandas as pd
import pysam
from pathlib import Path


def gather_study_mutations(mut_path):
    mut_dir = Path(mut_path)
    try:
        dfs = []
        for mut_file in mut_dir.glob("*_prep.mut.tsv"):
            df = pd.read_csv(mut_file, sep="\t", usecols=["CHROM", "POS", "REF", "ALT"])
            dfs.append(df)
    except Exception as e:
        raise ValueError(f"Error reading mutation prep files: {e}")

    mut_master = (
        pd.concat(dfs, ignore_index=True)
        .drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
        .sort_values(["CHROM", "POS"])
        .reset_index(drop=True)
    )
    return mut_master


def get_bam_files(bam_path):
    bam_dir = Path(bam_path)
    bam_dict = {
        bam_file.stem: bam_file
        for bam_file in bam_dir.glob("*.bam")
    }
    return bam_dict


def get_counts(bam, chrom, pos, ref, alt, min_bq=20, min_mq=20):   
    ref_count, alt_count = 0, 0
    for pileup_col in bam.pileup(chrom, pos - 1, pos,
                                  min_base_quality=min_bq,
                                  min_mapping_quality=min_mq,
                                  truncate=True):
        for read in pileup_col.pileups:
            if not read.is_del and not read.is_refskip:
                base = read.alignment.query_sequence[read.query_position]
                if base == ref:
                    ref_count += 1
                elif base == alt:
                    alt_count += 1
    return ref_count, alt_count


def process_bam_info(bam_dict, mut_master, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_dfs = {}

    for sample, bam_path in bam_dict.items():
        df = mut_master.copy()

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            counts = df.apply(
                lambda row: pd.Series(
                    get_counts(bam, row["CHROM"], row["POS"], row["REF"], row["ALT"]),
                    index=["ref_counts", "alt_counts"]
                ),
                axis=1
            )

        df[["ref_counts", "alt_counts"]] = counts
        df["VAF"] = df["alt_counts"] / (df["ref_counts"] + df["alt_counts"])
        df["mutation_id"] = (                         
            df["CHROM"].astype(str) + ":" +
            df["POS"].astype(str) + ":" +
            df["REF"].astype(str) + ":" +
            df["ALT"].astype(str)
        )
        sample_dfs[sample] = df

    for sample, df in sample_dfs.items():
        df.to_csv(output_dir / f"{sample}_check.mut.tsv", sep="\t", index=False)


if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()          
    input_parser.add_argument("--mut_path", action='store', required=True)
    input_parser.add_argument("--bam_path", action='store', required=True)
    input_parser.add_argument("--output_dir", action='store', required=True)
    args = input_parser.parse_args()

    mut_master = gather_study_mutations(args.mut_path)
    bam_dict = get_bam_files(args.bam_path)
    process_bam_info(bam_dict, mut_master, args.output_dir)
