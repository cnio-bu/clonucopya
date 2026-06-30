from pyliftover import LiftOver
import pandas as pd
import argparse
import os
import re

chain = None

def liftover_id(mutation_id):
    chrom, pos, ref, alt = mutation_id.split(":")
    result = chain.convert_coordinate(f"chr{chrom}", int(pos) - 1)
    if result:
        new_chrom = result[0][0].replace("chr", "")
        new_pos = result[0][1] + 1
        return f"{new_chrom}:{new_pos}:{ref}:{alt}"
    return None



def get_hg38_id(pvi_input, reference, liftover, output_file):

    global chain

    pvi = pd.read_csv(pvi_input, sep="\t")

    # Upload Liftover file
    chain = LiftOver(liftover)
    
    if reference == 'GRCh38':
        pvi.to_csv(output_file, sep='\t', index=False)
     
    else:
        pvi["id38"] = pvi["mutation_id"].map(liftover_id)

        # Drop GRCh37 mutation_id version and replace it with GRCh38 upgraded ids
        pvi = pvi.drop(columns="mutation_id").rename(columns={"id38": "mutation_id"})

        # Reorder dataframe columns
        cols = ["mutation_id"] + [c for c in pvi.columns if c != "mutation_id"]
        pvi = pvi[cols]

        # Filter contigs or non-chromosomal mutations
        valid_chroms = re.compile(r'^(1[0-9]|2[0-2]|[1-9]|X|Y):')
        pvi = pvi[pvi['mutation_id'].notna() & pvi['mutation_id'].str.match(valid_chroms)]
        pvi.to_csv(output_file, sep='\t', index=False)



if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--pvi_input", action='store', required=True)
    input_parser.add_argument("--genome_reference", action='store', required=True, choices=["GRCh37", "GRCh38"], help="Genome used to anotate mutations (SNVs and CNAs).")
    input_parser.add_argument("--liftOver", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)

    args = input_parser.parse_args()

    get_hg38_id(args.pvi_input, args.genome_reference, args.liftOver, args.output_file)
