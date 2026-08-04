import pandas as pd
import argparse
import os


def apply_alias(pvi_path, metadata_path):
    pvi = pd.read_csv(pvi_path, sep="\t")

    if os.path.exists(metadata_path):
        meta = pd.read_csv(metadata_path)
        alias_dict = dict(zip(meta["sample_id"], meta["alias"]))

        pvi["sample_id"] = pvi["sample_id"].map(alias_dict).fillna(pvi["sample_id"])
    
    return pvi


def pvi_filtering(pvi_aliased, just_snv, output_path):

    # Indels filter
    if just_snv:
        alleles = pvi_aliased['mutation_id'].str.split(':')
        ref_allele = alleles.str[2].str.strip()
        alt_allele = alleles.str[3].str.strip()
        pvi_aliased = pvi_aliased[
            (ref_allele.str.len() == 1) &
            (alt_allele.str.len() == 1)
        ]

    # Replace character of absent allele from . to - 
    mutation_map = {}
    for mut in pvi_aliased['mutation_id'].unique():
        try:
            chrom, pos, ref, alt = mut.split(':')
        except ValueError:
            print(f"Warning: malformed mutation_id skipped: {mut}")
            continue
    
        new_ref = '-' if ref == '.' else ref
        new_alt = '-' if alt == '.' else alt
    
        if new_ref != ref or new_alt != alt:
            mutation_map[mut] = f"{chrom}:{pos}:{new_ref}:{new_alt}"
    
    pvi_aliased['mutation_id'] = pvi_aliased['mutation_id'].replace(mutation_map)
                

    # Remove absent data at major_cn, minor_cn, and normal_cn columns
    pvi_mut_filt = pvi_aliased.dropna(subset=['major_cn', 'minor_cn', 'normal_cn']).copy()


    pvi_mut_filt.to_csv(output_path, sep="\t", index=False)


if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--pvi_input", action='store', required=True)
    input_parser.add_argument("--metadata", action='store', required=True)
    input_parser.add_argument("--output_file", action='store', required=True)
    input_parser.add_argument("--just_snv", type=lambda x: x.lower() in ('true', '1', 'yes'), required=True)
    

    args = input_parser.parse_args()

    pvi_alias = apply_alias(args.pvi_input, args.metadata)
    pvi_filtering(pvi_alias, args.just_snv, args.output_file)