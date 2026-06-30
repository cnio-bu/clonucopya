import pandas as pd
import argparse
import os


def apply_alias(pvi_path, metadata_path, output_path):
    pvi = pd.read_csv(pvi_path, sep="\t")

    if os.path.exists(metadata_path):
        meta = pd.read_csv(metadata_path)
        alias_dict = dict(zip(meta["sample_id"], meta["alias"]))

        pvi["sample_id"] = pvi["sample_id"].map(alias_dict).fillna(pvi["sample_id"])

    pvi.to_csv(output_path, sep="\t", index=False)


if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--pvi_input",    action='store', required=True)
    input_parser.add_argument("--metadata",     action='store', required=True)
    input_parser.add_argument("--output_file",  action='store', required=True)

    args = input_parser.parse_args()

    apply_alias(args.pvi_input, args.metadata, args.output_file)