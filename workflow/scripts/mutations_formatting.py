import pandas as pd
import argparse
import pysam


def get_vcf_columns(input_vcf):
    """Read CHROM line and returns all column names"""
    with open(input_vcf) as f:
        for line in f:
            if line.startswith('#CHROM'):
                return line.strip().lstrip('#').split('\t')
    raise ValueError("Header not found at VCF file")


def extract_format_field(genotype_series, format_series, field):
    """
    Extracts genotype field using FORMAT column.
    """
    results = []
    for gt, fmt in zip(genotype_series, format_series):
        keys = str(fmt).split(':')
        vals = str(gt).split(':')
        fmt_dict = dict(zip(keys, vals))
        results.append(fmt_dict.get(field, pd.NA))
    return pd.Series(results, index=genotype_series.index)

    

def process_vcf_mutations(input_vcf, just_snv, output_file, sample):
    """
    Process VCF file to extract mutation information and read counts.
    
    Args:
        input_vcf (str): Path to input VCF file
        just_snv (bool): Filters out indels if True
        output_file (str, optional): Path to output TSV file.
        sample (str): Sample name to extract ('tumor' by default, 'normal' also accepted).
                      If not found, uses the last sample column.
        
    Returns:
        pandas.DataFrame: Processed mutations data
    """
    # Check vcf integrity
    try:
        vcf = pysam.VariantFile(input_vcf)
        # Check every snv entry
        for rec in vcf.fetch():
            pass
        vcf.close()
    except (ValueError, OSError, RuntimeError) as e:
        raise ValueError(f"Truncated VCF file: {input_vcf} -> {e}")


    # Read header to get column names
    vcf_cols = get_vcf_columns(input_vcf)

    # Fixed VCF columns
    fixed_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    sample_cols = [c for c in vcf_cols if c not in fixed_cols]

    if not sample_cols:
        raise ValueError("No sample column found at VCF file.")

    # Select the target sample column
    if sample in sample_cols:
        target_sample = sample
    else:
        target_sample = sample_cols[-1]
        print(f"Sample '{sample}' not found. Using '{target_sample}'.")

    target_sample_idx = vcf_cols.index(target_sample)

    #Load VCF (skip comment lines)
    try:
        mut_vcf = pd.read_csv(
            input_vcf, sep='\t', comment='#', header=None,
            usecols=range(len(vcf_cols))
        )
        mut_vcf.columns = vcf_cols
    except Exception as e:
        raise ValueError(f"Error reading VCF file: {e}")

    # Select relevant columns
    mut_vcf_filt = mut_vcf[['CHROM', 'POS', 'REF', 'ALT']].copy()

    # Adapt format to VEP standard
    mut_vcf_filt['REF'] = mut_vcf_filt['REF'].str.replace('.', '-', regex=False)
    mut_vcf_filt['ALT'] = mut_vcf_filt['ALT'].str.replace('.', '-', regex=False)

    # Add genotype and FORMAT columns
    mut_vcf_filt['_genotype'] = mut_vcf[target_sample].astype(str).values
    mut_vcf_filt['_format']   = mut_vcf['FORMAT'].astype(str).values

    # Indels filter
    if just_snv:
        mut_vcf_filt = mut_vcf_filt[
            (mut_vcf_filt['REF'].str.len() == 1) &
            (mut_vcf_filt['ALT'].str.len() == 1)
        ].copy()

    # Extract read counts from FORMAT
    try:
        ad_raw = extract_format_field(mut_vcf_filt['_genotype'], mut_vcf_filt['_format'], 'AD')
        mut_vcf_filt['ref_counts'] = pd.to_numeric(
            ad_raw.str.split(',').str[0], errors='coerce')
        mut_vcf_filt['alt_counts'] = pd.to_numeric(
            ad_raw.str.split(',').str[1], errors='coerce')
    except Exception as e:
        raise ValueError(f"Error parsing read counts: {e}")

    # Drop support columns
    mut_vcf_filt = mut_vcf_filt.drop(columns=['_genotype', '_format'])

    # Build mutation ID
    mut_vcf_filt['mutation_id'] = (
        mut_vcf_filt['CHROM'].astype(str) + ':' +
        mut_vcf_filt['POS'].astype(str) + ':' +
        mut_vcf_filt['REF'] + ':' +
        mut_vcf_filt['ALT']
    )

    mut_vcf_filt.to_csv(output_file, sep='\t', index=False)

    return mut_vcf_filt


    
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input_vcf", action='store', required=True)
    input_parser.add_argument("--just_snv", type=lambda x: x.lower() in ('true', '1', 'yes'), required=True)
    input_parser.add_argument("--output_file", action='store', required=True)
    input_parser.add_argument("--sample", action='store', default="tumour")
    args = input_parser.parse_args() 
    
    process_vcf_mutations(args.input_vcf, args.just_snv, args.output_file, args.sample)
