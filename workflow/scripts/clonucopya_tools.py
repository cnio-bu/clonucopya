def chr_to_num(chr_str):
    """Convierte el cromosoma a un número para ordenamiento"""
    chr_str = chr_str.replace('chr', '')
    if chr_str == 'X':
        return 23
    elif chr_str == 'Y':
        return 24
    else:
        return int(chr_str)



def is_file_empty(file):
    try:
        df = pd.read_csv(file)
        return df.empty
    except Exception:
        return True



def json_to_csv(json_origin, csv_destination):

    """
    Convert PanDrugs JSON output  and extracts the geneDrugInfo data into a CSV format.

    Args:
        json_origin (str): Path to the input JSON file containing PanDrugs results
        csv_destination (str): Path where the output CSV file will be written

    Returns:
        None
    """
