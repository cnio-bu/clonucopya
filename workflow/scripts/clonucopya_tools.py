import json
import csv

def chr_to_num(chr_str):
    """Convert X and Y chromosome into a number to sort chromosomes from 1 to 24"""
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

    # Opening JSON file and loading the data
    with open(fr'{json_origin}') as json_file:
            data = json.load(json_file)

    geneDrugGroup = data['geneDrugGroup']
    data_file = open(fr'{csv_destination}', 'w')
    csv_writer = csv.writer(data_file)
    count = 0

    for gdg in geneDrugGroup:
        geneDrugInfo = gdg['geneDrugInfo']
        for gdi in geneDrugInfo:
            if count == 0:
                # Writing headers of CSV file
                header = gdi.keys()
                csv_writer.writerow(header)
                count += 1

            # Writing data of CSV file
            csv_writer.writerow(gdi.values())

    data_file.close()
