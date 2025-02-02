import requests
import os
import argparse
import time
import pandas as pd
import json
import csv


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



def pandrugs_query(vep_vcf,output_dir):
    """
    Query PanDrugs database (https://www.pandrugs.org) with a VEP-annotated VCF file and retrieve drug-gene interactions.

    Args:
        vep_vcf (str): Path to the input VEP-annotated VCF file
        output_dir (str): Directory where output files will be stored

    Returns:
        None

    Outputs:
        - {sample_id}_vscore.vcf: Variant scores file
        - {sample_id}_gene-drug.json: Drug-gene interactions in JSON format
        - {sample_id}_gene-drug.csv: Drug-gene interactions in CSV format
        - {sample_id}_computation.tsv: Computation tracking information
    """

    if os.path.exists(vep_vcf):
        os.makedirs(output_dir, exist_ok=True)
        sample_id = os.path.splitext(os.path.basename(vep_vcf))[0]
        dict_computations = {}

        # POST - computation
        ff = open(vep_vcf, 'r')
        data = ff.read()
        ff.close()

        url_post_computation = "https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/guest"

        querystring = {"name": "TestOfComputation"}

        payload = data
        headers = {
            'content-type': "text/plain",
            'authorization': "Basic Z3Vlc3Q6Z3Vlc3Q=",
            'cache-control': "no-cache"
        }

        response = requests.request("POST", url_post_computation, data=payload, headers=headers, params=querystring)

        print('variantsanalysis --> status_code = {}'.format(response.status_code))

        if response.status_code == 201:
            location = response.headers.get('location')
            print('variantsanalysis --> location = {}'.format(location))
            dict_computations[str(sample_id)] = location
            url_variant_analysis = location

            url_words = url_variant_analysis.split('/')
            computation_id = url_words[len(url_words) - 1]

            headers = {
                'authorization': "Basic Z3Vlc3Q6Z3Vlc3Q=",
                'cache-control': "no-cache"
            }

            sta_code = 500

            # GET - Status generation of the VScore file
            url_status = 'https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/guest/'
            exit_loop = False
            finished = False
            failed = False

            while not exit_loop:
                try:
                    print('(GET) Status generation of the VScore file --> {}'.format(location))
                    response = requests.request("GET", location, headers=headers)
                    sta_code = response.status_code
                    print('      Status generation of the VScore file --> status_code = {}'.format(response.status_code))
                    r_json = response.json()
                    failed = r_json["failed"]
                    finished = r_json["finished"]
                    print('      Status generation of the VScore file --> failed={} and finished={}'.format(failed, finished))

                    if sta_code == 401 or sta_code == 403:
                        time.sleep(3)
                    elif sta_code == 500:
                        time.sleep(3)
                    elif sta_code == 200:
                        if finished:
                            exit_loop = True
                        else:
                            time.sleep(3)
                except requests.exceptions.RequestException as error:
                    print("Error: ", error)

            if sta_code == 200 and finished and not failed:

                # GET - VScore file
                url = 'https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/files/guest/{}/vscorefile/'.format(
                    computation_id)

                headers = {
                    'authorization': "Basic Z3Vlc3Q6Z3Vlc3Q=",
                    'cache-control': "no-cache",
                    'Accept': "application/octet-stream"
                }

                print('(GET) Download of the VScore file --> {}'.format(url))
                r = requests.get(url, stream=True, allow_redirects=False)
                sta_code = r.status_code
                print('      Download of the VScore file --> response.status_code = {}'.format(r.status_code))

                filename = os.path.join(output_dir,  f"{sample_id}_vscore.vcf")
                print('      Download of the VScore file --> {}'.format(filename))
                f_output = open(filename, 'wb').write(r.content)

                # GET - drug-gene
                url = "https://www.pandrugs.org/pandrugs-backend/api/genedrug/fromComputationId?computationId={}&cancerDrugStatus=APPROVED&cancerDrugStatus=CLINICAL_TRIALS&nonCancerDrugStatus=CLINICAL_TRIALS&nonCancerDrugStatus=EXPERIMENTAL&nonCancerDrugStatus=APPROVED&biomarker=true&pathwayMember=true&directTarget=true".format(computation_id)

                headers = {
                    'Accept': "application/json",
                    'cache-control': "no-cache"
                }

                print('(GET) Download of the Gene-Drug file --> {}'.format(url))
                r = requests.get(url, stream=True, allow_redirects=False)
                sta_code = r.status_code
                print('      Download of the Gene-Drug file --> response.status_code = {}'.format(r.status_code))
                print(r)

                if sta_code == 200:
                    filename = os.path.join(output_dir, f"{sample_id}_gene-drug.json")
                    print('      Download of the Gene-drug file --> {}'.format(filename))
                    f_output = open(filename, 'wb').write(r.content)

                    # write 'geneDrugInfo' tags into a csv
                    json_to_csv(filename, filename.replace('json', 'csv'))
            df_computations = pd.DataFrame.from_dict(dict_computations, orient='index', columns=['url_computation'])
            df_computations.to_csv(os.path.join(output_dir, f"{sample_id}_computation.tsv"), sep='\t', index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep_vcf", required=True, help="Input VEP-annotated VCF file to analyze")
    parser.add_argument("--out_dir", required=True, help="Output directory for results (will be created if it doesn't exist)")

    args = parser.parse_args()

    pandrugs_query(args.vep_vcf, args.out_dir)
