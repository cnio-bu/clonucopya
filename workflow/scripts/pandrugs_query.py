import requests
import os
import argparse
import time
import pandas as pd
import json_to_csv


def pandrugs_query(vep_vcf,output_dir):
     dict_computations = {}
     sample_id = os.path.splitext(vep_vcf)[0]

     if os.path.exists(vep_vcf):
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
                    json_to_csv.write(filename, filename.replace('json', 'csv'))
            df_computations = pd.DataFrame.from_dict(dict_computations, orient='index', columns=['url_computation'])
            df_computations.to_csv(os.path.join(output_dir, f"{sample_id}_computation.tsv")),
                                           sep='\t',
                                           index=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep_dir", action='store', required=True)
    parser.add_argument("--out_dir", action='store', required=True)

    args = parser.parse_args()

    # Iterate over VCF files in the directory ESTO VA EN LA REGLA
    for clone in os.listdir(args.vep_dir):
        if clone.endswith('.vcf'):
            vcf_path = os.path.join(args.vep_dir, clone)
            pandrugs_query(vcf_path, args.out_dir)
