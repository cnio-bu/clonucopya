# It is also valid for the AML's dataset
from pathlib import Path
import pandas as pd
import time
import requests
import re
import argparse
from requests.auth import HTTPBasicAuth
import numpy as np
import os

import json_to_csv

sample_id = ''
output_file = ''
dict_computations = {}
filename = ''
cluster_id = ''


def get_filename_from_cd(cd):
    """
    Get filename from content-disposition
    """
    if not cd:
        return None
    fname = re.findall('filename=(.+)', cd)
    if len(fname) == 0:
        return None
    return fname[0]


if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--clusters", nargs='+', action='append')
    input_parser.add_argument("--patient", action='store')
    input_parser.add_argument("--vep_dir", action='store')
    input_parser.add_argument("--out_dir", action='store')
    args = input_parser.parse_args()

    # i.e: /home/cleon/vep_data/18T59-I/tables/cluster.tsv /home/cleon/vep_data/16T320/tables/cluster.tsv
    clusters_tsv = args.clusters

    # i.e: '/home/cleon/vep_data'
    vep_dir = args.vep_dir

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/pandrugs'
    output_dir = args.out_dir

    arr = np.array(clusters_tsv)
    for idx_cluster_tsv, row_cluster_tsv in np.ndenumerate(arr):

        # i.e: /home/cleon/vep_data/18T59-I/tables/cluster.tsv
        cluster_tsv = pd.read_csv(row_cluster_tsv, sep='\t')
        arr = np.array(cluster_tsv['cluster_id'].unique())
        for idx, x in np.ndenumerate(arr):
            df_cluster = pd.DataFrame(cluster_tsv.loc[cluster_tsv['cluster_id'] == x])

            for index, row in df_cluster.iterrows():
                sample_id = row['sample_id']
                cluster_id = row['cluster_id']

                # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/16T320
                if not os.path.exists(output_dir + '/' + str(sample_id)):
                    os.makedirs(output_dir + '/' + str(sample_id))

                # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/16T320/0
                if not os.path.exists(output_dir + '/' + str(sample_id) + "/" + str(cluster_id)):
                    os.makedirs(output_dir + '/' + str(sample_id) + "/" + str(cluster_id))

                # i.e: ensembl_16T320.0.vcf
                vcf_file = 'ensembl_{}.{}.vcf'.format(sample_id, cluster_id)
                # print('>>>>>>>>>>>>>>{}'.format(vcf_file))

                input_param = vep_dir + "/" + sample_id + "/" + vcf_file

                # #######################################################################################
                # Requests to Pandrugs individually
                # #######################################################################################

                # POST - computation
                ff = open(input_param, 'r')
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
                            response = requests.request("GET", location, auth=HTTPBasicAuth('guest', 'guest'))
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

                        filename = str(output_dir) + '/{}/{}/{}'.format(sample_id, cluster_id, str(get_filename_from_cd(r.headers.get('Content-Disposition'))).replace('"', ''))
                        print('      Download of the VScore file --> {}'.format(filename))
                        f_output = open(filename, 'wb').write(r.content)

                        # GET - drug-gene
                        url = "https://www.pandrugs.org/pandrugs-backend/api/genedrug/fromComputationId?computationId={}&cancerDrugStatus=APPROVED&cancerDrugStatus=CLINICAL_TRIALS&nonCancerDrugStatus=CLINICAL_TRIALS&nonCancerDrugStatus=EXPERIMENTAL&nonCancerDrugStatus=APPROVED&biomarker=true&pathwayMember=true&directTarget=true".format(computation_id)

                        headers = {
                            'Accept': "application/json"
                        }

                        print('(GET) Download of the Gene-Drug file --> {}'.format(url))
                        r = requests.get(url, stream=True, allow_redirects=False)
                        sta_code = r.status_code
                        print('      Download of the Gene-Drug file --> response.status_code = {}'.format(r.status_code))
                        print(r)

                        if sta_code == 200:
                            filename = str(output_dir) + '/{}/{}/{}'.format(sample_id, cluster_id, 'gene-drug.json')
                            print('      Download of the Gene-drug file --> {}'.format(filename))
                            f_output = open(filename, 'wb').write(r.content)

                            # write 'geneDrugInfo' tags into a csv
                            json_to_csv.write(filename, filename.replace('json', 'csv'))

                    df_computations = pd.DataFrame.from_dict(dict_computations, orient='index', columns=['url_computation'])
                    df_computations.to_csv(str(output_dir + '/{}/{}/computation.tsv').format(sample_id, cluster_id),
                                           sep='\t',
                                           index=True)

