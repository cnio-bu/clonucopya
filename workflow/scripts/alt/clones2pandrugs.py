import io
import sys
import numpy as np
import pandas as pd
import time
import requests
import re
from requests.auth import HTTPBasicAuth

sample_id = ''
frames = []
dict_computations = {}


def read_vcf(path):
    d_vcf = {}
    arr_headers = []
    lines = []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('##'):
                lines.append(line)
            else:
                arr_headers.append(line.replace('\n', ''))

    d_vcf['headers'] = pd.DataFrame(data=arr_headers)

    d_vcf['body'] = pd.read_csv(io.StringIO(''.join(lines)),
                                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                                       'QUAL': str, 'FILTER': str, 'INFO': str},
                                sep='\t'
                                )
    return d_vcf


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

    sample_id = sys.argv[1]

    dict_final = {}
    dict_vcf = read_vcf(sys.argv[2])
    header_vcf = dict_vcf['headers']
    df_filtered = dict_vcf['body']

    '''df_clusters = pd.read_csv(f'/home/cleon/extdata/tumor_{sample_id}_clustered.tsv',
                              sep="\t",
                              low_memory=False)'''
    df_clusters = pd.read_csv(sys.argv[3],
                              sep="\t",
                              low_memory=False)

    arr = np.array(df_clusters['cluster_id'].unique())
    for idx, x in np.ndenumerate(arr):
        df_cluster = pd.DataFrame(df_clusters.loc[df_clusters['cluster_id'] == x])

        for index, row in df_cluster.iterrows():
            mutation_id = row['mutation_id']
            chro = mutation_id.split(':')[1]
            pos = mutation_id.split(':')[2]

            rslt_df = df_filtered.loc[df_filtered['POS'] == int(pos)]
            if pd.DataFrame(rslt_df).size > 0:
                frames.append(rslt_df)

        param_output_split = str(sys.argv[4]).split('/')
        param_output_split[len(param_output_split) - 1] = '{}.{}.vcf'.format(sample_id, x)
        outputfile = '/'.join(param_output_split)

        '''header_vcf.to_csv(f'/home/cleon/extdata/input_pandrugs_{sample_id}.vcf',
                          sep='\t',
                          index=False,
                          header=False)'''
        '''header_vcf.to_csv(outputfile,
                          sep='\t',
                          index=False,
                          header=False)'''
        if len(frames) > 0:
            df_frames = pd.concat(frames)
            '''df_frames.to_csv(f'/home/cleon/extdata/input_pandrugs_{sample_id}.vcf',
                             sep='\t',
                             index=False,
                             header=True,
                             mode='a')'''
            df_frames.to_csv(outputfile,
                             sep='\t',
                             index=False,
                             header=True,
                             mode='w')

            # #######################################################################################
            # Requests to Pandrugs individually
            # #######################################################################################

            # POST - computation
            # ff = open(f'/home/cleon/extdata/input_pandrugs_{sample_id}.vcf', 'r')
            ff = open(outputfile, 'r')
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
                dict_computations[str(x)] = location
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
                    url = 'https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/files/guest/{}/vscorefile/'.format(computation_id)

                    headers = {
                        'authorization': "Basic Z3Vlc3Q6Z3Vlc3Q=",
                        'cache-control': "no-cache",
                        'Accept': "application/octet-stream"
                    }

                    print('(GET) Download of the VScore file --> {}'.format(url))
                    r = requests.get(url, stream=True, allow_redirects=False)
                    sta_code = r.status_code
                    print('      Download of the VScore file --> response.status_code = {}'.format(r.status_code))

                    param_output_split[len(param_output_split) - 1] = str(get_filename_from_cd(r.headers.get('Content-Disposition'))).replace('"', '')
                    filename = '/'.join(param_output_split)
                    print('      Download of the VScore file --> {}'.format(filename))
                    f_output = open(filename, 'wb').write(r.content)

                df_computations = pd.DataFrame.from_dict(dict_computations, orient='index', columns=['url_computation'])
                df_computations.to_csv(str(sys.argv[4]), sep='\t', index=True)
