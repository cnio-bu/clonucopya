import requests
import sys

sample_id = ''


if __name__ == '__main__':

    sample_id = sys.argv[1]

    # POST - computation
    # ff = open(f'/home/cleon/extdata/input_pandrugs_{sample_id}.vcf', 'r')
    ff = open(sys.argv[2], 'r')
    data = ff.read()
    ff.close()

    url_post_computation = "https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/guest"

    querystring = {"name": "TestOfComputation"}

    payload = data
    headers = {
        'content-type': "text/plain",
        'authorization': "Basic Z3Vlc3Q6Z3Vlc3Q=",
        'cache-control': "no-cache",
        'postman-token': "8a493f45-ec51-41d4-0fff-45cca212c478"
    }

    response = requests.request("POST", url_post_computation, data=payload, headers=headers, params=querystring)

    print(f'status={response.status_code}')

    if response.status_code == 201:
        # GET - variant analysis
        url_variant_analysis = response.headers.get('location')

        url_words = url_variant_analysis.split('/')
        computation_id = url_words[len(url_words)-1]

        headers = {
            'authorization': "Basic Z3Vlc3Q6Z3Vlc3Q=",
            'cache-control': "no-cache",
            'postman-token': "08d47c68-ce25-6491-7b78-6911acf45e4c"
        }

        sta_code = 500

        while sta_code != 200:
            response = requests.request("GET", url_variant_analysis, headers=headers)
            sta_code = response.status_code
            print(response.text)

        if response.status_code == 200:
            # GET - VScore file
            url = f'https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/files/guest/{computation_id}/vscorefile/'

            headers = {
                'authorization': "Basic Z3Vlc3Q6Z3Vlc3Q=",
                'cache-control': "no-cache",
                'postman-token': "85554b98-fb26-b7fa-a0a4-9ed6f4ce8674"
            }

            response = requests.request("GET", url, headers=headers)

            print(response.text)

            file = open(sys.argv[3], "wb")

            file.write(response.content)

            file.close()
