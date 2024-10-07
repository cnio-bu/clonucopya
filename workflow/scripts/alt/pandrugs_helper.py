import requests

# defining the api-endpoint
API_ENDPOINT_REGISTRATION = "https://www.pandrugs.org/pandrugs-backend/api/registration/"

data = '''{ 'password': 'TFM', \
          'email': 'clramos314@gmail.com', \
          'login': 'clramos314', \
          'uuid': 'clramos314' }'''

# sending post request and saving response as response object
r = requests.post(url=API_ENDPOINT_REGISTRATION, data=data)

# extracting response text
pastebin_url = r.text
print("The pastebin URL is:%s"%pastebin_url)
