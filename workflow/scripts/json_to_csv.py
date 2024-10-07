import json
import csv


def write(json_origin, csv_destination):

	# Opening JSON file and loading the data
	# into the variable data
	with open(fr'{json_origin}') as json_file:
		data = json.load(json_file)

	geneDrugGroup = data['geneDrugGroup']

	# now we will open a file for writing
	data_file = open(fr'{csv_destination}', 'w')

	# create the csv writer object
	csv_writer = csv.writer(data_file)

	# Counter variable used for writing
	# headers to the CSV file
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
