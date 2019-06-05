import pandas as pd
import csv


def catalog_maker(row):
    with open('/home/ranjan/python/edibles/catalog/catalog_exmp.csv', 'a') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(row)

    csvFile.close()
    



