import pandas as pd
import csv


def catalog_maker(row):
    with open('catalog_exmp1.csv', 'a') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(row)

    csvFile.close()
    



