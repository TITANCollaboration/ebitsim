import csv

def getElementAbv(z):
    with open('PeriodicTable.csv', newline='') as csvfile:
        elements = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in elements:
            if row[0] == str(z):
                abv = row[2]
    return abv


