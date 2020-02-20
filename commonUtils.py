import csv


def getElementAbv(z):
    with open('PeriodicTable.csv', newline='') as csvfile:
        elements = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in elements:
            if row[0] == str(z):
                abv = row[2]
    return abv


def column(matrix, index):
    # Simple function to be able to pull out a single column of data from a list of lists, this lets us plot easily
    # as we can extract out our x,y values
    mycolumn = [i[index] for i in matrix]
    return mycolumn
