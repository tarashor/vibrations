import csv

data=[1, 2, 3]

with open('eggs.csv', 'w') as csvfile:
     writer = csv.writer(csvfile, delimiter=',')
     for line in data:
         writer.writerow(line)