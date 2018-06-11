import csv

import pickle

CSV_EXT = ".csv"


def save_in_file(file_name, data):
    with open(file_name + CSV_EXT, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for line in data:
            writer.writerow(line)


def read_from_file(file_name):
    data = []
    with open(file_name + CSV_EXT, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for line in reader:
            data.append(line)

    return data


def save_results(filename, results):
    with open(filename + '.res', 'wb') as f:
        pickle.dump(results, f)

def load_results(filename):
    with open(filename + '.res', 'rb') as f:
        return pickle.load(f)
