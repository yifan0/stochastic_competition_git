import argparse
import os
import csv
import matplotlib.pyplot as plt

def read_csv_file(csv_file):
    data = {}
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if i < 6:
                continue  # Ignore the first 6 lines
            if i == 6:
                labels = row[2:]  # Labels for columns
                for label in labels:
                    data[label] = []
            else:
                # row[0] is the run number, this graph is for a single run
                step = int(row[1])
                for j, value in enumerate(row[2:]):
                    label = labels[j]
                    data[label].append((step, float(value)))
    return data

def plot_data(data):
    for label, values in data.items():
        steps, values = zip(*values)
        plt.figure(figsize=(10, 6))
        plt.plot(steps, values, label=label)
        plt.xlabel('Thousands of steps')
        plt.ylabel(label)
        plt.title(f'Plot of {label}')
        plt.grid(True)
        plt.legend()
        plt.savefig(f'{label}.png')

def main():
    parser = argparse.ArgumentParser(description='Generate plots from CSV data')
    parser.add_argument('csv_file', type=str, help='CSV file name')
    args = parser.parse_args()

    if not os.path.isfile(args.csv_file):
        print(f"Error: The file '{args.csv_file}' does not exist.")
        return

    data = read_csv_file(args.csv_file)
    plot_data(data)

if __name__ == "__main__":
    main()
