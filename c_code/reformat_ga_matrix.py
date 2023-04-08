import sys
import re
import numpy as np

# Get the file name from the command line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <file_name>")
    sys.exit(1)

file_name = sys.argv[1]

# Read in the file
try:
    with open(file_name, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"File {file_name} not found.")
    sys.exit(1)

# Find the first non-blank line
header_line = None
for line in lines:
    if line.strip() != "":
        header_line = line.strip()
        break

# Search for square brackets in the non-blank line
match = re.search(r'\[(\d+):(\d+),(\d+):(\d+)\]', header_line)

start_row = end_row = start_col = end_col = 0
if match:
    # Extract the integers from the square brackets
    start_row = int(match.group(1))
    end_row = int(match.group(2))
    start_col = int(match.group(3))
    end_col = int(match.group(4))
    
    print(f"Found header: [{start_row}:{end_row},{start_col}:{end_col}]")
else:
    print("No header found.")
    sys.exit(1)

matrix = np.zeros((end_row-start_row+1, end_col-start_col+1))

# Find start of data
first_line = 0
for i in range(len(lines)):
    if(len(lines[i].strip()) > 0 and lines[i].strip()[0].isdigit()):
        first_line = i
        break

# Look at file data
for i in range(first_line,len(lines),end_row-start_row+4):
    # i+0 is column labels
    columns = [eval(i) for i in lines[i].split()]
    # i+1 is dashes
    # i+2 through i+ end_row-start_row+3 is data
    for j in range(2,end_row-start_row+3):
        row = [eval(i) for i in lines[i+j].split()]
        row_idx = row[0]
        row = row[1:]
        for k in range(len(columns)):
            matrix[row_idx-1][columns[k]-1] = row[k]

np.savetxt(f"{file_name[:-4]}.csv", matrix, delimiter=",", fmt='%.5f')
print("Saved to", f"{file_name[:-4]}.csv")

