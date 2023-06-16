import sys
import matplotlib.pyplot as plt
import pandas as pd

# Check if a filename is provided
if len(sys.argv) < 2:
    print("Please provide a filename")
    sys.exit(1)

# Get the filename from the command line argument
filename = sys.argv[1]

# Read the CSV file using pandas
data = pd.read_csv(filename)

# Create a heatmap using Matplotlib
#plt.imshow(data, cmap='coolwarm', interpolation='antialiased')
plt.pcolor(data, cmap=plt.cm.seismic, vmin=0.8, vmax=1.2)
plt.colorbar()

# Set the title and labels

# Save the figure as a PNG file with the same name as the CSV file
plt.savefig(f"{filename[:-4]}.png")
print(f"Saved figure: {filename[:-4]}.png")

# Show the plot
plt.show()

