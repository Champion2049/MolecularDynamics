import os
import matplotlib.pyplot as plt

def read_xvg_file(filename):
    """
    Reads an .xvg file, skipping comment lines, and returns the data as a list of lists.
    """
    data = []
    try:
        with open(filename, 'r') as file:
            for line in file:
                # Skip comment lines starting with '#' or '@'
                if line.startswith('#') or line.startswith('@'):
                    continue
                # Split the line into columns and convert them to floats
                columns = line.strip().split()
                data.append([float(col) for col in columns])
        return data
    except FileNotFoundError:
        print(f"File {filename} not found.")
        return None

def plot_rmsd_comparison(rmsd_data, rmsd_xtal_data):
    """
    Plots RMSD and RMSD_xtal on the same graph.
    """
    if not rmsd_data or not rmsd_xtal_data:
        print("One or both datasets are missing.")
        return
    
    # Extract x and y values from the data
    rmsd_x, rmsd_y = zip(*rmsd_data)
    rmsd_xtal_x, rmsd_xtal_y = zip(*rmsd_xtal_data)
    
    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.plot(rmsd_x, rmsd_y, label="RMSD (Equilibrium)", color="blue")
    plt.plot(rmsd_xtal_x, rmsd_xtal_y, label="RMSD (Crystal)", color="red", linestyle="--")
    
    # Add labels, title, and legend
    plt.title("Comparison of RMSD (Equilibrium vs. Crystal)")
    plt.xlabel("Time (ps)")
    plt.ylabel("RMSD (nm)")
    plt.legend()
    plt.grid(True)
    
    # Show the plot
    plt.tight_layout()
    plt.show()

def main():
    # Define the filenames for RMSD and RMSD_xtal
    rmsd_file = "rmsd.xvg"
    rmsd_xtal_file = "rmsd_xtal.xvg"

    # Read the data from the files
    rmsd_data = read_xvg_file(rmsd_file)
    rmsd_xtal_data = read_xvg_file(rmsd_xtal_file)

    # Check if both files were successfully read
    if rmsd_data and rmsd_xtal_data:
        # Plot the comparison
        plot_rmsd_comparison(rmsd_data, rmsd_xtal_data)
    else:
        print("Could not read one or both RMSD files.")
plt.savefig("rmsd_comparison.png", dpi=1200)

if __name__ == "__main__":
    main()