import os
import subprocess
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

def open_in_xmgrace(filename):
    """
    Opens an .xvg file in xmgrace using subprocess.
    """
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return
    
    try:
        # Run xmgrace with the specified .xvg file
        print(f"Opening {filename} in xmgrace...")
        subprocess.run(["xmgrace", filename], check=True)
    except FileNotFoundError:
        print("Error: xmgrace is not installed or not found in your PATH.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running xmgrace: {e}")

def plot_and_save_xvg_data(data, filename, output_dir):
    """
    Plots the data from an .xvg file and saves it as an image file.
    """
    if not data:
        return
    
    # Convert data to a list of columns
    columns = list(zip(*data))
    
    # The first column is usually the x-axis (e.g., time)
    x = columns[0]
    
    # Plot each subsequent column against the x-axis
    num_plots = len(columns) - 1  # Number of y-columns to plot
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # If there are too many columns, use subplots for clarity
    if num_plots > 5:
        fig, axes = plt.subplots(num_plots, 1, figsize=(10, 6 * num_plots), sharex=True)
        for i, ax in enumerate(axes, start=1):
            ax.plot(x, columns[i], label=f"Column {i}")
            ax.set_title(f"Plot of Column {i} from {os.path.basename(filename)}")
            ax.set_ylabel("Y-axis (Values)")
            ax.legend()
            ax.grid(True)
        axes[-1].set_xlabel("X-axis (e.g., Time)")
    else:
        plt.figure(figsize=(10, 6))
        for i, column in enumerate(columns[1:], start=1):
            plt.plot(x, column, label=f"Column {i}")
        
        # Add labels and title
        plt.title(f"Plot of {os.path.basename(filename)}")
        plt.xlabel("X-axis (e.g., Time)")
        plt.ylabel("Y-axis (Values)")
        plt.legend()
        plt.grid(True)
    
    # Save the plot as an image file
    output_filename = os.path.join(output_dir, os.path.splitext(os.path.basename(filename))[0] + ".png")
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()  # Close the figure to free memory
    print(f"Saved plot: {output_filename}")

def main():
    # List of .xvg files to process
    xvg_files = [
        "density.xvg",
        "gyrate.xvg",
        "potential.xvg",
        "pressure.xvg",
        "rmsd_xtal.xvg",
        "rmsd.xvg",
        "temperature.xvg"
    ]

    # Directory to save the plots
    output_directory = "xvg_plots"

    # Process each file
    for xvg_file in xvg_files:
        if os.path.exists(xvg_file):
            print(f"Processing {xvg_file}...")
            
            # Step 1: Open in xmgrace
            open_in_xmgrace(xvg_file)
            
            # Step 2: Read the data and save the plot
            data = read_xvg_file(xvg_file)
            if data:
                plot_and_save_xvg_data(data, xvg_file, output_directory)
        else:
            print(f"File {xvg_file} does not exist.")

if __name__ == "__main__":
    main()