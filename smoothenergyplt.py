import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def read_xvg(file_path):
    """Reads an .xvg file and returns the data as a NumPy array."""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(('@', '#')):
                continue  # Skip comment lines
            data.append([float(x) for x in line.split()])
    return np.array(data)

def smooth_data(data, window_length=51, polyorder=2, passes=2):
    """
    Smooths the y-values of the data using a Savitzky-Golay filter.
    :param data: A 2D NumPy array where the first column is x and the second is y.
    :param window_length: The length of the smoothing window (must be odd).
    :param polyorder: The order of the polynomial used for smoothing.
    :param passes: Number of smoothing passes to apply.
    :return: Smoothed y-values.
    """
    x = data[:, 0]
    y = data[:, 1]
    
    # Apply smoothing multiple times
    for _ in range(passes):
        y = savgol_filter(y, window_length=window_length, polyorder=polyorder)
    
    return x, y

def plot_and_save_data(original_x, original_y, smoothed_x, smoothed_y, output_png):
    """
    Plots the original and smoothed data and saves the plot as a PNG file.
    :param original_x: Original x-values.
    :param original_y: Original y-values.
    :param smoothed_x: Smoothed x-values.
    :param smoothed_y: Smoothed y-values.
    :param output_png: Path to save the PNG file.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(original_x, original_y, label='Original Data', color='lightgray', linewidth=1)
    plt.plot(smoothed_x, smoothed_y, label='Smoothed Data', color='red', linewidth=2)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Smoothed Curve')
    plt.legend()
    plt.grid(True)
    
    # Save the plot as a PNG file
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_png}")
    
    # Optionally display the plot
    plt.show()

# Main script
if __name__ == "__main__":
    # Path to your .xvg file
    input_file = "potential.xvg"
    
    # Read the data
    data = read_xvg(input_file)
    
    # Smooth the data with increased smoothing parameters
    smoothed_x, smoothed_y = smooth_data(data, window_length=51, polyorder=2, passes=2)
    
    # Save the smoothed data as a PNG file
    output_png = "smoothed_potential.png"
    plot_and_save_data(data[:, 0], data[:, 1], smoothed_x, smoothed_y, output_png)