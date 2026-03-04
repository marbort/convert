import json
import argparse
import numpy as np
import matplotlib.pyplot as plt

def hist_avg_angle(input_file, output_file,coord,color=0,title="title",legend="legend",limy=100):
    """
    Read average angle data from a JSON file and plot the histogram of angles.
    
    Parameters:
    input_file (str): Path to the input JSON file containing average angle data.
    output_file (str): Path to save the histogram plot.
    """
    
    pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    
    print(legend[1])
    plt.figure(figsize=(10, 6))
    plt.subplot(211)
    plt.hist(np.array(data[coord]["data"]["type"])*57.2598, bins=180, color=pres_colors[color], alpha=0.7, edgecolor='black')
    plt.axvline(data[coord]["avg"]["type"]*57.2598, color='red', linestyle='dashed', linewidth=1)
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Frequency')
    plt.ylim(0, limy)
    plt.xlim(0, 180)
    plt.legend(["avg",legend[0]], loc='upper right')
    plt.title(title)
    
    plt.subplot(212)
    plt.hist(np.array(data[coord]["data"]["atoms"])*57.2598, bins=180, color=pres_colors[color], alpha=0.7, edgecolor='black', hatch='//')
    plt.axvline(data[coord]["avg"]["atoms"]*57.2598, color='red', linestyle='dashed', linewidth=1)
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Frequency')
    plt.ylim(0, limy)
    plt.xlim(0, 180)
    
    plt.title(title)   
    plt.tight_layout()
    
    plt.legend(["avg",legend[1]], loc='upper right')
    plt.savefig(output_file, dpi=150)
    plt.close()
    
def main(): 
    parser = argparse.ArgumentParser(description='Plot histogram of average angles')
    parser.add_argument('--input', type=str, required=True, help='Input JSON file with average angle data')
    parser.add_argument('--output', type=str, required=True, help='Output file for the histogram plot')
    parser.add_argument('--coord', type=str, required=True, help='Coordinate to plot (e.g., "A2_A1_A3")')
    parser.add_argument('--color', type=int, default=0, help='Color index for the histogram (default: 0)')
    parser.add_argument('--title', type=str, default="Average Angles Histogram", help='Title for the histogram plot')
    parser.add_argument('--legend',  nargs='+', default="Average Angles", help='Legend for the histogram plot')
    parser.add_argument('--limy', type=int, default=100, help='Y-axis limit for the histogram (default: 100)')
    
    args = parser.parse_args()

    hist_avg_angle(args.input, args.output, args.coord, args.color,args.title, args.legend, args.limy) 
    
if __name__ == "__main__":
    main()