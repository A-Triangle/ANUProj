import matplotlib.pyplot as plt
import argparse
import numpy
import json

def main():
    
    #Providing CLI
    
    
    parser = argparse.ArgumentParser(description='Alignes input structures and calculates structure wide RMSD, by residue RMSD, and by reside delta c alpha for analysis. Outputs aligned model.')

    parser.add_argument('aligned_error_filename', metavar='aligned_error_filename', type=str, help='enter the filename of the predicted aligned error json')
    parser.add_argument('-xresidues', type=str, help='enter the x range to examine, sperated by a comma')
    parser.add_argument('-yresidues', type=str, help='enter the y range to examine, sperated by a comma')
    args = parser.parse_args()

    filename = args.aligned_error_filename.removesuffix(".json")
    xmin, xmax = args.xresidues.split(',')
    ymin, ymax = args.yresidues.split(',')
    
    xmin = int(xmin); xmax = int(xmax); ymin = int(ymin); ymax = int(ymax)
    print(xmin,xmax,ymin,ymax)
  
    #Attemtping to open structure files
    
    
    try: 
        
        with open(f'input/{filename}.json','r') as handle:
            data = json.load(handle)
            
    except NameError:
        exit('Error: filenotfound - Ensure target files are in "input" folder, check spelling.')
        
    #cropping numpy-array
    

    matrix = numpy.array(data['predicted_aligned_error'], dtype=numpy.float64)
    matrix = matrix[xmin:xmax, ymin:ymax]
    
    #producing heat-plot
    
    # Show all ticks and label them with the respective list entries
    
    ax.set_xticks(np.arange(len(farmers)), labels=farmers)
    ax.set_yticks(np.arange(len(vegetables)), labels=vegetables)

    # Rotate the tick labels and set their alignment.
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")
    
    fig, ax = plt.subplots()
    ax.imshow(matrix)

    ax.set_title(f"Predicted aligned error in range {xmin}, {xmax} - {ymin}, {ymax}")
    fig.tight_layout()
    
    plt.savefig(f'output/PAE/{filename}')

if __name__ == '__main__':
    main()