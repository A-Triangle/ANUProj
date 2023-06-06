import json 
import matplotlib as plt
import os 
import argparse

def main():
    
    #Providing CLI
    
    
    parser = argparse.ArgumentParser(description='Alignes input structures and calculates structure wide RMSD, by residue RMSD, and by reside delta c alpha for analysis. Outputs aligned model.')

    parser.add_argument('', metavar='reference_filename', type=str, help='enter the filename of .fasta target sequence')
    parser.add_argument('sample_filename', metavar='sample_filename', type=str, help='enter the filename of the .txt mut_list file')
    args = parser.parse_args()

    reference_filename = args.reference_filename.removesuffix(".pdb")
    sample_filename = args.sample_filename.removesuffix(".pdb")


    #Organising filesystem
    

    if os.path.isdir(f'output/{sample_filename}'):
        shutil.rmtree(f'output/{sample_filename}')
        os.mkdir(f'output/{sample_filename}')
    else: 
        os.mkdir(f'output/{sample_filename}')
        
        
    #Attemtping to open structure files
    
    
    try: 
        
        with open(f'input/{sample_filename}.pdb', 'r') as handle:
            structure_sample = Bio.PDB.PDBParser().get_structure('Sample',handle)

        with open(f'input/{reference_filename}.pdb','r') as handle:
            structure_ref = Bio.PDB.PDBParser().get_structure('Ref',handle)
            
    except NameError:
        exit('Error: filenotfound - Ensure target files are in "input" folder, check spelling.')

    model_ref = structure_ref[0]; model_sample = structure_sample[0]
    

if __name__ == '__main__':
    main()