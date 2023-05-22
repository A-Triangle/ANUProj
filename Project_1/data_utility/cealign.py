from matplotlib import pyplot as plt
import Bio.PDB
import argparse
from numpy import array
from numpy.linalg import norm
import os
import shutil

def main():
    
    parser = argparse.ArgumentParser(description='Alignes input structures and calculates structure wide RMSD, by residue RMSD, and by reside delta c alpha for analysis. Outputs aligned model.')

    parser.add_argument('reference_filename', metavar='reference_filename', type=str, help='enter the filename of .fasta target sequence')
    parser.add_argument('sample_filename', metavar='sample_filename', type=str, help='enter the filename of the .txt mut_list file')
    args = parser.parse_args()

    reference_filename = args.reference_filename
    sample_filename = args.sample_filename

    if os.path.isdir('output'):
        shutil.rmtree('output')
        os.mkdir('output'); os.mkdir(f'output/{reference_filename.removesuffix(".pdb")}')
    else: 
        os.mkdir('output'); os.mkdir(f'output/{reference_filename.removesuffix(".pdb")}')
    
    try: 
        
        with open(f'input/{sample_filename}') as handle:
            structure_sample = Bio.PDB.PDBParser().get_structure('Sample',handle)

        with open(f'input/{reference_filename}','r') as handle:
            structure_ref = Bio.PDB.PDBParser().get_structure('Ref',handle)
            
    except NameError:
        exit('Error: filenotfound - Ensure target files are in "input" folder, check spelling.')

    model_ref = structure_ref[0]; model_sample = structure_sample[0]
    
    error, model_ref, model_sample = align_structures(model_ref, model_sample)
    
    parser = Bio.PDB.PDBIO() 
    parser.set_structure(structure_sample) 
    parser.save(f'output/{reference_filename.removesuffix(".pdb")}/align-{reference_filename}-sample')
    
    parser.set_structure(structure_sample) 
    parser.save(f'output/{reference_filename.removesuffix(".pdb")}/align-{reference_filename}-reference')

    residues, CA_comparison, CF_comparison= c_alpha_rmsd(model_ref, model_sample)   

    plt.plot(residues,CA_comparison)
    
    plt.title('C alpha Atomic distances per residue')
    plt.xlabel('Amino acid residue')
    plt.ylabel('alpha carbon distance (Angstroms)')
    
    plt.savefig('output/output_fig')
    
    plt.plot(residues, CF_comparison)
    
    plt.title('Distal atom Atomic distances per residue')
    plt.xlabel('Amino acid residue')
    plt.ylabel('alpha carbon distance (Angstroms)')
    
    plt.savefig('output/output_fig_CF')

def align_structures(model_ref, model_sample):

    struct_align = Bio.PDB.cealign.CEAligner(window_size=1)
    struct_align.set_reference(model_ref)
    struct_align.align(model_sample, transform=True)
    
    return struct_align.rms, model_ref, model_sample

def c_alpha_rmsd(model_ref, model_sample):
    
    sample_data=[]
    for chain_sample in model_sample:
        for resid_sample in chain_sample:
            sample_data.append({'residue_id':resid_sample.id[1],'residue_CA':array(resid_sample['CA'].coord),'residue_CF':array(resid_sample.child_list[-1].coord)})
        
    reference_data=[]    
    for chain_ref in model_ref:
        for resid_ref in chain_ref:
            reference_data.append({'residue_id':resid_ref.id[1],'residue_CA':array(resid_ref['CA'].coord),'residue_CF':array(resid_ref.child_list[-1].coord)})
    
    #aligning residues for comparison
 
    i = 0 
    init_id = 0
    
    residues = []
    CA_deviations = []
    CF_deviations = []
    
    while i < len(reference_data):
        
        magnitude_comparison_CA = 100
        j = i
        
        while j < len(reference_data):
            comparison_CA = (reference_data[i]['residue_CA'] - sample_data[j]['residue_CA']) 
            magnitude_CA=norm(comparison_CA, 2)
            
            if magnitude_comparison_CA > magnitude_CA:
                magnitude_comparison_CA = magnitude_CA
                init_id = j
                
            j += 1 
        
        comparison_CF = (reference_data[init_id]['residue_CF'] - sample_data[init_id]['residue_CF']) 
        magnitude_comparison_CF=norm(comparison_CF, 2)
        
        CA_deviations.append(float(magnitude_comparison_CA))
        CF_deviations.append(float(magnitude_comparison_CF))
        residues.append(int(sample_data[i]['residue_id']))
        
        i += 1
        
    return residues, CA_deviations, CF_deviations


if __name__ == '__main__':
    main()
