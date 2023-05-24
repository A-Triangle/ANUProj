from matplotlib import pyplot as plt
import Bio.PDB
import argparse
from numpy import array
from numpy.linalg import norm
import os
import shutil
import json

def main():
    
    #Providing CLI
    
    parser = argparse.ArgumentParser(description='Alignes input structures and calculates structure wide RMSD, by residue RMSD, and by reside delta c alpha for analysis. Outputs aligned model.')

    parser.add_argument('reference_filename', metavar='reference_filename', type=str, help='enter the filename of .fasta target sequence')
    parser.add_argument('sample_filename', metavar='sample_filename', type=str, help='enter the filename of the .txt mut_list file')
    args = parser.parse_args()

    reference_filename = args.reference_filename
    sample_filename = args.sample_filename

    #Organising filesystem

    if os.path.isdir('output'):
        shutil.rmtree('output')
        os.mkdir('output'); os.mkdir(f'output/{reference_filename.removesuffix(".pdb")}')
    else: 
        os.mkdir('output'); os.mkdir(f'output/{reference_filename.removesuffix(".pdb")}')
        
    #Attemtping to open structure files
    
    try: 
        
        with open(f'input/{sample_filename}', 'r') as handle:
            structure_sample = Bio.PDB.PDBParser().get_structure('Sample',handle)

        with open(f'input/{reference_filename}','r') as handle:
            structure_ref = Bio.PDB.PDBParser().get_structure('Ref',handle)
            
    except NameError:
        exit('Error: filenotfound - Ensure target files are in "input" folder, check spelling.')

    model_ref = structure_ref[0]; model_sample = structure_sample[0]
    
    #aligining structures, calculating outputs and generating graphs
    
    error, model_ref, model_sample = align_structures(model_ref, model_sample)
    
    output = c_alpha_rmsd(model_ref, reference_filename, model_sample, sample_filename)   
    
    make_graphs(output, error, sample_filename, reference_filename)
    
    #Writing output files 
    
    parser = Bio.PDB.PDBIO() 
    
    parser.set_structure(structure_ref) 
    with open (f'output/{reference_filename.removesuffix(".pdb")}/align-{reference_filename}', 'a+') as handle: 
        parser.save(handle)
    
    parser.set_structure(structure_sample) 
    with open (f'output/{reference_filename.removesuffix(".pdb")}/align-{sample_filename}','a+') as handle: 
        parser.save(handle)
    
    with open(f'output/{reference_filename.removesuffix(".pdb")}-{sample_filename.removesuffix(".pdb")}.json', 'w+') as file:
        file.write(json.dumps(output, indent = 2))


def align_structures(model_ref, model_sample):
    
    #using CEalgn to transform structures to best alignment assumes that the first model of each pdb file is the target of comparison

    struct_align = Bio.PDB.cealign.CEAligner(window_size=1)
    struct_align.set_reference(model_ref)
    struct_align.align(model_sample, transform=True)
    
    return struct_align.rms, model_ref, model_sample

def c_alpha_rmsd(model_ref, reference_filename, model_sample, sample_filename):
    
    #calling dssp to extract SASA and secondary structure
    
    dssp_sample = Bio.PDB.DSSP(model_sample, f'input/{sample_filename}', dssp='/home/lw/anaconda3/envs/biochem/bin/mkdssp')
     
    #looping over all of the sample structure residues 
      
    sample_data=[]
    for chain_sample in model_sample:
        for resid_sample in chain_sample:
            
            res_data = dssp_sample[chain_sample.id, resid_sample.id[1]]
            
            #branching aa's in pdb files notated by number reffering to the branch. 
            #reversing the the residues child list and excluding numbered atom id's gives
            #the most distal unambigous atom for comparison.
            
            for atom in reversed(resid_sample.child_list):
                if atom.id.isalpha():
                    CF=atom.coord 
                    break
                break
    
            #extracting relevant data
            
            residue =({
            'residue_id':resid_sample.id[1],
            'residue_CA':array(resid_sample['CA'].coord),
            'residue_CF':array(CF),
            'SASA':res_data[3]})
            
            sample_data.append(residue)
            
    #same as above but for the reference structure
        
    dssp_ref = Bio.PDB.DSSP(model_ref, f'input/{reference_filename}', dssp='/home/lw/anaconda3/envs/biochem/bin/mkdssp')  
    
    reference_data=[]    
    for chain_ref in model_ref:
        for resid_ref in chain_ref:
            
            res_data = dssp_ref[chain_ref.id, resid_ref.id[1]]
            
            for atom in reversed(resid_ref.child_list):
                if atom.id.isalpha():
                    CF=atom.coord 
                    break
                break
            
            #editing secondary structure assignments to only include beta sheets, alpha helices and loops
            
            if res_data[2] in ['-','T','S','B']:
                structure = 'L'
            elif res_data[2] in ['E','B'] :
                structure = 'E'
            else:
                structure = res_data[2]
                
            residue =({
            'residue_id':resid_ref.id[1],
            'residue_CA':array(resid_ref['CA'].coord),
            'residue_CF':array(CF),
            'secondary_structure':structure,
            'SASA':res_data[3]})
            
            reference_data.append(residue)
            
    #generating empty output datastructure to populate, this method makes matplotlib geneartion simple
    
    output = {
    'H':[[],[],[],[],[],[]],
    'G':[[],[],[],[],[],[]],
    'I':[[],[],[],[],[],[]],
    'E':[[],[],[],[],[],[]],
    'L':[[],[],[],[],[],[]]
    }
    
    i = 0
    
    while i < len(reference_data):
        
        #arbitrary initation constant for magnitude comparison
        
        magnitude_comparison_CA = 100
        j = i
        
        
        while j < (j + 10) :
            
            #extracting l2 norms from the difference of 3x1 atom coordinate matrices
            #finds distance between any two CA atoms
            
            comparison_CA = (reference_data[i]['residue_CA'] - sample_data[j]['residue_CA']) 
            magnitude_CA=norm(comparison_CA, 2)
            
            #finds the closest CA atom in opposing structure for each residue to allow for 
            #pairwise comparisons in the case of small indels - the inner loop position is
            #extracted via init_id to allow for singleton calculations for remaining values
            
            if magnitude_comparison_CA > magnitude_CA:
                magnitude_comparison_CA = magnitude_CA
                init_id = j
                
            j += 1 
            
        print(i, i-init_id)    
            
        #calculating distance of furthest unambigous atom and the delta of the solvant accesable surface area
        
        comparison_CF = (reference_data[i]['residue_CF'] - sample_data[init_id]['residue_CF']) 
        magnitude_comparison_CF=norm(comparison_CF, 2)
        
        delta_SASA = reference_data[i]['SASA'] - sample_data[init_id]['SASA']
        
        #populating emtpy datapoints with the residue id and null to allow for plotting of secondary structure
        #with matplotlib
        
        types = ['H','G','I','E','L']
        types.remove(reference_data[i]['secondary_structure'])
        
        for type in types:
            output[type][0].append(int(sample_data[i]['residue_id']))
            output[type][1].append(None)
            output[type][2].append(None)
            output[type][3].append(None)   
            output[type][4].append(None)
            output[type][5].append(None)
        #populating genuine data
        
        output[reference_data[i]['secondary_structure']][0].append(int(sample_data[i]['residue_id']))
        output[reference_data[i]['secondary_structure']][1].append(float(magnitude_comparison_CA))
        output[reference_data[i]['secondary_structure']][2].append(float(magnitude_comparison_CF))
        output[reference_data[i]['secondary_structure']][3].append(float(delta_SASA))
        output[reference_data[i]['secondary_structure']][4].append(float(reference_data[i]['SASA']))
        output[reference_data[i]['secondary_structure']][5].append(float(sample_data[i]['SASA']))
        
        i += 1
        
    return output

def make_graphs(output, error, sample_filename, reference_filename):
    
    #making output graphs
    
    fig, ax = plt.subplots(dpi=400)
    
    ax.plot(output['H'][0], output['H'][1], label='Helix')
    ax.plot(output['E'][0], output['E'][1], label='Beta-Sheet')
    ax.plot(output['L'][0], output['L'][1], label='Loop')
    ax.legend()
   
    plt.savefig(f'output/{reference_filename.removesuffix(".pdb")}-{sample_filename.removesuffix(".pdb")}_CA')
    
    fig2, ax = plt.subplots(dpi=400)

    ax.plot(output['H'][0], output['H'][2], label='Helix')
    ax.plot(output['E'][0], output['E'][2], label='Beta-Sheet')
    ax.plot(output['L'][0], output['L'][2], label='Loop')
    ax.legend()
   
    plt.savefig(f'output/{reference_filename.removesuffix(".pdb")}-{sample_filename.removesuffix(".pdb")}_CF')
    
    fig3, ax = plt.subplots(dpi=400)
    
    ax.plot(output['H'][0], output['H'][3], label='Helix')
    ax.plot(output['E'][0], output['E'][3], label='Beta-Sheet')
    ax.plot(output['L'][0], output['L'][3], label='Loop')
    ax.legend()
   
    plt.savefig(f'output/{reference_filename.removesuffix(".pdb")}-{sample_filename.removesuffix(".pdb")}_delta-sasa')
    
    fig4, ax = plt.subplots(dpi=400)
    
    ax.plot(output['H'][0], output['H'][4], label='Helix')
    ax.plot(output['E'][0], output['E'][4], label='Beta-Sheet')
    ax.plot(output['L'][0], output['L'][4], label='Loop')
    ax.legend()
    
    plt.savefig(f'output/{reference_filename.removesuffix(".pdb")}-{sample_filename.removesuffix(".pdb")}_reference-sasa')
    
    fig5, ax = plt.subplots(dpi=400)

    ax.plot(output['H'][0], output['H'][5], label='Helix')
    ax.plot(output['E'][0], output['E'][5], label='Beta-Sheet')
    ax.plot(output['L'][0], output['L'][5], label='Loop')
    ax.legend()
   
    plt.savefig(f'output/{reference_filename.removesuffix(".pdb")}-{sample_filename.removesuffix(".pdb")}_sample-asa')

if __name__ == '__main__':
    main()
