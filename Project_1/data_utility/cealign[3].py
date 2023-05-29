from matplotlib import pyplot as plt
import Bio.PDB
import Bio.SeqUtils
from Bio import Align
import argparse
from numpy import array
from numpy.linalg import norm
import os
import shutil
import math
import json

def main():
    
    #Providing CLI
    
    
    parser = argparse.ArgumentParser(description='Alignes input structures and calculates structure wide RMSD, by residue RMSD, and by reside delta c alpha for analysis. Outputs aligned model.')

    parser.add_argument('reference_filename', metavar='reference_filename', type=str, help='enter the filename of .fasta target sequence')
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
    
    
    #aligining structures, calculating outputs and generating graphs
    
    
    error, model_ref, model_sample = align_structures(model_ref, model_sample)
    
    sample_parameters = get_data(model_sample, sample_filename)
    reference_parameters = get_data(model_ref, reference_filename)
    
    output = c_alpha_rmsd(reference_parameters, sample_parameters)  
    
    data_processing(output, error, sample_filename, reference_filename) 
    
    
   #Writing output files 
   
    
    parser = Bio.PDB.PDBIO() 
    
    parser.set_structure(structure_ref) 
    with open (f'output/{sample_filename}/align-{reference_filename}.pdb', 'a+') as handle: 
        parser.save(handle)
    
    parser.set_structure(structure_sample) 
    with open (f'output/{sample_filename}/align-{sample_filename}.pdb','a+') as handle: 
        parser.save(handle)
    
    with open(f'output/{sample_filename}/{reference_filename}-{sample_filename}.json', 'w+') as file:
        file.write(json.dumps(output, indent = 2))
        

def align_structures(model_ref, model_sample):
    
    #using CEalgn to transform structures to best alignment, assumes that the first model of each pdb file is the target of comparison

    struct_align = Bio.PDB.cealign.CEAligner(window_size=1)
    struct_align.set_reference(model_ref)
    struct_align.align(model_sample, transform=True)
    
    return struct_align.rms, model_ref, model_sample

def c_alpha_rmsd(reference_parameters, sample_parameters):
    
    reference_sequence, reference_data = reference_parameters; sample_sequence, sample_data = sample_parameters
    
    #generating empty output datastructure 
    
    output = {
    'residue':[],
    'CA':[],
    'CF':[],
    'DSASA':[],
    'Sec_structure':[]
    }

    #generating pariwise gloabl alignment to ensure correct amino acid comparisons
    #(important where seq. length differes or in comparison of xray with potentially 
    # missing disordered segments)
    
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(reference_sequence, sample_sequence)

    reference_alignment = alignments[0][0,:]
    sample_alignment = alignments[0][1,:]
    
    h = j = i = total_magnitude_CA = total_magnitude_CF = total_magnitude_dSASA = comparisons = 0
 
    while h < len(sample_alignment):

        # this method isnt super readable so the explination: '-' represents an alignment gap
        # we index up j (for reference) or i(for sample) in these cases to ensure that pairwise alignments occour,
        # because pdb files simply omit the missing residues (as we are indexing by list number not residue number
        # the loop needs to skip missing segments). H index's wrt the global alignment. Note that this method omits
        # a direct comparison of the mutated segments implicitly.
        
        if reference_alignment[h] == '-':
            j += 1; h += 1
            
        elif sample_alignment[h] == '-':
            i += 1; h += 1
            
        # populating output data
    
        else:
            
            comparisons = comparisons + 1
            
            comparison_CA = (reference_data[i]['residue_CA'] - sample_data[j]['residue_CA']) 
            magnitude_comparison_CA=norm(comparison_CA, 2)
            total_magnitude_CA = total_magnitude_CA + magnitude_comparison_CA
            output['CA'].append(magnitude_comparison_CA)
        
            comparison_CF = (reference_data[i] ['residue_CF'] - sample_data[j]['residue_CF']) 
            magnitude_comparison_CF=norm(comparison_CF, 2)
            total_magnitude_CF = total_magnitude_CF + magnitude_comparison_CF
            output['CF'].append(magnitude_comparison_CF)
        
            delta_SASA = reference_data[i]['SASA'] - sample_data[j]['SASA']
            total_magnitude_dSASA = total_magnitude_dSASA + delta_SASA
            output['DSASA'].append(delta_SASA)
            
            output['Sec_structure'].append(reference_data[i]['secondary_structure'])
            
            output['residue'].append(reference_data[i]['residue_id'])
            
            i += 1; j += 1; h += 1
        
        mean_CA = total_magnitude_CA/comparisons
        mean_CF = total_magnitude_CF/comparisons
        mean_SASA = total_magnitude_dSASA/comparisons
        
    return output, mean_CA, mean_CF, mean_SASA

def get_data(model, filename):
    
    sequence = ''
    data = []
    
    #calling dssp to extract SASA and secondary structure
    
    dssp = Bio.PDB.DSSP(model, f'input/{filename}.pdb', dssp='mkdssp')
     
    #looping over all of the sample structure residues 
      
    for chain in model:
        for resid in chain:
            
            res_data = dssp[chain.id, resid.id[1]]
            
            #branching aa's in pdb files notated by number reffering to the branch. 
            #reversing the the residues child list and excluding numbered atom id's gives
            #the most distal unambigous atom for comparison.
            
            for atom in resid.child_list:
                if atom.id.isalpha():
                    CF=atom.coord 

            #setting residue counting unabigously to 0, building protein sequence.

            numbered = resid.id[1] - model.child_list[0].child_list[0].id[1]
            sequence = sequence + Bio.SeqUtils.IUPACData.protein_letters_3to1[resid.get_resname().capitalize()]
            
            #filtering dssp secondary structure where pertinant
                        
            if res_data[2] in ['-','T','S']:
                structure = 'L'
            elif res_data[2] in ['E','B'] :
                structure = 'E'
            else:
                structure = res_data[2]
            
            residue =({
            'residue_id':numbered,
            'residue_CA':array(resid['CA'].coord),
            'residue_CF':array(CF),
            'secondary_structure': structure,
            'SASA':res_data[3]
            })
            
            data.append(residue)
            
        return (sequence,data)
            

def data_processing(output, error, sample_filename, reference_filename):

    #doing data analysis to identify outlier residues for analysis.
    
    data, mean_CA, mean_CF, mean_SASA = output

    #making output graphs
    
    for datapoint in ['CA','CF','DSASA']:

        mean = 'mean_'+datapoint

        sum=0
        for value in data[datapoint]:
            sum = sum + (value - mean)**2

        standard_deviation = math.sqrt(sum/len(data[datapoint]))
        print(standard_deviation)

        Outlier_Y = []
        Outlier_X = []
        for i, value in enumerate(data[datapoint]):
            if value > 3*standard_deviation:
                Outlier_Y.append(value)
                Outlier_X.append(i)
        
        fig, ax = plt.subplots(dpi=400)
        ax.plot(data['residue'], data[datapoint])

        for outlier in Outlier_Y:
            plt.annotate()

        #funky way of including a correct legend
        
        plt.axvspan(0,0, color ='red', alpha=0.1, label = 'Loop')
        plt.axvspan(0,0, color ='blue', alpha=0.1, label = 'Beta-Sheet')
        plt.axvspan(0,0, color ='green', alpha=0.1, label = 'Alpha-Helix')
        plt.axvspan(0,0, color ='yellow', alpha=0.1, label = '310-Helix')
        plt.axvspan(0,0, color ='pink', alpha=0.1, label = 'Pi-Helix')
        
        for i, structure in enumerate(data['Sec_structure']):
            if structure == 'L':
                plt.axvspan(i,i+1, color ='red', alpha=0.1)
            elif structure == 'E':
                plt.axvspan(i,i+1, color ='blue', alpha=0.1)
            elif structure == 'H':
                plt.axvspan(i,i+1, color ='green', alpha=0.1)
            elif structure == 'G':
                plt.axvspan(i,i+1, color ='yellow', alpha=0.1)
            elif structure == 'I':
                plt.axvspan(i,i+1, color ='pink', alpha=0.1)

        plt.title(sample_filename + ' ' + reference_filename + ' ' + datapoint)
        plt.xlabel('Residue Number')
        plt.ylabel('Distance(Angstroms)')

        plt.legend()

        plt.savefig(f'output/{sample_filename}/{reference_filename}-{sample_filename}_{datapoint}')

if __name__ == '__main__':
    main()