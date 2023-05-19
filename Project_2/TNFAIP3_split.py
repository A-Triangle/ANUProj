from Bio.PDB import *
import pymol
import requests
import sys

def main():
    pAE_data = downloadmodeldata(sys.argv[1])
#    folded_residues = unpackalignederror(pAE_data, sys.argv[2])
    
    #using biopython to copy structure data to new files containing regions of defined fold
    
def downloadmodeldata(PDB):
    #findinguniprotacc.for given pdbfile
    
    requestURL = f'https://www.ebi.ac.uk/proteins/api/proteins/PDB:{PDB}?'
    data = requests.get(requestURL, headers={"Accept" : "application/json"}).text
    
    print(data)
   
    
    #usinguniprogaccesion to download .pdb and pAE files
    
    #parsing pAE as json object and returning relevant data, closing .pdb file
    
#def unpackalignederror(pAE, cutoff):
    #somehow seelct regions with low predicted aligned error (tailorable cutoff)
    #return a list of ranges of amino acids with defined structures
    
if __name__ == '__main__':
    main()
