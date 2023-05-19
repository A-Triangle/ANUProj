import Bio.PDB
import csv
import re 

def main():
    select = input('specify to generate domains from predicted templates (0) or from crystalographic fragments (1)')

    if select == 0:
        get_predicted_clusters()
        get_predicted_fragments()
    elif select == 1:
        get_experimental_fragments()
        

def get_predicted_clusters():
     
    #generating a list of lists containing predicted clusters from pae_to_domains...
    domains = []
    with open('/home/acube/anuproj/ANUProj/Project_2/pae_to_domains/clusters.csv') as f:
        output = csv.reader(f)
        for line in output:
            domains.append(line)
    return output


def get_predicted_domains():

    #opening structure file to split
    with open('AF-P21580-F1-model_v4.pdb') as f:
        structure = Bio.PDB.PDBParser().get_structure('a20',f)

    #filtering clusers.csv file
    for domain in domains:
        i = 0
        while i < len(domain):
            if domain[i] == '':
                del domain[i:len(domain)]
                break
            i+=1
    #generating domain pdb files
        Bio.PDB.Dice.extract(structure, 'A', int(domain[0]), int(domain[-1]), f'output/{domain[0]}-{domain[-1]}_domain.PDB')

def get_experimental_fragments(): 

if __name__ == '__main__':
     main()