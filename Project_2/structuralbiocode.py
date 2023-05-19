import Bio.PDB
import csv

domains = []

#generating a list of lists containing predicted clusters from pae_to_domains...
with open('/home/acube/anuproj/ANUProj/Project_2/pae_to_domains/clusters.csv') as f:
        output = csv.reader(f)
        for line in output:
            domains.append(line)

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
    Bio.PDB.Dice.extract(structure, 'A', int(domain[0]), int(domain[-1]), f'output/{domain[0]}-{domain[-1]}_domain.PDB')