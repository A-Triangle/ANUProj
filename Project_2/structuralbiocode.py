import Bio.PDB
import csv

domains = []

#generating a list of lists containing predicted clusters from pae_to_domains...
with open('/home/acube/anuproj/ANUProj/Project_2/tristanic/pae_to_domains/clusters.csv') as f:
        output = csv.reader(f)
        for line in output:
            domains.append(line)

with open('AF-P21580-F1-model_v4.pdb') as f:
    structure = Bio.PDB.PDBParser().get_structure('a20',f)

model = structure[0]
for chain in model:
    print(f'chain {chain}, chainid: {chain.id}')

chain_A = model['A']
for res in chain_A:
    print(f'residue name {res.resname}, resid: {res.id[2]}')

    