import sys
import requests
import os

#finding string preferred identifier for query proten
ProteinMatch = requests.get(f'https://string-db.org/api/json/get_string_ids?identifiers={sys.argv[1]}').json()
StringID=ProteinMatch[0]["stringId"]

#finding query protein fasta
InitialSeq=requests.get(f'https://rest.ensembl.org/sequence/id/{StringID.removeprefix("9606.")}?content-type=text/x-fasta;species=homo_sapiens;type=cds').text 

#finding interacting proteins
interacting = requests.get(f'https://string-db.org/api/json/network?identifiers={StringID}').json()

#scrubbing output of repeat proteins
output=[]
for entry in interacting:
    if entry['stringId_A'] not in output:
        output.append({'EMBL':entry['stringId_A'],'GENE_NAME':entry['preferredName_A']})

#making new directory for output files
if not os.path.exists('output/'+sys.argv[1]):
    os.mkdir('output/'+sys.argv[1])
        
#finding interacting sequences and assembling comparison files      
for entry in output:
    entry['EMBL'] = entry['EMBL'].removeprefix('9606.')
    FindSeq = requests.get(f'https://rest.ensembl.org/sequence/id/{entry}?content-type=text/x-fasta;species=homo_sapiens;type=cds').text 
    
    file = open('output/'+sys.argv[1]+'/'+entry["GENE_NAME"]+'.fasta', 'w')
    
    file.write(InitialSeq)
    file.write(FindSeq)
    
    file.close()