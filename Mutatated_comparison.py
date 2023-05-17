import sys
import requests
import os
import shutil
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import MutableSeq

def main():  
    
    try:
        assert len(sys.argv) == 4
    except AssertionError:
        sys.exit('Incorrect number of arguments given')
    
    path = 'output/'+sys.argv[1]
    
    if os.path.exists(path):
        shutil.rmtree(path)
        os.mkdir(path)
    else:
        os.mkdir(path)
        
    query = queryprotein(sys.argv, path)
    interactingproteins = assembleinteracting(query)

    for gene in interactingproteins:
        FindSeq = requests.get(f'https://rest.ensembl.org/sequence/id/{gene["EMBL"]}?content-type=text/x-fasta;species=homo_sapiens;type=cds').text 
    
        with open(path+'/Comp-'+gene["GENE_NAME"]+'-'+sys.argv[1]+'.fasta', 'w') as handle:
            SeqIO.write(query, handle, 'fasta')
            handle.write(FindSeq)

def queryprotein(input, path):
    ID = input[1]; position = input[2]; variant = input[3]

    ProteinMatch = requests.get(f'https://string-db.org/api/json/get_string_ids?identifiers={ID}').json()
    embl_id=ProteinMatch[0]["stringId"].removeprefix("9606.")

    InitialSeq=requests.get(f'https://rest.ensembl.org/sequence/id/{embl_id}?content-type=text/x-fasta;species=homo_sapiens;type=cds').text
    
    with open(f'{path}/{ID}.fasta', 'w') as handle:
        handle.write(InitialSeq)

    with open(f'{path}/{ID}.fasta', 'r') as handle:
        query = SeqIO.read(handle,'fasta')

    query.id = embl_id
    aminoacids = MutableSeq(query.seq)
    aminoacids[int(position)-1] = variant.upper()
    query.seq = aminoacids

    with open(f'{path}/{ID}.fasta', 'w') as handle:
        SeqIO.write(query, handle, 'fasta')
        
    return query

def assembleinteracting(query):
    
    interacting = requests.get(f'https://string-db.org/api/json/network?identifiers={query.id}').json()

    output=[]
    for entry in interacting:
        if entry['stringId_A'] not in output:
            output.append({'EMBL':entry['stringId_A'].removeprefix("9606."),'GENE_NAME':entry['preferredName_A']})

    return output

if __name__ == '__main__':
    main()