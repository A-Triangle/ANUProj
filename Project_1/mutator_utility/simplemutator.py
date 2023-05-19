from Bio import SeqIO                                             
from Bio.Seq import MutableSeq
import sys

with open(f'input_data/{sys.argv[1]}', 'r') as handle:
        query = SeqIO.read(handle,'fasta')

aminoacids = MutableSeq(query.seq)
aminoacids[int(sys.argv[2])-1] = sys.argv[3].upper()
query.seq = aminoacids

if sys.argv[4].lower() == 'yes':
    data = input('provide aa range for intended domain (#-#)')
    start, end = data.split('-')

    start = int(start); end = int(end)
    domain = query.seq[start:end]
    query.seq = domain 

    query.id = query.id + str(start) + '-' + str(end) 

    print(query.seq)

    with open(f'truncated/{sys.argv[1]}-{sys.argv[2]+sys.argv[3]}-{data}', 'w') as handle:
        SeqIO.write(query, handle, 'fasta')

elif sys.argv[4].lower() == 'no':

    print(query.seq)
    with open(f'mutated/{sys.argv[1]}-{sys.argv[2]+sys.argv[3]}', 'w') as handle:
        SeqIO.write(query, handle, 'fasta')

else:
     sys.exit('invalid command')