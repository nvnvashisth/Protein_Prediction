import sys
import re
import string
from collections import Counter
fasta = []
test = []
count = 0
with open("tests.fasta") as file_one:
    for line in file_one:
        count=count+1
        if(count < 20):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    test.append(''.join(fasta))
                    fasta = []
                continue
            fasta.append(line)
        else:
            break
print(count)
if fasta:
    test.append(''.join(fasta))
    test = test[1:]
    removetable = str.maketrans('', '', '*')
    out_list = [s.translate(removetable) for s in test]
#print(out_list)
check = []
for index, seq in enumerate(out_list):
    print("Read %d Sequence"%index)
    if(index < 100000):
        for j in range(0, len(seq)):
            if(j<len(seq)-3):
                check.append(seq[j]+seq[j+1]+seq[j+2]+seq[j+3])

with open("output.txt",'w') as f:
    print('Filename:', Counter(check), file=f)

#print(Counter(check))