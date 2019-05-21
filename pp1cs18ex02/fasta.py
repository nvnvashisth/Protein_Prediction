import sys
import re
import string
fasta = []
test = []
with open("tests.fasta") as file_one:
    for line in file_one:
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

if fasta:
    test.append(''.join(fasta))
    test = test[1:]
    removetable = str.maketrans('', '', '*')
    out_list = [s.translate(removetable) for s in test]
print(out_list)
'''
seq = []
# Print the test list
for i, row in enumerate(test):
    #print(i)
    seq.append(row)

print(seq)'''