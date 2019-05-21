##############
# Exercise 2.6
##############

# You can use tests.fasta for your own testing. Good luck!
from collections import Counter
import sys
import re
import string



def get_counts():
    
    with open('test.txt') as f:
        aa_seq = f.read().splitlines()
    count = len(aa_seq)
    return count
          
#get_counts()

def get_average_length():
    sum = 0
    with open('test.txt') as f:
        aa_seq = f.read().splitlines()
    count = len(aa_seq)
    for line in aa_seq:
        sum = sum + len(line)
    return sum/count

#get_average_length()

def read_fasta(filename):
    fasta = []
    test = []
    head = []    
    with open(filename) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                head.append(active_sequence_name)
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
    file = open("test.txt","w")
    for item in out_list:
        file.write("%s\n" % item)
    file.close()
    #print(out_list)
    return (head,out_list)

#read_fasta("tests.fasta")

def get_abs_frequencies():
    #aa_seq = read_fasta("tests.fasta")
    total = []
    with open('test.txt') as f:
        aa_seq = f.read().splitlines()
    for line in aa_seq:
        for char in line:
            total.append(char)
    return dict((x,total.count(x)) for x in set(total))

#get_abs_frequencies()



def get_av_frequencies():
    
    total = []
    with open('test.txt') as f:
        aa_seq = f.read().splitlines()
    for line in aa_seq:
        for char in line:
            total.append(char)

    a = dict((x,total.count(x)) for x in set(total))
    sum = 0
    b = {}
    for key, value in a.items():
        sum = sum + value
    for key, value in a.items():
        value = value / sum
        b.setdefault(key)
        b[key] = value 
    return b

#get_av_frequencies()