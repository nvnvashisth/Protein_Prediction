import string

aa = {
    "TTT" : "F", "TTC" : "F",
    "TTA" : "L", "TTG" : "L", "CTT" : "L", "CTC" : "L", "CTA" : "L", "CTG" : "L",
    "ATT" : "I", "ATC" : "I", "ATA" : "I",
    "ATG" : "M",
    "GTT" : "V", "GTC" : "V", "GTA" : "V", "GTG" : "V",
    "TCT" : "S", "TCC" : "S", "TCA" : "S", "TCG" : "S",
    "CCT" : "P", "CCC" : "P", "CCA" : "P", "CCG" : "P",
    "ACT" : "T", "ACC" : "T", "ACA" : "T", "ACG" : "T",
    "GCT" : "A", "GCC" : "A", "GCA" : "A", "GCG" : "A",
    "TAT" : "Y", "TAC" : "Y",
    "TAA" : "Ochre", "TAG" : "Amber",
    "CAT" : "H", "CAC" : "H",
    "CAA" : "Q", "CAG" : "Q",
    "AAT" : "N", "AAC" : "N",
    "AAA" : "K", "AAG" : "K",
    "GAT" : "D", "GAC" : "D",
    "GAA" : "E", "GAG" : "E",
    "TGT" : "C", "TGC" : "C",
    "TGA" : "Opal",
    "TGG" : "W",
    "CGT" : "R", "CGC" : "R", "CGA" : "R", "CGG" : "R",
    "AGT" : "S", "AGC" : "S",
    "AGA" : "R", "AGG" : "R",
    "GGT" : "G", "GGC" : "G", "GGA" : "G", "GGG" : "G"
}

str1 = "ttacgcagtaccgccaacagtcaggttatcaactttcagcgttggctggcccacgccaaccggcaaactttgcccttctttaccgcagacacccacgccgttatccagtttcaggtcgttgccaaccatcgaaatctgctgcatggtttcgataccggaaccaatcaacgttgcgcctttcaccggcttcgttactttaccgttttcaatcagatatgcttctgaagtggagaaaacgaatttgccggaggtgatatccacctgaccgccaccaaagttcggtgcatagataccgtactcaacggattcaataatttcctgcggggtcgatttacccggcagcatataggtgttggtcatacgcggcatgggcagatgggcgtaggattcacggcgaccgttgccagtcggcgtcatccccatcaaacgcgcgttgagtttatcctgcatgtagcctttcagaatgccgttctcaatcagcacgttgtactggcctggcgtaccttcgtcatcaatcgccaccgaaccacggcgatcgaccatcgtaccgtcatcaaccacggtacacagttctgaagccaccagctccccgacctgtccactaaatactgaagtgccacggcggttgaagtcgccttccagaccgtgacctaccgcttcatgcaacagcacgcccggccaacctgcgccaagtactaccggcatggtgcccgctggtgcagcaacggcagaaagattgaccagcgccatacgtactgcttcttttgcccatgcatccgcacggacttcgccgtcgagatcggcaaggaagaactcataaccaaaacgaccgccgccgccactggcaccgcgttcgcgtttgccatcttcttcgatcagaacgctcacggaaagacgcaccagcggacggacatccgccgctagcgtgccgtcggtggccgcaaccagaattaattcatagacaccactgaggctggcagtcacttcctgtacgcgcttgtccgcttcgcgggcaaccttatcgacgcgacgcaggatatccagcttctcttcacggctcatgctttgcagcggatctaccgaggtatacaacgggctatgctctaccgcgcccagcgtctgtactttgccatcaccgctatcacgggcgatggtgcgcgccgcttgcgcactctgttccagcgccagcaggctgatttggtcagcgtaagcaaatccggttttttcaccgctgattgcacgcacaccaacgccctgatcgatgttgtaagaaccatctttaataatgcggtcttctaaaacccaggattcgtgatagctcgactgaaaatagagatcgccataatcaaggcgacgttcggccagttgaccgaggatcacgaacaagtcctgatgtttcaggccgttcgccgctagcaattgttcacttaccaggttaagactcatcg"
#str1 = ""
orfs = []

def translate(str1):
    str1 = str1.upper()
    l = len(str1)
    for i in str1:
        if i != 'A' and i != 'C' and i != 'G' and i != 'T':
            raise TypeError('Given sequence is not a DNA sequence')
    
    s1 = 0
    ctr = 0
    for i in range(0, 3):
        str2 = ""
        for j in range(i, l-2, 3):
            str3 = ""
            str3 = str1[j] + str1[j + 1] + str1[j + 2]
            str3 = aa[str3]
            if str3 == "Ochre" or str3 == "Amber" or str3 == "Opal":
                ctr = 0
                s2 = j + 2
                if len(str2) > 30:
                    orfs.append(tuple((s1, s2, str2)))
                str2 = ""
                j = j + 1
            if ctr == 1:
                str2 = str2 + str3
            if str3 == "M" and ctr == 0:
                s1 = j
                print (j)
                ctr = 1
                str2 = str3

    str1 = str1[::-1]
    str1 = complementary(str1)

    for i in range(0, 3):
        str2 = ""
        for j in range(i, l-2, 3):
            str3 = ""
            str3 = str1[j] + str1[j + 1] + str1[j + 2]
            str3 = aa[str3]
            if str3 == "Ochre" or str3 == "Amber" or str3 == "Opal":
                ctr = 0
                s2 = j + 2
                if len(str2) > 30:
                    s1 = len(str1) - s1 - 1
                    s2 = len(str1) - s2 - 1
                    orfs.append(tuple((s1, s2, str2)))
                str2 = ""
                j = j + 1
            if ctr == 1:
                str2 = str2 + str3
            if str3 == "M" and ctr == 0:
                s1 = j
                ctr = 1
                str2 = str3 
            
    return orfs

def complementary(strand):
    
    strand2 = ""
    
    for i in range(0, len(strand)):
        if strand[i] == "T":
            strand2 += "A"
        elif strand[i] == "A":
            strand2 += "T"
        elif strand[i] == "C":
            strand2 += "G"
        elif strand[i] == "G":
            strand2 += "C"

    return strand2



orfs2 = translate(str1)
print(orfs2)