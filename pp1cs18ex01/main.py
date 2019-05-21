# The following method should take a strand as an argument and
# return a complementary strand

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def complementary(strand):
     return ''.join([complement[character] for character in strand.upper()])
