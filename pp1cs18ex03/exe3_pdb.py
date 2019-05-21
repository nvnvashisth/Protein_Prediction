# -*- coding: utf-8 -*-
"""
IMPORTANT!:
Before writing an email asking questions such as
'What does this input has to be like?' or 
'What return value do you expect?' PLEASE read our
exercise sheet and the information in this template
carefully.
If something is still unclear, PLEASE talk to your
colleagues before writing an email!

If you experience technical issues or if you find a
bug we are happy to answer your questions. However,
in order to provide quick help in such cases we need 
to avoid unnecessary emails such as the examples
shown above.
"""

from Bio.PDB.MMCIFParser import MMCIFParser  # Tip: This module might be useful for parsing...

from Bio.PDB.Polypeptide import PPBuilder


############# Exercise 2: Protein Data Bank #############
# General remark: In our exercise every structure will have EXACTLY ONE model.
# This is true for nearly all X-Ray structures. NMR structures have several models.
class PDB_Parser:
    parser = MMCIFParser()  # parser object for reading in structure in CIF format
    ppb = PPBuilder()

    def __init__(self, path):
        '''
            Initialize every PDB_Parser with a path to a structure-file in CIF format.
            An example file is included in the repository (7ahl.cif).
            Tip: Store the parsed structure in an object variable instead of parsing it
            again & again ...
        '''
        self.structure = self.parser.get_structure('PHA-L',
                                                   path)  # Parse the structure once and re-use it in the functions below

    # 3.8 Chains
    def get_number_of_chains(self):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
            Return:
                Number of chains in this structure as integer.
        '''
        n_chains = 0
        for model in self.structure:
            for chain in model:
                n_chains = n_chains + 1
        return n_chains

    # 3.9 Sequence
    def get_sequence(self, chain_id):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the amino acid sequence (single-letter alphabet!) of a given chain (chain_id)
                in a Biopython.PDB structure as a string.
        '''
        sequence = ''

        model_nr = 1
        polypeptide_list = self.ppb.build_peptides(self.structure, model_nr)
        # for polypeptide in polypeptide_list:
        #     print(polypeptide.get_sequence())
        sequence = polypeptide_list[0].get_sequence()

        # chain = self.structure[0][chain_id]
        # print(chain)
        # for residue in chain:
        # print(residue)
        # sequence = sequence + self.aa[residue.get_resname()]

        # for model in self.structure:
        #     for chain in model:
        #         if chain == chain_id:
        #             for residue in chain:
        #                 if is_aa(residue.get_resname(), standard=True):
        #                     sequence.append(three_to_one(residue.get_resname()))
        #                 else:
        #                     sequence.append("X")
        # print(chain_id)
        # sequence = 'SEQWENCE'
        # print(sequence)
        return sequence

    # 3.10 Water molecules
    def get_number_of_water_molecules(self, chain_id):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the number of water molecules of a given chain (chain_id)
                in a Biopython.PDB structure as an integer.
        '''
        n_waters = 0
        print("Hello")
        # print(self.structure)

        chain = self.structure[0][chain_id]
        for residue in chain:
            # print(residue)
            if residue.get_resname() == "HOH":
                n_waters = n_waters + 1
        # n_waters = 12
        return n_waters

    # 3.11 C-Alpha distance
    def get_ca_distance(self, chain_id_1, index_1, chain_id_2, index_2):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id_1 : String (usually in ['A','B', 'C' ...]. The number of chains
                                depends on the specific protein and the resulting structure)
                index_1    : index of a residue in a given chain in a Biopython.PDB structure
                chain_id_2 : String (usually in ['A','B', 'C' ...]. The number of chains
                            depends on the specific protein and the resulting structure)
                index_2    : index of a residue in a given chain in a Biopython.PDB structure

                chain_id_1 and index_1 describe precisely one residue in a PDB structure,
                chain_id_2 and index_2 describe the second residue.

            Return:
                Return the C-alpha (!) distance between the two residues, described by
                chain_id_1/index_1 and chain_id_2/index_2. Round the returned value via int().

            The reason for using two different chains as an input is that also the distance
            between residues of different chains can be interesting.
            Different chains in a PDB structure can either occur between two different proteins
            (Heterodimers) or between different copies of the same protein (Homodimers).
        '''
        ca1 = 0
        ca2 = 0
        doNothing = 0

        chain1 = self.structure[0][chain_id_1]
        chain2 = self.structure[0][chain_id_2]
        residue1 = chain1[index_1]
        residue2 = chain2[index_2]
        try:
            ca_distance = residue1['CA'] - residue2['CA']
        except KeyError:
            doNothing = 1

        # chain1 = self.structure[0][chain_id_1]
        # residue1 = chain1[index_1]
        # # print(residue1)
        # ca1 = residue1['CA']
        #
        # chain2 = self.structure[0][chain_id_2]
        # residue2 = chain2[index_2]
        # ca2 = residue2['CA']

        # chain = self.structure[0][chain_id_1]
        # for residue in chain:
        #     # print(residue)
        #     if residue == chain[index_1]:
        #         ca1 = 0
        #
        # chain = self.structure[0][chain_id_2]
        # for residue in chain:
        #     # print(residue)
        #     if residue == chain[index_2]:
        #         ca2 = 0

        # for model in self.structure:
        #     for chain in model:
        #         if chain == chain_id_1:
        #             for residue in chain:
        #                 # if residue.has_id(("", index_1, “ ”)):
        #                 print(residue)
        #                 ca1 = residue['CA']
        #
        # for model in self.structure:
        #     for chain in model:
        #         if chain == chain_id_2:
        #             for residue in chain:
        #                 # if residue.has_id(("", index_2, “ ”)):
        #                 ca2 = residue['CA']

        # ca_distance = ca1 - ca2
        return int(ca_distance)

    # 3.12 Contact Map
    def get_contact_map(self, chain_id):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return a complete contact map (see description in exercise sheet)
                for a given chain in a Biopython.PDB structure as numpy array.
                The values in the matrix describe the c-alpha distance between all residues
                in a chain of a Biopython.PDB structure.
                Only integer values of the distance have to be given (see below).
        '''

        import numpy as np

        length = 0
        residue_list = []
        chain = self.structure[0][chain_id]
        for residue in chain:
            if residue.get_resname() != "HOH":
                residue_list.append(residue)
                length = length + 1

        # print(residue_list)
        # print(length)

        contact_map = np.zeros((length, length), dtype=np.float32)

        for i in range(0, length):
            residue1 = residue_list[i]
            for j in range(0, length):
                residue2 = residue_list[j]
                try:
                    contact_map[i][j] = residue1['CA'] - residue2['CA']
                except KeyError:
                    doNothing = 1


        # def calc_residue_dist(residue_one, residue_two):
        #     """Returns the C-alpha distance between two residues"""
        #     diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
        #     return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

        return contact_map.astype(np.int)  # return rounded (integer) values

    # 3.13 B-Factors
    def get_bfactors(self, chain_id):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the B-Factors for all residues in a chain of a Biopython.PDB structure.
                The B-Factors describe the mobility of an atom or a residue.
                In a Biopython.PDB structure B-Factors are given for each atom in a residue.
                Calculate the mean B-Factor for a residue by averaging over the B-Factor
                of all atoms in a residue.
                Sometimes B-Factors are not available for a certain residue;
                (e.g. the residue was not resolved); insert np.nan for those cases.

                Finally normalize your B-Factors using Standard scores (zero mean, unit variance).
                You have to use np.nanmean, np.nanvar etc. if you have nan values in your array.
                The returned data structure has to be a numpy array rounded again to integer.
        '''
        import numpy as np

        length = 0
        residue_list = []
        chain = self.structure[0][chain_id]
        for residue in chain:
            if residue.get_resname() != "HOH":
                residue_list.append(residue)
                length = length + 1

        # print(length)

        arr = np.zeros(length, dtype=np.float32)
        b_factors = np.array(arr, dtype=np.float32)
        # print(b_factors.shape)

        # sum = 0
        for i in range(0, length):
            residue1 = residue_list[i]
            sum = 0
            ctr = 0
            # print(i)
            for atom in residue1:
                ctr = ctr + 1
                sum = sum + atom.get_bfactor()
            # print(ctr)
            # print(sum)
            b_factors[i] = sum / ctr

        # print(b_factors)
        # print(b_factors.astype(np.int))
        # print(b_factors.mean(axis=0))
        # print(b_factors.std())
        # print(b_factors.var())
        # print((b_factors - b_factors.mean()) / b_factors.var())
        b_factors = (b_factors - b_factors.mean()) / b_factors.std()

        print(b_factors)
        print(b_factors.astype(np.int))

        return b_factors.astype(np.int)  # return rounded (integer) values


def main():
    print('PDB parser class.')
    # print("xyz")
    protein = PDB_Parser("7ahl.cif")
    # protein.get_number_of_chains()
    # protein.get_sequence('A')
    # print("Hello1")
    # protein.get_number_of_water_molecules('A')
    # protein.get_ca_distance('A', 0, 'A', 0)
    # print("hello2")
    # protein.get_contact_map('A')
    # protein.get_bfactors('A')
    return None


if __name__ == '__main__':
    main()
