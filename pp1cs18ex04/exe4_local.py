import numpy as np
from numpy import unravel_index


class LocalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns, string2 is used to index rows
        """
        self.string1 = string1
        self.string2 = string2
        self.gap_penalty = gap_penalty
        self.substitution_matrix = matrix
        self.score_matrix = np.zeros((len(self.string2) + 1, len(self.string1) + 1), dtype=np.int)
        #print(self.score_matrix)
        self.index1 = []
        self.index2 = []
        #self.globalalignB = ""
        self.max = 0
        self.pos = None
        #print(self.align())
        self.align()

    def align(self):
        """
        Align given strings using the Smith-Waterman algorithm.
        NB: score matrix and the substitution matrix are different matrices!
        """
        for col in range(len(self.string1)+1):
            #print(self.score_matrix[0][col])
            self.score_matrix[0][col] = 0
        for row in range(len(self.string2)+1):
            #print(self.score_matrix[row][0])
            self.score_matrix[row][0] = 0

        for i in range(len(self.string2) + 1):
            for j in range(len(self.string1) + 1):
                if i != 0 and j != 0:
                    self.score_matrix[i][j] = max(0, self.score_matrix[i - 1][j - 1] + self.substitution_matrix[self.string2[i-1]][self.string1[j-1]],\
                                              self.score_matrix[i - 1][j] + self.gap_penalty, \
                                              self.score_matrix[i][j - 1] + self.gap_penalty)
                    if self.score_matrix[i][j] > self.max:
                        self.max = self.score_matrix[i][j]
                        self.pos = (i,j)

        #print(self.score_matrix)
        #print(self.pos)


    def has_alignment(self):
        """
        :return: True if a local alignment has been found, False otherwise
        """

        c, a = self.get_alignment()
        if c == "" and a == "":
            return False
        else:
            return True


    def get_alignment(self):
        """
        :return: alignment represented as a tuple of aligned strings


        """
        (index_row, index_col) = unravel_index(self.score_matrix.argmax(), self.score_matrix.shape)
        #alignments = []
        AlignmentA = ""
        AlignmentB = ""
        i = index_row
        j = index_col
        while i > 0 or j > 0:
            # print( self.substituion_matrix[self.string2[i]][self.string1[j]])
            if i > 0 and j > 0 and self.score_matrix[i][j] == self.score_matrix[i - 1][j - 1] + \
                    self.substitution_matrix[self.string2[i - 1]][self.string1[j - 1]]:
                AlignmentA = self.string2[i - 1] + AlignmentA
                AlignmentB = self.string1[j - 1] + AlignmentB
                self.index2.append(i-1)
                self.index1.append(j-1)
                i = i - 1
                j = j - 1
            elif i > 0 and self.score_matrix[i][j] == self.score_matrix[i - 1][j] + self.gap_penalty:
                AlignmentA = self.string2[i - 1] + AlignmentA
                #self.index1.append(i-1)
                AlignmentB = "-" + AlignmentB
                i = i - 1
            else:
                AlignmentA = "-" + AlignmentA
                AlignmentB = self.string1[j - 1] + AlignmentB
                #self.index2.append(j - 1)
                j = j - 1



        count1 = len(AlignmentB)
        index = 0
        #index_1 = 0
        while count1>0:
           while AlignmentB[0] == '-':

                AlignmentB = AlignmentB[:index] + AlignmentB[index + 1:]
                AlignmentA = AlignmentA[:index] + AlignmentA[index + 1:]
                #print(AlignmentB,AlignmentA)
           count1 -= 1


        count = len(AlignmentA)
        #index_1 = 0

        while count>0:
            while AlignmentA[0] == '-':
                AlignmentB = AlignmentB[:index] + AlignmentB[index + 1:]
                AlignmentA = AlignmentA[:index] + AlignmentA[index + 1:]
                #print(AlignmentB,AlignmentA)
            count -= 1


        return (AlignmentB,AlignmentA)

    def is_residue_aligned(self, string_number, residue_index):
        """
        :param string_number: number of the string (1 for string1, 2 for string2) to check
        :param residue_index: index of the residue to check
        :return: True if the residue with a given index in a given string has been aligned
                 False otherwise
        """
        self.get_alignment()
        #print(self.index1)
        #print(self.index2)

        if string_number == 1:
            index = self.index1
        else:
            index = self.index2

        if residue_index in index:
            return True
        else:
            return False



