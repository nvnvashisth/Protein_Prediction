import numpy as np


class GlobalAlignment:
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
        self.substituion_matrix = matrix
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int)
        self.align()

    def align(self):
        """
        Align given strings using the Needleman-Wunsch algorithm,
        store the alignments and the score matrix used to compute those alignments.
        NB: score matrix and the substitution matrix are different matrices!

        """
        for x in range(len(self.string1)+1):
            self.score_matrix[0][x] = self.gap_penalty * x
        for y in range(len(self.string2)+1):
            self.score_matrix[y][0] = self.gap_penalty * y

        for a in range(len(self.string2)+1):
            for b in range(len(self.string1)+1):
                if a != 0 and b != 0:
                    self.score_matrix[a][b] = max(self.score_matrix[a - 1][b - 1] + self.substituion_matrix[self.string2[a-1]][self.string1[b-1]],
                                                         self.score_matrix[a - 1][b] + self.gap_penalty, self.score_matrix[a][b - 1] + self.gap_penalty)

        #print(self.score_matrix)

    def get_best_score(self):
        """
        :return: the highest score for the aligned strings, int

        """
        return self.score_matrix[len(self.string2)][len(self.string1)]

    def get_number_of_alignments(self):
        """
        :return: number of found alignments with the best score

        """
        if self.get_alignments():
            num = (len(self.get_alignments()))
            return num

        else:
            return 0

    def get_alignments(self):
        """
        :return: list of alignments, where each alignment is represented
                 as a tuple of aligned strings
        """
        a = len(self.string2)
        b = len(self.string1)
        if a == 1 or b == 1:
            if self.string1[0] != self.string2[0]:
                return None
        stack_list = []
        stack_list.append(("", "", a, b))
        final_answers = []

        while len(stack_list) > 0:
            # print(stack_list)
            AlignmentA, AlignmentB, m, n = stack_list.pop()
            if m == 0 and n == 0:
                final_answers.append((AlignmentB, AlignmentA))
                continue
            # print( self.substituion_matrix[self.string2[m]][self.string1[n]])
            if m > 0 and n > 0 and self.score_matrix[m][n] == self.score_matrix[m - 1][n - 1] + \
                    self.substituion_matrix[self.string2[m - 1]][self.string1[n - 1]]:
                AlignmentA_1 = self.string2[m - 1] + AlignmentA
                AlignmentB_1 = self.string1[n - 1] + AlignmentB
                stack_list.append((AlignmentA_1, AlignmentB_1, m - 1, n - 1))

            if n > 0 and self.score_matrix[m][n] == self.score_matrix[m][n - 1] + self.gap_penalty:
                AlignmentA_2 = "-" + AlignmentA
                AlignmentB_2 = self.string1[n - 1] + AlignmentB
                stack_list.append((AlignmentA_2, AlignmentB_2, m, n - 1))

            if m > 0 and self.score_matrix[m][n] == self.score_matrix[m - 1][n] + self.gap_penalty:
                AlignmentA_3 = self.string2[m - 1] + AlignmentA
                AlignmentB_3 = "-" + AlignmentB
                stack_list.append((AlignmentA_3, AlignmentB_3, m - 1, n))

        #print(len(final_answers))
        return final_answers

    def get_score_matrix(self):
        """
        :return: matrix built during the alignment process as an np.array

        """
        return self.score_matrix
        # return np.array([
        #     [0, -1, -2, -3, -4, -5, -6],
        #     [-1, 1, 0, -1, -2, -3, -4],
        #     [-2, 0, 2, 1, 0, -1, -2],
        #     [-3, -1, 1, 3, 2, 1, 0],
        #     [-4, -2, 0, 2, 4, 3, 2],
        #     [-5, -3, -1, 1, 3, 4, 3],
        #     [-6, -4, -2, 0, 2, 3, 4],
        #     [-7, -5, -3, -1, 1, 2, 4]
        # ])
identity_matrix = {
    'A': {'A': 1, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'R': {'A': 0, 'R': 1, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'N': {'A': 0, 'R': 0, 'N': 1, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'D': {'A': 0, 'R': 0, 'N': 0, 'D': 1, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'C': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 1, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'E': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 1, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'Q': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 1, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'G': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 1, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'H': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 1, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'I': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 1, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'L': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 1, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'K': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'M': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 1, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'F': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 1, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'P': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 1, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'S': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 1, 'T': 0, 'W': 0, 'Y': 0, 'V': 0},
    'T': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 1, 'W': 0, 'Y': 0, 'V': 0},
    'W': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 1, 'Y': 0, 'V': 0},
    'Y': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 1, 'V': 0},
    'V': {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 1}
}


ga = GlobalAlignment(
        'A','C',-1,
        identity_matrix)
print(ga.get_alignments())
#print(ga.is_residue_aligned(1,2))


