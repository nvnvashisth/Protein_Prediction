import numpy as np


"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix or array provided as
           parameters. Further, use those indices when generating or returning
           any matrices or arrays. Failure to do so will most likely result in
           not passing the tests.
EXAMPLE: To access the background frequency for alanine ('A') in the bg_array
         from the get_pssm_background() method, use bg_array[AA_TO_INT['A']].
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY-'
AA_TO_INT = {aa:index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index:aa for index, aa in enumerate(ALPHABET)}
GAP_INDEX = AA_TO_INT['-']


class MSA:


    def __init__(self, path):
        """
        Initialize the MSA class by reading the provided FASTA file. Check the
        sequences for correctness. Pre-calculate any statistics you seem fit.
        
        :param path: Path to the FASTA file containing the MSA sequences.
        """


    def get_pssm_basic(self):
        """
        Calculate a simple PSSM for the underlying MSA. Use only the observed
        amino acids and uniform background frequencies (i.e. P = 0.05).
        Do NOT use sequence weights, redistributed gaps, or pseudocounts!
        Every row in the resulting PSSM corresponds a non-gap position in the
        primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Values that would be -inf must be replaced by -20 in the final PSSM.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).
        
        :return: PSSM as numpy array of shape (L, 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        
        pssm = np.arange(100).reshape(5, 20)
        
        return np.rint(pssm).astype(np.int64)


    def get_pssm_weighted(self):
        """
        Calculate a weighted PSSM for the underlying MSA. Use the same
        statistics as in get_pssm_basic(), but also apply sequence weights.
        Every row in the resulting PSSM corresponds a non-gap position in the
        primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Values that would be -inf must be replaced by -20 in the final PSSM.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).
        
        :return: PSSM as numpy array of shape (L, 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        
        pssm = np.arange(100).reshape(5, 20)
        
        return np.rint(pssm).astype(np.int64)


    def get_pssm_background(self, bg_array):
        """
        Calculate a weighted PSSM for the underlying MSA. Use the same
        statistics as in get_pssm_weighted(), but instead of using uniform
        background frequencies use the correct frequencies provided. Further,
        redistribute the observed gap counts according to those frequencies.
        Every row in the resulting PSSM corresponds a non-gap position in the
        primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Values that would be -inf must be replaced by -20 in the final PSSM.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).
        
        :param bg_array: Amino acid background frequencies as a numpy array.
                         Access the array using the indices from AA_TO_INT.
        
        :return: PSSM as numpy array of shape (L, 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        
        pssm = np.arange(100).reshape(5, 20)
        
        return np.rint(pssm).astype(np.int64)


    def get_pssm_pseudocounts(self, bg_matrix, beta):
        """
        Calculate a weighted PSSM for the underlying MSA. Use the same
        statistics as in get_pssm_background(), but also add pseudocounts.
        Calculate the pseudocounts according to the provided matrix of
        amino acid pair frequencies. Further, calculate the background
        frequencies for all amino acids from the matrix of amino acid pairs.
        Every row in the resulting PSSM corresponds a non-gap position in the
        primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).
        
        :param bg_matrix: Amino acid pair frequencies as numpy array (20, 20).
                          Access the matrix using the indices from AA_TO_INT.
        :param beta: Beta value (float) used to weight the pseudocounts
                     against the observed amino acids in the MSA.
        
        :return: PSSM as numpy array of shape (L, 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        
        pssm = np.arange(100).reshape(5, 20)
        
        return np.rint(pssm).astype(np.int64)


    def get_size(self):
        """
        Return the number of sequences in the MSA and the MSA length, i.e.
        the number of columns in the MSA. This includes gaps.
        
        :return: Tuple of two integers. First element is the number of
                 sequences in the MSA, second element is the MSA length.
        """
        
        num_seq = 666
        msa_len = 1337
        
        return (num_seq, msa_len)


    def get_primary_sequence(self):
        """
        Return the primary sequence of the MSA. In this exercise, the primary
        sequence is always the first sequence in the MSA FASTA file. The
        returned sequence must NOT include gap characters.
        
        :return: String containing the ungapped primary sequence.
        """
        
        ungapped_primary_seq = 'THISISATEST'
        
        return ungapped_primary_seq


    def get_sequence_weights(self):
        """
        Return the calculated sequence weights for all sequences in the MSA.
        The order of weights in the array must be equal to the order of the
        sequences in the MSA FASTA file.
        
        :return: Numpy array (dtype=numpy.float64) containing the weights for
                 all sequences in the MSA.
        """
        
        weights = np.arange(666) + 42
        
        return weights.astype(np.float64)


    def get_number_of_observations(self):
        """
        Return the estimated number of independent observations in the MSA.
        
        :return: Estimate of independent observation (dtype=numpy.float64).
        """
        
        num_obs = 3.14159265358979323846264338327950288419716939937510
        
        return np.float64(num_obs)