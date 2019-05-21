import numpy as np


"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix or array provided as
           parameters. Further, use those indices when generating or returning
           any matrices or arrays. Failure to do so will most likely result in
           not passing the tests.
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_INT = {aa:index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index:aa for index, aa in enumerate(ALPHABET)}


def get_blast_seeds(pssm, min_score):
    """
    Calculate all 3-gram BLAST seeds that have at least the
    minimum score for a given query PSSM. Return a dictionary
    with the seeds as keys, and a tuple with the positions of
    the seeds in the query as values. The positions must be in
    increasing order.
    
    :param pssm: PSSM as a numpy array of shape (L, 20).
                 L = length of the query sequence.
                 Each row corresponds to one position in the query.
                 Each column to one of the 20 amino acids.
    :param min_score: Minimum score for the BLAST seed
    
    :return: Dictionary with seeds as keys, tuple of positions as values.
             Keys are strings, values are tuples of integers.
    """
    
    seeds = {'AAA': (0, 2, 4)}
    
    return seeds