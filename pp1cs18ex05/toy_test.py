import pytest
import blast
import numpy as np
from pssm import MSA


pssm_1 = np.array([[ 1, -1,  0,  0, -3,  0, -1, -3,  0, -3, -2,  0, -1,  0, -1,  4,  1, -2, -3, -2],
                   [-1, -4,  1,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -3, -3, -2],
                   [ 2, -2,  0,  4, -3, -1, -1, -2,  0, -2, -2, -1, -1,  1, -1,  0, -1, -1, -3, -2],
                   [ 0, -2,  1,  0, -3, -1,  0, -3,  0, -3, -2,  5, -2,  0, -1,  3,  0, -2, -3, -2]])
score_1 = 12
seeds_1 = {'EAN': (1,), 'EEN': (1,), 'EES': (1,), 'SEE': (0,)}


def test_get_blast_seeds():
    seeds = blast.get_blast_seeds(pssm_1, score_1)
    num_seeds = len(seeds)
    
    check_size = num_seeds == len(seeds_1)
    assert check_size, 'Incorrect number of BLAST seeds: %d' % num_seeds
    
    check_keys = len(set(seeds.keys()) ^ set(seeds_1.keys())) == 0
    assert check_keys, "Wrong BLAST seeds returned."
    
    check_equal = seeds == seeds_1
    assert check_equal, 'Wrong BLAST seed positions returned.'


msa_1 = MSA('msa.fasta')


def test_get_size():
    size = msa_1.get_size()
    
    check_type = (isinstance(size, tuple) and len(size) == 2)
    assert check_type, 'Return value is not a tuple with two elements.'
    
    check_seqs = size[0] == 4
    assert check_seqs, 'Incorrect number of sequences returned.'
    
    check_cols = size[1] == 5
    assert check_cols, 'Incorrect MSA length returned.'


def test_get_primary_sequence():
    seq = msa_1.get_primary_sequence()
    
    check_type = isinstance(seq, str)
    assert check_type, 'Return value is not a string.'
    
    check_equal = seq.upper() == 'SEAN'
    assert check_equal, 'Incorrect primary sequence returned.'


def test_get_sequence_weights():
    weights = msa_1.get_sequence_weights()
    
    check_type = (hasattr(weights, 'shape') and weights.dtype == np.float64)
    assert check_type, 'Return value not a numpy array of type numpy.float64.'
    
    check_shape = weights.shape == (4,)
    assert check_shape, 'Returned sequence weights have incorrect shape.'
    
    solution = np.array([0.66666667, 0.66666667, 1.0, 0.66666667])
    check_equal = np.allclose(weights, solution)
    assert check_equal, 'Incorrect sequence weights returned.'


def test_get_number_of_observations():
    num_obs = msa_1.get_number_of_observations()
    
    check_type = (hasattr(num_obs, 'dtype') and num_obs.dtype == np.float64)
    assert check_type, 'Return value not of type numpy.float64.'
    
    check_equal = np.isclose(num_obs, 1.6)
    assert check_equal, 'Incorrect estimate of independent observation.'


def check_pssm(test_pssm, master_pssm):
    check_type = (hasattr(test_pssm, 'shape') and test_pssm.dtype == np.int64)
    assert check_type, 'Return value not a numpy array of type numpy.int64.'
    
    check_shape = test_pssm.shape == master_pssm.shape
    assert check_shape, 'Returned PSSM has incorrect shape.'
    
    check_equal = np.array_equal(test_pssm, master_pssm)
    assert check_equal, 'Returned PSSM contains incorrect values.'


pssm_2 = np.array([[-20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,   9, -20, -20,  -20, -20],
                   [-20, -20, -20,   9, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,  -20, -20],
                   [  7, -20, -20,   7, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,  -20, -20],
                   [-20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,   7, -20, -20, -20,   7, -20, -20,  -20, -20]])
pssm_3 = np.array([[-20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,   9, -20, -20,  -20, -20],
                   [-20, -20, -20,   9, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,  -20, -20],
                   [  6, -20, -20,   7, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,  -20, -20],
                   [-20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20,   7, -20, -20, -20,   6, -20, -20,  -20, -20]])


def test_pssms():
    check_pssm(msa_1.get_pssm_basic(), pssm_2)
    check_pssm(msa_1.get_pssm_weighted(), pssm_3)
    """
    Since pssm_background and pssm_pseudocounts cannot be tested without exposing
    the background array (and calculating this array is part of the exercise),
    those will only be tested online.
    """