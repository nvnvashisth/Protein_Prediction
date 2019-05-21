import pytest
import numpy as np
import exe4_global
import exe4_local

"""
You can use this toy test to verify that your code doesn't have
syntax errors and adheres to the input-output specifications.
To execute this toy test, run `pytest` in the project directory.
To install `pytest`, run `conda install pytest` or
`pip install pytest` (depending on your python installation).
"""

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


def test_that_code_can_be_run():
    ga = exe4_global.GlobalAlignment(
        'ADMINS',
        'ADMIRES',
        -1,
        identity_matrix
    )

    assert ga.get_best_score() == 4, 'Wrong high score'
    assert ga.get_number_of_alignments() == 2, 'Wrong number of alignments'
    assert set(ga.get_alignments()) == {('ADMI-NS', 'ADMIRES'), ('ADMIN-S', 'ADMIRES')}, 'Wrong alignments'
    assert np.array_equal(ga.get_score_matrix(),
                          np.array([
                              [0, -1, -2, -3, -4, -5, -6],
                              [-1, 1, 0, -1, -2, -3, -4],
                              [-2, 0, 2, 1, 0, -1, -2],
                              [-3, -1, 1, 3, 2, 1, 0],
                              [-4, -2, 0, 2, 4, 3, 2],
                              [-5, -3, -1, 1, 3, 4, 3],
                              [-6, -4, -2, 0, 2, 3, 4],
                              [-7, -5, -3, -1, 1, 2, 4]
                          ])), 'Wrong score matrix'

    la = exe4_local.LocalAlignment(
        'ADAC',
        'DAC',
        -1,
        identity_matrix
    )

    assert la.has_alignment(), 'Wrong value returned from has_alignments()'
    assert la.get_alignment() == ('DAC', 'DAC'), 'Wrong alignment returned'
    assert not la.is_residue_aligned(1, 0), 'Wrong value returned from is_residue_aligned'
