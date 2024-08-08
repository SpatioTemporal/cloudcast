#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""

bayes.py
~~~~~~~~

$ python bayes.py

"""
# Standard Imports
import os

# Third-Party Imports
# import numpy as np
# import matplotlib.pyplot as plt
# import scipy.stats as st

# Local Imports

##
# List of Public objects from this module.
__all__ = ['bayes']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

r"""
use_this = 0
    Sample size                                       : 210
    Ratio of B to C                                   : 0.047619047619047616

    Sample size B                                     : 10
    Sample size C                                     : 200

    P(H) == P(B) Prior Probability                    : 0.047619047619047616
    P(C) Marginal Probability                         : 0.9523809523809523

    Alpha size B                                      : 4
    Alpha size C                                      : 20

    P(alpha | B) = P(E|H) Likelihood Probability      : 0.4
    P(alpha | C) = P(E | not H)                       : 0.1
    P(E) == P(alpha) Evidence Probability             : 0.11428571428571428
    Alpha size                                        : 24

    P(B | alpha) = P(H|E) Posterior Probability       : 0.16666666666666666
    P(B | alpha) = P(H|E) Alt                         : 0.16666666666666666
    P(B | alpha) = P(H|E) Alt                         : 0.16666666666666669

use_this = 1

    Sample size                                       : 210
    Ratio of B to C                                   : 0.5238095238095238

    Sample size B                                     : 110
    Sample size C                                     : 100

    P(H) == P(B) Prior Probability                    : 0.5238095238095238
    P(C) Marginal Probability                         : 0.47619047619047616

    Alpha size B                                      : 55
    Alpha size C                                      : 10

    P(alpha | B) = P(E|H) Likelihood Probability      : 0.5
    P(alpha | C) = P(E | not H)                       : 0.1
    P(E) == P(alpha) Evidence Probability             : 0.30952380952380953
    Alpha size                                        : 65

    P(B | alpha) = P(H|E) Posterior Probability       : 0.8461538461538461
    P(B | alpha) = P(H|E) Alt                         : 0.8461538461538461
    P(B | alpha) = P(H|E) Alt                         : 0.8461538461538461
"""


##
# Setup Sample distribution
n_sample = 210

##
# Create 2 sub-populations
use_this = 1
ratio_bc = [1 / 21, 110/210][use_this]
ratio_cb = [20 / 21, 100/210][use_this]
n_sample_b = int(ratio_bc * n_sample)
n_sample_c = n_sample - n_sample_b
print(f"{'Sample size':<50s}: {n_sample}")
print(f"{'Ratio of B to C':<50s}: {ratio_bc}")
print(f"\n{'Sample size B':<50s}: {n_sample_b}")
print(f"{'Sample size C':<50s}: {n_sample_c}")

##
# Prior/Marginal Probability    P(H) == P(B)
p_b = ratio_bc
p_c = ratio_cb
print(f"\n{'P(H) == P(B) Prior Probability':<50s}: {p_b}")
print(f"{'P(C) Marginal Probability':<50s}: {p_c}")

##
# Property alpha incidence
#   P(alpha | B) = P(E|H) Likelihood Probability
p_alpha_b = [0.4, 55/110][use_this]
n_alpha_b = int(p_alpha_b * n_sample_b)
#  P(alpha | C) == P(Alpha | not B) = P(E | not H)
p_alpha_c = [0.1, 10/100.][use_this]
n_alpha_c = int(p_alpha_c * n_sample_c)
# Evidence Probability == P(E)
p_alpha = (n_alpha_b + n_alpha_c) / n_sample
n_alpha = n_alpha_b + n_alpha_c
print(f"\n{'Alpha size B':<50s}: {n_alpha_b}")
print(f"{'Alpha size C':<50s}: {n_alpha_c}")
print(f"\n{'P(alpha | B) = P(E|H) Likelihood Probability':<50s}: {p_alpha_b}")
print(f"{'P(alpha | C) = P(E | not H)':<50s}: {p_alpha_c}")
print(f"{'P(E) == P(alpha) Evidence Probability':<50s}: {p_alpha}")
print(f"{'Alpha size':<50s}: {n_alpha}")

##
# Posterior Probability P(H|E) = P(B | alpha)
p_b_alpha_a = n_alpha_b / n_alpha
print(f"\n{'P(B | alpha) = P(H|E) Posterior Probability':<50s}: {p_b_alpha_a}")

p_b_alpha_num = p_b * p_alpha_b
p_b_alpha_dem = p_b_alpha_num + (p_c * p_alpha_c)
p_b_alpha_b = p_b_alpha_num / p_b_alpha_dem
print(f"{'P(B | alpha) = P(H|E) Alt':<50s}: {p_b_alpha_b}")

p_b_alpha_c = p_b_alpha_num / p_alpha
print(f"{'P(B | alpha) = P(H|E) Alt':<50s}: {p_b_alpha_c}")

