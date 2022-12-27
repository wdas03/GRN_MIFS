import math
import pandas as pd

from sympy import *

import itertools

def get_conj_disj_combos_numeric(n):
    nums = [0,1]

    comb = list(itertools.product(nums, repeat=n))
    return comb

def get_conj_disj_combos(inp_symbols):
    combs = []
    for c in get_conj_disj_combos(len(inp_symbols) - 1):
        exp = inp_symbols[0]
        for idx, i in enumerate(c):
            if i == 1:
                exp = exp & inp_symbols[idx+1]
            else:
                exp = exp | inp_symbols[idx+1]
        combs.append(exp)

    combs[0].sub

def gene_wise_dynamics_consistency(pred, actual):
    return sum([1 if i == j else 0 for (i, j) in zip(pred, actual)]) / len(pred)

# Get entropy of list
def entropy(x):
    x_len = len(x)
    x_prob = [x.count(xs) / x_len for xs in x]
    return -sum([xs * math.log(xs, 10) for xs in x_prob])

# Get mutual info between two lists
def mutual_info(x, y):
    return entropy(x) + entropy(y) - entropy(list(zip(x, y)))
    
    