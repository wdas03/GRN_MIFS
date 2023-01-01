import math
import pandas as pd

from sympy import *

import itertools

# Get conjunctive and disjunctive expression combinations
def get_conj_disj_combos_numeric(n):
    nums = [0,1]

    comb = list(itertools.product(nums, repeat=n))
    return comb

def get_conj_disj_combos(inp_symbols):
    combs = []
    for c in get_conj_disj_combos_numeric(len(inp_symbols) - 1):
        exp = inp_symbols[0]
        for idx, i in enumerate(c):
            if i == 1:
                exp = exp & inp_symbols[idx+1]
            else:
                exp = exp | inp_symbols[idx+1]
        combs.append(exp)

    return combs

# Get time-series prediction based on boolean expression and target node
def get_pred(target_node, bool_exp, data):
    pred = []
    for i in range(len(data['time']) - 1):
        inp = bool_exp.atoms()
        inp_dict = {}
        for j in inp:
            str_j = str(j)
            inp_dict[str_j] = True if data[str_j][i] == 1 else False
        pred.append(1 if bool_exp.subs(inp_dict) else 0)
        #print(inp_dict, bool_exp, bool_exp.subs(inp_dict))
        
    return pred

# Get max consistency given target node and inp symbols
def get_max_consistency(target_node, inp, data):
    max_consistency = 0
    for bool_exp in get_conj_disj_combos(inp):
        pred = get_pred(target_node, bool_exp, data)
        consistency = gene_wise_dynamics_consistency(pred, data[target_node][1:])

        max_consistency = max(max_consistency, consistency)
    
    return max_consistency

# Gene-wise dynamics consistency
def gene_wise_dynamics_consistency(pred, actual):
    return sum([1 if i == j else 0 for (i, j) in zip(pred, actual)]) / len(pred)

# Get entropy of list
def entropy(x):
    x_len = len(x)
    x_prob = [x.count(xs) / x_len for xs in x]
    return -sum([xs * math.log(xs, 2) for xs in x_prob])

# Get mutual info between two lists
def mutual_info(x, y):
    return entropy(x) + entropy(y) - entropy(list(zip(x, y)))

# Run MIFS SWAP routine
def mifs_swap(data, k):
    node_inp = {}
    final_inp = {}

    gene_names = list(data.keys())[1:]

    for gene in gene_names:
        inp_genes = [i for i in gene_names if i != gene]
        
        curr_inp = []

        max_mi = max ([(i, mutual_info(data[gene], data[i])) for i in inp_genes], key = lambda x: x[1])
        curr_inp.append(max_mi)
        inp_genes = [g for g in inp_genes if g != max_mi[0]]

        while len(curr_inp) < k:
            new_max = max ([(w, mutual_info(data[gene], data[w]) - sum([mutual_info(data[w], data[s[0]]) for s in curr_inp])) for w in inp_genes], 
                            key = lambda x: x[1])
            curr_inp.append(new_max)
            inp_genes = [g for g in inp_genes if g != new_max[0]]
        
        node_inp[gene] = curr_inp

    for gene in gene_names:
        selected_nodes = [i[0] for i in sorted(node_inp[gene], key=lambda x: x[1])]
        unselected_nodes = list(set(gene_names) - set(selected_nodes) - set([gene]))

        #print(gene, selected_nodes, unselected_nodes)

        max_consistency = get_max_consistency(gene, symbols(selected_nodes), data)
        for i, s_node in enumerate(selected_nodes):
            for j, u_node in enumerate(unselected_nodes):
                new_inp = [s for s in selected_nodes if s != s_node]
                new_inp.append(u_node)

                #print(new_inp)

                new_consistency = get_max_consistency(gene, symbols(new_inp), data)
                if new_consistency > max_consistency:
                    selected_nodes = new_inp
                    unselected_nodes = [u for u in unselected_nodes if u != u_node]
                    unselected_nodes.append(s_node)
                    max_consistency = new_consistency

        print(gene, selected_nodes, max_consistency)
        final_inp[gene] = selected_nodes
    
    return final_inp
    
    