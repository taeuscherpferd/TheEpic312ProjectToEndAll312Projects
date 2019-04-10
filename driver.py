#!/usr/bin/python3

import numpy as np
import time

from TSPSolver import TSPSolver


def fancy(matrix):
    results = {}
    ncities = len(matrix[0])

    original_matrix = matrix
    cost_lookup = [{} for i in range(ncities)]

    # Inital cost lookup
    for i in range(1, ncities):
        cost_lookup[0][(i, frozenset())] = original_matrix[i, 0]

    # 1st layer
    for i in range(1, ncities):
        for key in cost_lookup[0]:
            if key[0] == i:
                pass
            else:
                cost_lookup[1][(i, key[1].union(
                    [key[0]]))] = original_matrix[i][key[0]] + cost_lookup[0][key]

    # >1 layers
    for i in range(2, ncities):
        for j in range(1, ncities):
            values = np.array([])
            for key in cost_lookup[i-1]:
                if key[0] == j:
                    pass
                else:
                    cost = original_matrix[j][key[0]
                                                ] + cost_lookup[i-1][key]
                    np.append(values, cost)

    end_time = time.time()
    results['cost'] = np.inf
    results['time'] = end_time - start_time
    results['count'] = np.inf
    results['soln'] = None
    results['max'] = None
    results['total'] = None
    results['pruned'] = None
    return results

if __name__ == "__main__":
    A = np.array([[np.inf, 7, 3, 12],
                  [3, np.inf, 6, 14],
                  [5, 8, np.inf, 6],
                  [9, 3, 5, np.inf]])

    print(A)
    print("*(**(*(**********")

    print(fancy(A))
