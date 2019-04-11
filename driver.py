#!/usr/bin/python3

import numpy as np
import time

from TSPSolver import TSPSolution

def find_not_visited(visited, ncities): # O(n)
    results = []
    if len(visited) == ncities:
        if len(visited) > 0:
            results = [visited[0]]
        return results
    elif len(visited) > ncities:
        return results
    else:
        for i in range(ncities):
            if i not in visited:
                results.append(i)
        return results

def fancy(matrix):
    start_time = time.time()
    results = {}
    ncities = len(matrix[0])

    original_matrix = matrix

    layers = []
    next_layer = []
    origin_tuple = (0,)
    not_visited = find_not_visited(origin_tuple, ncities)
    if len(not_visited) > 0:
        for i in not_visited:
            path_dict = {}
            path_dict[origin_tuple] = []
            if original_matrix[0][i] != np.inf:
                path_dict[origin_tuple].append((i, original_matrix[0][i]))
            if len(path_dict) > 0:
                next_layer.append(path_dict)
    if len(next_layer) > 0:
        layers.append(next_layer)

    # layers[0] = [ {[A]:[(B,7)]}, {[A]:[(C,3)]}, {[A]:[(D,12)]} ]
    # layers[1] = [ {[AB]:[(C,13),(D,21)]}, {[AC]:[(B,11),(D,9)]}, {[AD]:[(B,15),(C,17)]} ]
    # layers[2] = [ {[ABC]:[(D,19)]}, {[ACD]:[(B,12)]}, {[ADB]:[(C,21)]} ]
    # layers[3] = [ {[ABCD]:[(A,28)]}, {[ACDB]:[(A,15)]}, {[ADBC]:[(A,26)]} ]
    #  end layers

    while len(layers) > 0:
        next_layer = []
        for a_dict in layers[-1]:
            for path in a_dict.keys():
                min_tuple = min(a_dict[path], key = lambda a_tuple: a_tuple[1])
                if min_tuple[1] == np.inf:
                    continue
                next_path = path + (min_tuple[0],)
                path_dict = {}
                path_dict[next_path] = []
                not_visited = find_not_visited(next_path, ncities)
                for i in not_visited:
                    path_dict[next_path].append((i, min_tuple[1] + original_matrix[next_path[-1]][i]))
                if len(path_dict[next_path]) > 0:
                    next_layer.append(path_dict)
        if len(next_layer) == 0:
            break
        layers.append(next_layer)

    solution = 'no solution'
    foundTour = False
    if len(layers) > 0:
        min_path = []
        min_length = np.inf
        for a_dict in layers[-1]:
            for path in a_dict.keys():
                a_min_tuple = min(a_dict[path], key=lambda a_tuple: a_tuple[1])
                if a_min_tuple[1] < min_length:
                    min_length = a_min_tuple[1]
                    min_path = path + (a_min_tuple[0],)
                    # min_path.append(a_min_tuple[0])
        if len(min_path) > 0 and min_path[0] == min_path[-1]:
            foundTour = True
            solution = TSPSolution(min_path)


    end_time = time.time()
    results['cost'] = solution.cost if foundTour else np.inf
    results['time'] = end_time - start_time
    results['count'] = len(layers)
    results['soln'] = solution if foundTour else None
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
