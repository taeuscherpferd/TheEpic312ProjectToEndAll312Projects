#!/usr/bin/python3
# .-------------'```'----....,,__                        _,
# |                               `'`'`'`'-.,.__        .'(
# |                                             `'--._.'   )
# |                                                   `'-.<
# \               .-'`'-.                            -.    `\
#  \               -.o_.     _                     _,-'`\    |
#   ``````''--.._.-=-._    .'  \            _,,--'`      `-._(
#     (^^^^^^^^`___    '-. |    \  __,,..--'                 `
#      `````````   `'--..___\    |`
#              jgs           `-.,’

import itertools
import heapq
from TSPClasses import *
from copy import deepcopy, copy
import numpy as np
import time
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))
try:
    import Queue as pq
except ImportError:
    import queue as pg


class TSPSolver:
    def __init__(self, gui_view):
        self.q = []
        self._scenario = None
        self.nodeCount = 0
        self.nodesPruned = 0
        self.cityList = []

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time()-start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''
    def findLowest( self, cities, route, city ):
        closest = None
        cost = np.inf
        for c in cities:
            if c in route:
                pass
            elif city.costTo( c ) < cost:
                cost = city.costTo( c )
                closest = c
        return closest

    def greedy(self, time_allowance=60.0):
        # Get all the cities and the number of cities
        cities = self._scenario.getCities()
        ncities = len( cities )

        approx = self.defaultRandomTour()
        cost = np.inf
        route = approx['soln']

        # Decides whether we skip checking the solution or not
        skip = False

        # Start from the first city in the list
        start_city = 0

        start_time = time.time()

        # Let each of the cities be the starting place
        # Get the best of all the solution we found
        while (start_city < ncities) and (time.time() - start_time < time_allowance):
            # Keep searching for any cities that we have not been to
            # and add it to the route
            tmp_route = [ cities[start_city] ]
            while len( tmp_route ) < ncities:
                # Find the city with the shortest distance from our current city
                # and of which we have not yet been to.
                closest = self.findLowest( cities, tmp_route, tmp_route[ -1 ] )
                # Change the starting city if we hit a dead end
                if closest == None:
                    skip = True
                    break
                else:
                    tmp_route.append( closest )

            if skip:
                skip = False
            else:
                # If the route we found is a better solution, replace it
                if TSPSolution( tmp_route )._costOfRoute() < cost:
                    route = TSPSolution( tmp_route )
                    cost = TSPSolution( tmp_route )._costOfRoute()

            # Change the starting city to the next one in the list
            start_city += 1


        end_time = time.time()

        results = {}
        results['cost'] = cost
        results['time'] = end_time - start_time
        results['count'] = 0
        results['soln'] = route
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

    def branchAndBound( self, time_allowance=60.0 ):
        # Reset the global values
        self.max_queued = 0
        self.total_nodes = 0
        self.soln_cnt = 0
        self.pruned = 0

        # Get the cities and tyhe number of cities
        cities = self._scenario.getCities()
        ncities = len(cities)

        # Use default calculation for the
        # inital solution approximation
        approx = self.greedy()
        bssf = approx['soln']
        cost = approx['cost']

        # Make the original cost matrix
        ocost = self.makeCostMatrix( cities )
        # Copy the cost matrix, then make the
        # reduced-cost matrix, get the lower
        # bound at the same time
        rcost = ocost.copy()
        lb = self.reduceMatrix( rcost )

        # Use heapq package for the implementation
        queue = []

        start_time = time.time()

        # Push the root problem onto the queue
        heapq.heappush( queue, Node( rcost, lb, [], cities[0] ) )
        self.total_nodes += 1

        while queue and time.time()-start_time < time_allowance:
            # Get the subproblem that has the most potential
            # to be a solution
            node = heapq.heappop( queue )

            # Check to see if we have found a potential
            if node.depth() == ncities:
                # Calculate the cost
                potential = TSPSolution( node.path() ).cost
                # If the potential is better than what we
                # currently have, make potential main
                if potential < cost:
                    bssf = TSPSolution( node.path() )
                    cost = potential
                    self.soln_cnt += 1
            # If the subproblem is still partial,
            # expand the subproblem
            else:
                self.expand( queue, cities, node, cost )

        end_time = time.time()

        # Construct the result struct
        results = {}
        results['cost'] = cost
        results['time'] = end_time - start_time
        results['count'] = self.soln_cnt
        results['soln'] = bssf
        results['max'] = self.max_queued
        results['total'] = self.total_nodes
        results['pruned'] = self.pruned
        return results

    def makeCostMatrix( self, cities ):
        size = len(cities)
        # Construct a n by n matrix
        matrix = np.zeros( (size, size) )
        # Input the values
        for row in range( size ):
            for col in range( size ):
                matrix[ row ][ col ] = cities[ row ].costTo( cities[ col ] )
        return matrix

    def reduceMatrix( self, matrix ):
        lb = 0
        # Get all the indices of all the minimum
        # values for each row, then minus each row
        # according to it's minimum values.
        # Also increment the lower bound accordingly.
        # Skip if the row is all infinity.
        rows = np.argmin( matrix, axis = 1 )
        for r_idx in range( len( matrix ) ):
            min_value = matrix[ r_idx ][ rows[ r_idx ] ]
            if min_value != np.inf:
                lb += min_value
                for c_idx in range( len( matrix[0] ) ):
                    matrix[ r_idx ][ c_idx ] -= min_value

        # Get all the indices of all the minimum
        # values for each column, then minus each
        # column according to it's minimum values.
        # Also increment the lower bound accordingly.
        # Skip if the column is all infinity.
        cols = np.argmin( matrix, axis = 0 )
        for c_idx in range( len( matrix[0] ) ):
            min_value = matrix[ cols[ c_idx ] ][ c_idx ]
            if min_value != np.inf:
                lb += min_value
                for r_idx in range( len( matrix ) ):
                    matrix[ r_idx ][ c_idx ] -= min_value
        return lb

    def expand( self, queue, cities, node, cost ):
        # If the lower bound is higher than bssf cost
        # DROP IT LIKE IT'S HOT
        if node.lowerBound() > cost:
            self.pruned += 1
            return

        # Get the current city information because
        # we are exploring from there.
        # Also get the reduced cost matrix because
        # we need to check whether the next city
        # is reachable from where we are.
        start_city = node.last()
        rmatrix = node.rMatrix()

        # Loop through almost all cities' indices
        # (Skip city 1 because we start from there)
        # We only explore the cities that are reachable
        # from where we are.
        for i in range( 1, len( rmatrix ) ):
            if node.rMatrix()[ start_city._index, i ] != np.inf:
                new_node = self.explore( cities, node, start_city._index, i )
                self.total_nodes += 1

                # We keep all possibilities as long as
                # its lowerbound is smaller than our current cost.
                # Put it in the queue biew biew
                if new_node.lowerBound() <= cost:
                    heapq.heappush( queue, new_node )
                    if len( queue ) > self.max_queued:
                        self.max_queued = len( queue )
                else:
                    self.pruned += 1

    def explore( self, cities, node, start_idx, end_idx ):
        # Now we really need a copy of the matrix
        rmatrix = node.rMatrix().copy()

        # Get the lower bound plus any additional
        # cost to explore the particular city
        lb = node.lowerBound() + rmatrix[ start_idx, end_idx ]
        # Copy copy the path
        path = node.path().copy()

        # Make the row all infinity
        for col in range( len(rmatrix) ):
            rmatrix[ start_idx ][ col ] = np.inf
            # Make the column all infinity
        for row in range( len(rmatrix) ):
            rmatrix[ row ][ end_idx ] = np.inf

        # No more turnning back
        rmatrix[ end_idx ][ start_idx ] = np.inf
        # Reduce the matrix and add on
        # additional lower bound
        lb += self.reduceMatrix( rmatrix )
        return Node( rmatrix, lb, path, cities[ end_idx ] )


    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

    def fancy(self, time_allowance=60.0):
        time_allowance = time_allowance - .5
        start_time = time.time()
        bssf = self.greedy(time_allowance)
        if (time.time() - start_time > time_allowance):
            return self.createResultsDictionary(bssf['cost'], time.time() - start_time, bssf['count'], bssf['soln'], bssf['max'], bssf['total'], bssf['pruned'])

        cities = self._scenario.getCities()
        ncities = len(cities)

        original_matrix = self.make_matrix()

        #print( "---- ORIGINAL ----" )
        #print( original_matrix )

        cost_lookup = [ {} for i in range(ncities) ]

        # Inital cost lookup
        for i in range( 1, ncities ):
            cost_lookup[0][ (i, frozenset()) ] = Subproblem( original_matrix[ i, 0 ], i )

        # 1st layer
        for i in range( 1, ncities ):
            for key, value in cost_lookup[0].items():
                if (time.time() - start_time > time_allowance):
                    return self.createResultsDictionary(bssf['cost'], time.time() - start_time, bssf['count'], bssf['soln'], bssf['max'], bssf['total'], bssf['pruned'])
                    
                if key[0] == i:
                    pass
                else:
                    cost_lookup[1][ (i, key[1].union([key[0]]) ) ] = Subproblem( original_matrix[i][key[0]] + value.cost, i, value.route )

                if (time.time() - start_time > time_allowance):
                    return self.createResultsDictionary(bssf['cost'], time.time() - start_time, bssf['count'], bssf['soln'], bssf['max'], bssf['total'], bssf['pruned'])


        # >1 layers
        for i in range( 2, ncities ):
            for j in range( 1, ncities ):
                # Store each value and key for comparison later
                value_array = []
                key_array = []
                # Loop through each item in the previous layer
                # Skip if the city is in route
                for key, value in cost_lookup[i-1].items():
                    if (time.time() - start_time > time_allowance):
                        return self.createResultsDictionary(bssf['cost'], time.time() - start_time, bssf['count'], bssf['soln'], bssf['max'], bssf['total'], bssf['pruned'])

                    if key[0] == j or j in key[1]:
                        pass
                    else:
                        # Compose the new key for the current layer
                        p = (j, key[1].union([key[0]]) )
                        key_array.append(p)
                        # Construct the new cost and route for the layer
                        sub = Subproblem(original_matrix[j][key[0]] + value.cost, j, value.route )
                        value_array.append(sub)
                        #print("p: {}, sub:{}".format(p, sub) )

                # Input each value
                # If there is already a value, choose whichever has a smaller cost
                for value_idx, key in enumerate(key_array):
                    if (time.time() - start_time > time_allowance):
                        return self.createResultsDictionary(bssf['cost'], time.time() - start_time, bssf['count'], bssf['soln'], bssf['max'], bssf['total'], bssf['pruned'])

                    if cost_lookup[i].get(key, None) == None:
                        cost_lookup[i][key] = value_array[value_idx]
                    else:
                        cost_lookup[i][key] = min( cost_lookup[i][key], value_array[value_idx] )

        # Last layer
        # We need to do this because we have added every city to the set
        # except for the starting city, so we need to add it "manually"
        # We are not really storing the last value in the table
        value_array = []
        for key, value in cost_lookup[ncities-2].items():
            if (time.time() - start_time > time_allowance):
                return self.createResultsDictionary(bssf['cost'], time.time() - start_time, bssf['count'], bssf['soln'], bssf['max'], bssf['total'], bssf['pruned'])

            sub = Subproblem(original_matrix[0][key[0]] + value.cost, 0, value.route)
            value_array.append(sub)
            #print("p: {}, sub:{}".format(p, sub) )

        # Get the index in which the optimal solution is in
        min_idx = np.argmin(value_array)
        optimal = value_array[min_idx]

        #self.printTable( cost_lookup )

        # Get the cost, since we have already computed
        cost = int(optimal.cost)

        # Convert from city indices to actual cities
        route = []
        for i in range( ncities ):
            route.append( cities[optimal.route[i]] )

        #print( "Cost: {}".format(TSPSolution(route)._costOfRoute()) )
        return self.createResultsDictionary(cost, time.time() - start_time, 1, TSPSolution(route), None, None, None)

    def createResultsDictionary(self, cost, time, count, soln, maxNodes, total, pruned):
        results = {}
        results['cost'] = cost
        results['time'] = time
        results['count'] = count
        results['soln'] = soln
        results['max'] = maxNodes
        results['total'] = total
        results['pruned'] = pruned
        return results

    # make the cities into a matrix
    def make_matrix(self):  # O(n^2) # O(n)
        whole_array = []
        cities = self._scenario.getCities()
        for city_from in cities:  # O(n^2) # O(n)
            from_to_array = []
            for city_to in cities:  # O(n) # O(n)
                from_to_array.append(city_from.costTo(city_to))  # O(1) # O(1)
            whole_array.append(from_to_array)
        return np.array(whole_array)  # O(n) # O(n)

    def printTable( self, table ):
        for row in table:
            for key, value in row.items():
                print( "{}: {}".format(key, value) )

class Subproblem:
    def __init__( self, cost, coming_from, route=[] ):
        self.cost = cost
        self.route = [ coming_from ] + route

    def __lt__( self, other ):
        return self.cost < other.cost

    def __str__( self ):
        return "cost: {}, route: {}".format( self.cost, self.route )

#My Node class
class SearchState:
    def __init__(self, matrix, lowerBound, depth, path):
        self.lowerBound = lowerBound
        self.matrix = matrix
        self.depth = depth
        self.path = path

    def __lt__(self, other):
        return (self.lowerBound / self.depth) < other.lowerBound / other.depth

class Node:
    def __init__( self, reduced_matrix, lb, path, dest ):
        self.rmatrix = reduced_matrix
        self.lb = lb
        self.pppath = path + [ dest ]    # Add the new city to the path

    # This is used for the heapq to sort
    # Prioritize depth over lower bound
    def __lt__( self, other ):
        if self.depth() == other.depth():
            return self.lowerBound() < other.lowerBound()
        else:
            return self.depth() > other.depth()

    def path( self ):
        return self.pppath

    def lowerBound( self ):
        return self.lb

    # Get the depth of the subproblem
    def depth( self ):
        return len( self.pppath )

    # Get the city that we are currently in
    def last( self ):
        return self.pppath[ -1 ]

    def rMatrix( self ):
        return self.rmatrix
