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
        cities = self._scenario.getCities().copy()
        ncities = len( cities )

        foundTour = False

        start_time = time.time()

        # while there are cities left to add
        while not foundTour and time.time() - start_time < time_allowance:
            for i in range( ncities ):
                tmp_route = [ cities[i] ]
                while len( tmp_route ) < ncities:
                    closest = self.findLowest( cities, tmp_route, tmp_route[ -1 ] )
                    if closest == None:
                        break
                    else:
                        tmp_route.append( closest )

                if TSPSolution( tmp_route )._costOfRoute() < np.inf:
                    route = TSPSolution( tmp_route )
                    cost = TSPSolution( tmp_route )._costOfRoute()
                    foundTour = True

        end_time = time.time()

        results = {}
        results['cost'] = cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = 0
        results['soln'] = route if foundTour else None
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

    def branchAndBound(self, time_allowance=60.0):
        results = {}
        qCounts = []
        numSolutionsFound = 0
        startTime = time.time()
        self.nodeCount = 0
        self.nodesPruned = 0
        self.cityList = deepcopy(self._scenario.getCities())
        cityCount = len(self.cityList)
        cityMatrix = [[0 for y in range(cityCount)] for x in range(cityCount)]

        #Use the random path to start off with
        bssf = self.defaultRandomTour(time_allowance)['soln']

        #Generate the starting matrix O(n^2)
        for row in range(cityCount):
            for col in range(cityCount):
                if col == row:
                    cityMatrix[row][col] = math.inf
                    continue
                cityMatrix[row][col] = self.cityList[row].costTo(self.cityList[col])

        #Calculate the initial reduced matrix
        #Time: O(n^2)
        #Space: O(n^2)
        lowerBoundMatrix, lowerBound = self.calcLowerBoundMatrix(np.array(cityMatrix), cityCount)

        startNode = SearchState(lowerBoundMatrix, lowerBound, 0, [0])

        self.nodeCount += 1
        numSolutionsFound += 1
        heapq.heappush(self.q, startNode)

        elapsedTime = time.time() - startTime

        #Quit on time out or if queue is empty.
        while elapsedTime < time_allowance and not len(self.q) == 0:
            qCounts.append(len(self.q))
            curLowestMatrix = heapq.heappop(self.q)

            #If the current path is complete update the BSSF
            if (len(curLowestMatrix.path) == cityCount):
                numSolutionsFound += 1
                actualCities = [self.cityList[x] for x in curLowestMatrix.path]
                bssf = TSPSolution(actualCities)

            #Expand the children
            self.expandNode(curLowestMatrix, cityCount, bssf)
            elapsedTime = time.time() - startTime

        results['cost'] = bssf.cost
        results['time'] = time.time() - startTime
        results['count'] = numSolutionsFound #num solutions found
        results['soln'] = bssf
        results['max'] = max(qCounts)
        results['total'] = self.nodeCount
        results['pruned'] = self.nodesPruned
        return results

    def expandNode(self, curLowestMatrix, size, bssf):
        matrix = curLowestMatrix.matrix
        curLowerBound = curLowestMatrix.lowerBound
        curPath = curLowestMatrix.path
        curDepth = curLowestMatrix.depth
        curLocation = curPath[-1]
        npArray = np.array(matrix)

        #Time: O(n^3)
        #Space O(n^2)
        for col in range(size):
            contenderMatrix = npArray.copy()

            #Make the current row infinity O(n)
            for x in range(size):
               contenderMatrix[curLocation][x] = math.inf

            curCost = matrix[curLocation][col]
            if curCost == math.inf:
               continue

            #Make the current column infinity O(n)
            for x in range(size):
               contenderMatrix[x][col] = math.inf

            #Calculate the reduced matrix
            #Time: O(n^2)
            #Space: O(n^2)
            contenderMatrix, contLowerbound = self.calcLowerBoundMatrix(contenderMatrix, size)
            contLowerbound += curCost
            contLowerbound += curLowerBound

            newPath = deepcopy(curPath)
            newPath.append(col)

            contNode = SearchState(contenderMatrix, contLowerbound, curDepth + 1, newPath)
            self.nodeCount += 1

            #Checks against the best solution so far. If better, add to queue else, prune
            if contLowerbound < bssf.cost:
                heapq.heappush(self.q, contNode)
            else:
                self.nodesPruned += 1

    def calcLowerBoundMatrix(self, matrix, size):
        lowerbound = 0

        #Time: O(n^2)
        #Space: O(n^2)
        for i in range(size):
              minim = np.min(matrix[i, :])
              if minim == math.inf:
                continue
              lowerbound += minim
              matrix[i, :] -= minim

        #Time: O(n^2)
        #Space: O(n^2)
        for i in range(size):
              minim = np.min(matrix[:, i])
              if minim == math.inf:
                continue
              lowerbound += minim
              matrix[:, i] -= minim

        return matrix, lowerbound

    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

    def fancy(self, time_allowance=60.0):
        results = {}
        start_time = time.time()

        cities = self._scenario.getCities()
        ncities = len(cities)

        original_matrix = self.make_matrix()
        cost_lookup = [ {} for i in range(ncities) ]

        ps = self.makePowerset(ncities)

        # Inital cost lookup
        for i in range( 1, ncities ):
            cost_lookup[0][ (i, frozenset()) ] = Subproblem( original_matrix[ i, 0 ], i )

        # 1st layer
        for i in range( 1, ncities ):
            for key in cost_lookup[0]:
                if key[0] == i:
                    pass
                else:
                    cost_lookup[1][ (i, key[1].union([key[0]]) ) ] = Subproblem( original_matrix[i][key[0]] + cost_lookup[0][key].cost, i, cost_lookup[0][key].route )


        for row in cost_lookup:
            print(row)
        # >1 layers
        for i in range( 2, ncities ):
            for j in range( 1, ncities ):
                value_array = []
                key_array = []
                for key in cost_lookup[i-1]:
                    if key[0] == j or j in key[1]:
                        pass
                    else:
                        p = (j, key[1].union([key[0]]) )
                        key_array.append(p)
                        sub = Subproblem(original_matrix[j][key[0]] + cost_lookup[i-1][key].cost, j, cost_lookup[i-1][key].route )
                        value_array.append(sub)

                min_idx = np.argmin(value_array)
                print("min_idx: {}, k:{}, v: {}".format(min_idx,key_array[min_idx], value_array[min_idx]) )
                cost_lookup[i][key_array[min_idx]] = value_array[min_idx]



        end_time = time.time()
        results['cost'] = math.inf
        results['time'] = end_time - start_time
        results['count'] = math.inf
        results['soln'] = None
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
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

    def makePowerset(self, ncities):
        powerset = [frozenset()]
        for i in range(ncities):
            length=len(powerset)
            for j in range(length):
                subset = powerset[j].union([i])
                powerset.append(subset)

        for subset in powerset:
            if 0 in subset:
                powerset.remove(subset)

        for subset in powerset:
            if len(subset) < 2:
                powerset.remove(subset)

        return sorted(powerset[1:], key=len)

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
