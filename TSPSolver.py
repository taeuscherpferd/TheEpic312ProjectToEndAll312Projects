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

    def greedy(self, time_allowance=60.0):
        start_time = time.time()
        results = {}
        count = 0
        cities = self._scenario.getCities().copy()
        bssf = None
        foundTour = False
        route = []
        route_cost = 0
        route.append(cities[0])
        cities.remove(cities[0])
        prev_city = cities[0]
        # while there are cities left to add
        while len(cities) > 0 and time.time() - start_time < time_allowance:
            # add the next city with the lowest path cost
            lowest_cost = math.inf
            lowest_cost_city = 0
            for city in cities:
                cost = prev_city.costTo(city)
                if cost <= lowest_cost:
                    lowest_cost = cost
                    lowest_cost_city = city
            cities.remove(lowest_cost_city)
            route.append(lowest_cost_city)
            route_cost += lowest_cost
            prev_city = lowest_cost_city
        bssf = TSPSolution(route)
        if bssf.cost < np.inf:
            # Found a valid route
            foundTour = True
            count += 1
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
        pass

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

#My Node class
class SearchState:
    def __init__(self, matrix, lowerBound, depth, path):
        self.lowerBound = lowerBound
        self.matrix = matrix
        self.depth = depth
        self.path = path

    def __lt__(self, other):
        return (self.lowerBound / self.depth) < other.lowerBound / other.depth
