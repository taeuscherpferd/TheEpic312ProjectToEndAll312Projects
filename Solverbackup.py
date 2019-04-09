#!/usr/bin/python3

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
        pass

    ''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

    # So this is all that I have to do for this project... right?
    # According to my understanding this is just the algorithm that uses all of the random matricies everywhere.
    # How to generate/find the matricies: Make a big matrix of all of of the cities.
    # The cost and the lowerbound are definitely different.
    def branchAndBound(self, time_allowance=60.0):
        results = {}
        qCounts = []
        numSolutionsFound = 0
        startTime = time.time()
        self.nodeCount = 0
        self.nodesPruned = 0
        self.cityList = self._scenario.getCities()
        cityCount = len(self.cityList)
        cityMatrix = [[0 for y in range(cityCount)] for x in range(cityCount)]

        #TODO: You'll need this later.... probably....
        bssf = self.defaultRandomTour(time_allowance)['soln']

        for row in range(cityCount):
            for col in range(cityCount):
                if col == row:
                    cityMatrix[row][col] = math.inf
                    continue
                cityMatrix[row][col] = self.cityList[row].costTo(self.cityList[col])

        lowerBoundMatrix, lowerBound = self.calcLowerBoundMatrix(cityMatrix, cityCount)

        startNode = SearchState(lowerBoundMatrix, lowerBound, 0, [self.cityList[0]])

        self.nodeCount += 1
        numSolutionsFound += 1
        heapq.heappush(self.q, startNode)

        #Quit on time out or if queue is empty. 
        elapsedTime = time.time() - startTime 
        while elapsedTime < time_allowance and not len(self.q) == 0:
            qCounts.append(len(self.q))
            curLowestMatrix = heapq.heappop(self.q)
            if (len(curLowestMatrix.path) == cityCount):
                numSolutionsFound += 1
                bssf = TSPSolution(curLowestMatrix.path)

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
        curLocation = curPath[-1]._index
        moddedContenderMatrix = deepcopy(matrix)
        for x in range(size):
            moddedContenderMatrix[curLocation][x] = math.inf

        for col in range(size):
            contenderMatrix = deepcopy(moddedContenderMatrix)
            curCost = matrix[curLocation][col] 
            if curCost == math.inf:
                continue
            
            for x in range(size):
                contenderMatrix[x][col] = math.inf

            contenderMatrix[col][curLocation] = math.inf
            
            contenderMatrix, contLowerbound = self.calcLowerBoundMatrix(contenderMatrix, size)
            contLowerbound += curCost
            contLowerbound += curLowerBound

            newPath = deepcopy(curPath)
            newPath.append(self.cityList[col])

            contNode = SearchState(contenderMatrix, contLowerbound, curDepth + 1, newPath)
            self.nodeCount += 1

            #TODO: Change this to check against the best solution so far
            if contLowerbound <= bssf.cost:
                heapq.heappush(self.q, contNode)
            else:
                self.nodesPruned += 1

    def calcLowerBoundMatrix(self, matrix, size):
        currentLowestVal = 0
        lowestRowVals = []
        lowestColVals = []
        lowerBoundMatrix = deepcopy(matrix)
        lowerBound = 0
        # Row Calculations
        for row in range(size):
            currentLowestVal = lowerBoundMatrix[row][0]
            for col in range(size):
                if lowerBoundMatrix[row][col] < currentLowestVal:
                    currentLowestVal = lowerBoundMatrix[row][col]
            
            if currentLowestVal == math.inf:
                currentLowestVal = 0 

            lowestRowVals.append(currentLowestVal)

            for col in range(size):
                if not lowerBoundMatrix[row][col] == math.inf:
                    lowerBoundMatrix[row][col] = lowerBoundMatrix[row][col] - currentLowestVal

        # Column Calculations
        for row in range(size):
            currentLowestVal = lowerBoundMatrix[0][row]
            for col in range(size):
                if lowerBoundMatrix[col][row] < currentLowestVal:
                    currentLowestVal = lowerBoundMatrix[col][row]

            if currentLowestVal == math.inf:
                currentLowestVal = 0 

            lowestColVals.append(currentLowestVal)

            for col in range(size):
                if not lowerBoundMatrix[col][row] == math.inf:
                    lowerBoundMatrix[col][row] = lowerBoundMatrix[col][row] - currentLowestVal

        lowerBound = sum(lowestRowVals) + sum(lowestColVals)
        return lowerBoundMatrix, lowerBound

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


#TODO: make a class to represent searchstates. 
# override the __lt__(self, other)
# is self less than other. Override this with lowerbound / depth
class SearchState:
    def __init__(self, matrix, lowerBound, depth, path):
        self.lowerBound = lowerBound
        self.matrix = matrix
        self.depth = depth
        self.path = path

    def __lt__(self, other):
        return (self.lowerBound / self.depth) < other.lowerBound / other.depth