#!/usr/bin/python3

import numpy as np

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
  # Skip if the row is all np.infinity.
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
  # Skip if the column is all np.infinity.
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

  # Get the current city np.information because
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
  if node.rMatrix()[ start_city.index(), i ] != np.inf:
    new_node = self.explore( cities, node, start_city.index(), i )
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

  # Make the row all np.infinity
  for col in range( len(rmatrix) ):
  rmatrix[ start_idx ][ col ] = np.inf
  # Make the column all np.infinity
  for row in range( len(rmatrix) ):
  rmatrix[ row ][ end_idx ] = np.inf

  # No more turnning back
  rmatrix[ end_idx ][ start_idx ] = np.inf
  # Reduce the matrix and add on
  # additional lower bound
  lb += self.reduceMatrix( rmatrix )
  return Node( rmatrix, lb, path, cities[ end_idx ] )


class Node:
  def __init__( self, reduced_matrix, lb, path, dest ):
    self.rmatrix = reduced_matrix
    self.lb = lb
    self.pppath = path + [ dest ]  # Add the new city to the path

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
