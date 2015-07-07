#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from clases import *
import util as ut
import numeric as num
LOGp = ['Report from %s.py'%(__name__)]
WARNINGp = ['WARNINGs from %s.py'%(__name__)]


def neighbours(base,max_neig=2,eps_dist=0.09,ncpus=1):
   """
      Look for neigbours as many times as needed until no new neighbours
    are found. This method is not efficient (checks 1 more neighbour
    than needed), but it is effective.
      latt_vec: Lattice vectors
      base_up: List of basis elements
      nvec: number of repetitions of the latt vecs to get
            all the max_neig neighbors
      eps_dist: distance quanta, anything closer than eps_dist is considered
                to be at the same point
      ncpus: number of cpus used to calculate
   """
   nvec =1
   base_up = base.elems
   latt_vec = base.latt
   repeat = True
   aux = 0
   while repeat:
      vecs_cells = ut.vec_neig(latt_vec,nvec)  # all possible combis of
      dists = np.array([0.])                   # lattice vectors
      for r in vecs_cells:
         # Matrix of distances (nx3)
         A = np.matrix(np.zeros((len(base_up),3),dtype = float))
         for i in range(len(base_up)): #x,y,z in zip(X,Y,Z):
            elem = base_up[i]
            A[i,0] = elem.position[0]
            A[i,1] = elem.position[1]
            A[i,2] = elem.position[2]
         # Matrix of distances displaced by vector r (nx3)
         B = np.matrix(np.zeros((len(base_up),3),dtype = float))
         for i in range(len(base_up)): #x,y,z in zip(X,Y,Z):
            elem = base_up[i]
            B[i,0] = elem.position[0] + r[0]
            B[i,1] = elem.position[1] + r[1]
            B[i,2] = elem.position[2] + r[2]
         dists = np.append(dists, num.dists(A,B))

      #   OJO!! if 2 atoms are eps_dist Angstroms away, they are the same atom
      # considering that the size of an H atom is 0.5 Angstroms... may be enough
      #
      # ordenado contains all the different distances in the problem
      ordenado = ut.sort(dists,eps=eps_dist)
      #print '-------',max_neig
      #print ordenado[0:5]

      bonds = [] # list of tuples (vector,neigh matrix)
      for r in vecs_cells:
         A = np.matrix(np.zeros((len(base_up),3),dtype = float))
         for i in range(len(base_up)): #x,y,z in zip(X,Y,Z):
            elem = base_up[i]
            A[i,0] = elem.position[0]
            A[i,1] = elem.position[1]
            A[i,2] = elem.position[2]
         B = np.matrix(np.zeros((len(base_up),3),dtype = float))
         for i in range(len(base_up)): #x,y,z in zip(X,Y,Z):
            elem = base_up[i]
            B[i,0] = elem.position[0] + r[0]
            B[i,1] = elem.position[1] + r[1]
            B[i,2] = elem.position[2] + r[2]
         neig = num.vecin(A,B,ordenado,eps_dist)
         min_neig = neig.min()
         if min_neig <= max_neig:
           bonds.append((r,neig))
      if aux == len(bonds): # No new neighbours stop looking
        repeat = False
      else: # New neighbours check next cells
        repeat = True
        nvec += 1
      aux = len(bonds)
   nvec = nvec -1 #  we tried with one more cell and checked
                  # that no new neighbour appeared
   # bonds, and number of repetitions needed to find all max_neig
   return bonds, nvec
