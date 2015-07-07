#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from imports import *
from vecin import neighbours

class Base(object):
     def __init__(self,sitios,latt,recip,nneig=1):
        self.elems = sitios
        base_atoms = []
        for i in range(len(sitios)):
           E = sitios[i]
           for at in E.atom:
              base_atoms.append(at+str(E.place))
        self.atoms = base_atoms
        base_orbs = []
        for i in range(len(sitios)):
           E = sitios[i]
           for orb in E.orbitals:
              base_orbs.append(orb+str(E.place))
        self.orbitals = base_orbs # list of orbitals. The Hamil matrix should
                                  # have its length
        self.nneig = nneig
        self.latt = latt
        self.recip = recip
        self.bonds,_ = neighbours(self, nneig)
        # Work on progress
        #def show(self):
        #   save_list(self,'base.pkl')
        #   os.system('python see_basis.py')
        def remove(self,i):
           selfx = deepcopy(self)
           aux = []
           for E in self.elems:
              if E.place != i:
                 aux.append(E)
           return Base(aux,selfx.latt,selfx.recip,selfx.nneig)

class geome():
    def __init__(self, latt):
       self.latt = latt
       self.recip = ut.reciprocal(latt)
