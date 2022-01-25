import pickle
from class_circuit_eq import *

import os
import sys
modulepath = os.path.expanduser('~/Documents/modelling/6eq/modules')  # os.path.expanduser(path) : return the argument with an initial component of ~ or ~user replaced by that user’s home directory.
sys.path.append(modulepath)
path = os.path.expanduser('~/Documents/modelling/6eq/parameter_space_search') #path of project folder: folder where code, results and parameterfiles are found.
from sympy import *

import pickle
from class_circuit_eq import *
import itertools

import numpy as np
from cellular_automaton import CellularAutomaton, MooreNeighborhood, CAWindow, EdgeRule
import random
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def ca():
   ''' Celluar automata with Python - K. Hong'''
   # 64 Boolean - True(1) : '*'
   #            - False(0): '-'
   # Rule - the status of current cell value is True
   # if only one of the two neighbors at the previous step is True('*')
   # otherwise, the current cell status is False('-')

   # list representing the current status of 64 cells
   ca = [
         0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  0,0,0,0]

   # new Cellular values
   ca_new = ca[:]

   # dictionary maps the cell value to a symbol
   dic = {0:'-', 1:'*'}

   # initial draw - step 0
   print  (''.join( [dic[e] for e in ca_new]))
   # additional 31 steps
   step = 1
   while(step < 32):
      ca_new = []
      # loop through 0 to 63 and store the current cell status in ca_new list
      for i in range(0,64):
         # inside cells - check the neighbor cell state
         if i > 0 and i < 63:
            if ca[i-1] == ca[i+1]:
               ca_new.append(0)
            else:
               ca_new.append(1)

         # left-most cell : check the second cell
         elif(i == 0):
            if ca[1] == 1:
               ca_new.append(1)
            else:
               ca_new.append(0)

         # right-most cell : check the second to the last cell
         elif(i == 63):
            if ca[62] == 1:
               ca_new.append(1)
            else:
               ca_new.append(0)

      # draw current cell state
      print  (''.join( [dic[e] for e in ca_new]))

      # update cell list
      ca = ca_new[:]

      # step count
      step += 1

if __name__ == '__main__':
   ca()
